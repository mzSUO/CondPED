## Stage 9 shared helpers: batched deterministic runner + Wilson CI
## utilities. Sourced by every run-stage9-*.R / analyze-stage9-*.R script.
## No production code here; orchestration only.

stage9_master_seed <- 20260826L

## Run a scenario grid in batches of `batch` reps with resume, partial
## summaries, one automatic retry, disk manifest and an RSS sampler.
## Mirrors the Stage 8.2/8.3/8.4 runner configuration exactly.
stage9_run_grid <- function(grid, out_dir, reps_total = 200L, batch = 50L,
                            workers = 4L,
                            master_seed = stage9_master_seed) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)

  rss_file <- file.path(out_dir, "peak_rss_kb.txt")
  sampler <- parallel::mcparallel({
    master <- Sys.getpid()
    peak <- 0
    repeat {
      rss <- tryCatch({
        pids <- c(master, as.integer(system2("pgrep", c("-P", master),
                                             stdout = TRUE)))
        vals <- as.integer(system2("ps", c("-o", "rss=", "-p",
                                           paste(pids, collapse = ",")),
                                   stdout = TRUE))
        sum(vals)
      }, error = function(e) 0)
      if (is.finite(rss) && rss > peak) {
        peak <- rss
        writeLines(format(peak, scientific = FALSE), rss_file)
      }
      Sys.sleep(5)
    }
  }, silent = TRUE)
  ## the sampler loops forever; collect its last value WITHOUT waiting and
  ## kill the child (mccollect(wait = TRUE) would hang on it)
  on.exit({
    try(parallel::mccollect(sampler, wait = FALSE), silent = TRUE)
    try(tools::pskill(sampler$pid), silent = TRUE)
  }, add = TRUE)

  t0 <- Sys.time()
  for (gi in seq_len(nrow(grid))) {
    row <- grid[gi, , drop = FALSE]
    scen <- row$scenario
    for (cum in seq(batch, reps_total, by = batch)) {
      run_condped_simulation(
        row, reps = cum, out_dir = out_dir,
        master_seed = master_seed, workers = workers, backend = "parallel",
        resume = TRUE, fail_policy = "save_unstable", progress = FALSE
      )
      key <- CondPED:::.scenario_dir_key(CondPED:::.canonical_scenario_id(
        c(as.list(row), list(master_seed = master_seed))))
      files <- list.files(file.path(out_dir, key), pattern = "^rep_.*\\.rds$",
                          full.names = TRUE)
      st <- vapply(files, function(f) {
        o <- tryCatch(readRDS(f), error = function(e) NULL)
        if (is.null(o)) "corrupt" else if (isTRUE(o$status$ok)) "ok" else "failed"
      }, character(1))
      partial <- data.frame(scenario = scen, reps_done = length(files),
                            n_ok = sum(st == "ok"),
                            n_failed = sum(st == "failed"),
                            batch_target = cum,
                            elapsed_s = difftime(Sys.time(), t0,
                                                 units = "secs"),
                            stringsAsFactors = FALSE)
      saveRDS(partial, file.path(out_dir,
                                 sprintf("partial_%s_%03d.rds", scen, cum)))
      cat(sprintf("PROGRESS %s %d/%d reps (ok %d, failed %d) elapsed %.0fs\n",
                  scen, length(files), reps_total, sum(st == "ok"),
                  sum(st == "failed"), partial$elapsed_s))
    }
  }

  ## one automatic retry of failed replicates
  all_files <- list.files(out_dir, pattern = "^rep_[0-9]+\\.rds$",
                          recursive = TRUE, full.names = TRUE)
  st0 <- vapply(all_files, function(f) {
    o <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(o)) "corrupt" else if (isTRUE(o$status$ok)) "ok" else "failed"
  }, character(1))
  retry_log <- data.frame()
  if (any(st0 != "ok")) {
    bad_files <- all_files[st0 != "ok"]
    retry_log <- data.frame(file = bad_files, first_status = st0[st0 != "ok"],
                            stringsAsFactors = FALSE)
    unlink(bad_files)
    run_condped_simulation(
      grid, reps = reps_total, out_dir = out_dir,
      master_seed = master_seed, workers = workers, backend = "parallel",
      resume = TRUE, fail_policy = "save_unstable", progress = FALSE
    )
    retry_log$retry_status <- vapply(bad_files, function(f) {
      o <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(o)) "failed" else if (isTRUE(o$status$ok)) "ok" else "failed"
    }, character(1))
  }
  write.csv(retry_log, file.path(out_dir, "retry_log.csv"), row.names = FALSE)

  all_files <- list.files(out_dir, pattern = "^rep_[0-9]+\\.rds$",
                          recursive = TRUE, full.names = TRUE)
  final <- data.frame(
    file = all_files,
    status = vapply(all_files, function(f) {
      o <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(o)) "corrupt" else if (isTRUE(o$status$ok)) "ok" else "failed"
    }, character(1)),
    stringsAsFactors = FALSE
  )
  saveRDS(final, file.path(out_dir, "manifest.rds"))
  wall <- difftime(Sys.time(), t0, units = "secs")
  write.csv(data.frame(metric = c("wall_seconds", "n_replicates", "n_ok",
                                  "n_failed", "n_retried",
                                  "n_retry_recovered"),
                       value = c(wall, nrow(final), sum(final$status == "ok"),
                                 sum(final$status != "ok"), nrow(retry_log),
                                 sum(retry_log$retry_status == "ok"))),
            file.path(out_dir, "runtime.csv"), row.names = FALSE)
  cat(sprintf("RUN DONE %s  wall %.0fs  ok %d/%d\n", out_dir, wall,
              sum(final$status == "ok"), nrow(final)))
  invisible(final)
}

## Wilson 95% CI helpers
stage9_wilson <- function(k, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- k / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(center - half, center + half)
}
stage9_prop <- function(x) {
  x <- x[is.finite(x)]
  ci <- stage9_wilson(sum(x), length(x))
  c(estimate = mean(x), ci_lo = ci[1], ci_hi = ci[2], n = length(x))
}
stage9_cont <- function(x) {
  x <- x[is.finite(x)]
  c(estimate = mean(x), ci_lo = mean(x) - 1.96 * stats::sd(x) / sqrt(length(x)),
    ci_hi = mean(x) + 1.96 * stats::sd(x) / sqrt(length(x)), n = length(x))
}

## scenario key helper for grid rows
stage9_key <- function(row, master_seed = stage9_master_seed) {
  CondPED:::.scenario_dir_key(CondPED:::.canonical_scenario_id(
    c(as.list(row), list(master_seed = master_seed))))
}

## Strict matching-based full/secondary recovery (Stage 8.2 convention):
## regenerate the identical dataset per rep (deterministic seed) and match
## truth vs estimated signals with the frozen .match_signals.
stage9_strict_recovery <- function(files, workers = 4L) {
  one <- function(f) {
    o <- readRDS(f)
    if (is.null(o) || !isTRUE(o$status$ok)) return(NULL)
    s <- o$settings
    args <- list(n = s$n, m = s$m, p = s$p, experiment = s$experiment,
                 scenario = s$architecture, locus_pve = s$locus_pve,
                 secondary_signal_pve = s$secondary_signal_pve,
                 correlation = s$correlation, tolerance = s$tolerance)
    for (fld in c("local_ld", "target_r2", "target_loss")) {
      if (!is.null(s[[fld]]) && !is.na(s[[fld]])) args[[fld]] <- s[[fld]]
    }
    set.seed(o$seed)
    sim <- do.call(simulate_condped_data, args)
    G <- sim$G
    truth <- sim$truth
    position <- stats::setNames(seq_len(ncol(G)) * 1000, colnames(G))
    mm <- CondPED:::.match_signals(
      truth, list(signals = o$estimates$signals, positions = position), G = G)
    n_truth <- nrow(truth$signals)
    matched_truth <- unique(mm$matches$truth_signal_id)
    full <- n_truth > 0L && nrow(mm$matches) >= n_truth &&
      all(truth$signals$signal_id %in% matched_truth)
    data.frame(full_recovery = full,
               secondary_recovery = n_truth >= 2L && full,
               n_extra = length(mm$extra), stringsAsFactors = FALSE)
  }
  if (.Platform$OS.type == "unix") {
    out <- parallel::mclapply(files, one, mc.cores = workers)
  } else {
    out <- lapply(files, one)
  }
  out <- out[!vapply(out, is.null, logical(1))]
  do.call(rbind, out)
}
