## Stage 8.2: Simulation I formal 500 reps (I-0 null, I-1 single,
## I-2 two_linked_trait_specific; freeze doc 05 模拟20260825.md).
##
## Configuration identical to the Stage 8.1 pilot: paired three-pipeline
## mode (Lead+ASSET / Resolved+ASSET / CondPED share the same
## Y/G/truth/seed within a replicate), deterministic per-replicate seeds,
## atomic per-rep RDS checkpoints, one automatic retry of failed reps,
## workers = 4. Per scenario the runner is called in batches of 50 reps
## (resume semantics); after each batch a partial summary RDS is written
## and progress is printed.
##
## Launch ONLY via:  bash run_sim.sh stage82 inst/formal-simu/run-stage82.R

devtools::load_all(quiet = TRUE)

master_seed <- 20260826L
reps_total <- 500L
batch <- 50L
workers <- 4L
out_dir <- "inst/formal-simu/output/stage82"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grid <- data.frame(
  experiment = "signal_resolution",
  scenario = c("null", "single_multi_trait", "two_linked_trait_specific"),
  analysis_mode = "full",
  comparison_pipeline = "paired",
  n = 1000L, m = 4L, p = 1000L,
  locus_pve = 0.02,
  secondary_signal_pve = 0.02,
  local_ld = NA_character_,
  target_r2 = c(NA_real_, NA_real_, 0.3),
  correlation = "block", target_loss = NA_real_,
  tolerance = 0.10, alpha_omnibus = 0.05,
  stringsAsFactors = FALSE
)
write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)

## ---- peak-RSS sampler -------------------------------------------------------
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
on.exit(parallel::mccollect(sampler), add = TRUE)

## ---- per-scenario batched run ------------------------------------------------
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
    ## partial summary from disk
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
    saveRDS(partial, file.path(
      out_dir, sprintf("partial_%s_%03d.rds", scen, cum)))
    cat(sprintf("PROGRESS %s %d/%d reps (ok %d, failed %d) elapsed %.0fs\n",
                scen, length(files), reps_total, sum(st == "ok"),
                sum(st == "failed"), partial$elapsed_s))
  }
}

## ---- one automatic retry of failed replicates --------------------------------
all_files <- list.files(out_dir, pattern = "^rep_[0-9]+\\.rds$",
                        recursive = TRUE, full.names = TRUE)
st0 <- vapply(all_files, function(f) {
  o <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(o)) "corrupt" else if (isTRUE(o$status$ok)) "ok" else "failed"
}, character(1))
retry_log <- data.frame()
if (any(st0 == "failed")) {
  bad_files <- all_files[st0 == "failed"]
  retry_log <- data.frame(file = bad_files, first_status = "failed",
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

## ---- final manifest from disk -------------------------------------------------
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
cat("STAGE82 RUN DONE  wall:", round(wall, 1), "s  ok:",
    sum(final$status == "ok"), "/", nrow(final), "\n")
