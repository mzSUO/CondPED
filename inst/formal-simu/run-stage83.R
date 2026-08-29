## Stage 8.3: Simulation II formal 500 reps.
## II-A (signal_oracle): A1 trait_specific, A2 two_trait_concordant,
## A3 two_trait_antagonistic, A4 broad_concordant.
## II-B (signal_trait_oracle): R1 highly_representable (rho_main = 0.05),
## R3 strongly_nonredundant. Freeze doc 05 模拟20260825.md sections
## 1.16 (fairness) / 1.18 (truth engine) / 1.30-1.35.
##
## Runner configuration identical to Stage 8.1/8.2: deterministic
## per-replicate seeds (master_seed 20260826, same seed family as 8.2),
## atomic per-rep RDS checkpoints, resume, one automatic retry of failed
## replicates, workers = 4, batches of 50 with partial RDS + progress.
##
## Launch ONLY via:  bash run_sim.sh stage83 inst/formal-simu/run-stage83.R

devtools::load_all(quiet = TRUE)

master_seed <- 20260826L
reps_total <- 500L
batch <- 50L
workers <- 4L
out_dir <- "inst/formal-simu/output/stage83"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mk <- function(scenario, mode, target_loss = NA_real_) {
  data.frame(
    experiment = "trait_representation", scenario = scenario,
    analysis_mode = mode, comparison_pipeline = "condped_full",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = NA_real_,
    correlation = "block", target_loss = target_loss,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}
grid <- rbind(
  mk("trait_specific", "signal_oracle"),          # A1
  mk("two_trait_concordant", "signal_oracle"),    # A2
  mk("two_trait_antagonistic", "signal_oracle"),  # A3
  mk("broad_concordant", "signal_oracle"),        # A4
  mk("strongly_nonredundant", "signal_trait_oracle")  # R3 (unconstrained)
)
## II-B (R1/R3) runs through the script-level conforming loop below, NOT the
## runner: the frozen registry generates R1/R3 effect vectors with mixed
## signs (direction "mixed"), while freeze 1.33 mandates concordant
## (+,+,+,+ up to allele flip). The runner cannot inject the constraint,
## so II-B replicates rejection-sample the effect vector within the
## replicate's own deterministic RNG stream (first conforming draw wins;
## n_sim_attempts recorded) and then call the SAME production functions in
## the SAME order as the runner's signal_trait_oracle path. II-A scenarios
## are fixed-directions by design and stay on the runner.
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

## ---- II-B R1: conforming regeneration loop (R3 runs on the runner) -------
## Freeze 1.16/1.33 (v2): R1's dominant feasible cone is ++++; the 500
## concordant R1 replicates were generated by rejection within each
## replicate's deterministic RNG stream (n_sim_attempts recorded). R3 needs
## no rejection: its runner-accepted draws already land in its feasible cone.
trait_names <- paste0("Trait", 1:4)

iib_row <- function(scenario, target_loss = NA_real_) {
  data.frame(
    experiment = "trait_representation", scenario = scenario,
    analysis_mode = "signal_trait_oracle", comparison_pipeline = "condped_full",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = NA_real_,
    correlation = "block", target_loss = target_loss,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}

iib_rep <- function(rep_id, row) {
  settings <- as.list(row)
  settings$master_seed <- master_seed
  scenario_id <- CondPED:::.canonical_scenario_id(settings)
  seed <- CondPED:::.seed_for_rep(master_seed, scenario_id, rep_id)
  file <- file.path(out_dir, CondPED:::.scenario_dir_key(scenario_id),
                    sprintf("rep_%04d.rds", rep_id))
  t0 <- proc.time()[["elapsed"]]
  result <- tryCatch({
    set.seed(seed)
    sim_args <- list(
      n = row$n, m = row$m, p = row$p, experiment = row$experiment,
      scenario = row$scenario, locus_pve = row$locus_pve,
      secondary_signal_pve = row$secondary_signal_pve,
      correlation = row$correlation, tolerance = row$tolerance
    )
    if (!is.na(row$target_loss)) sim_args$target_loss <- row$target_loss
    ## rejection until freeze-conforming truth (deterministic in-stream)
    attempt <- 0L
    sim <- NULL
    repeat {
      attempt <- attempt + 1L
      sim <- do.call(simulate_condped_data, sim_args)
      if (isTRUE(sim$status$ok)) {
        bq <- sim$truth$B_Q[1L, ]
        nz <- bq != 0
        concordant <- length(unique(sign(bq[nz]))) == 1L
        fullset <- all(vapply(sim$truth$candidate_traits, function(s) {
          setequal(s, trait_names)
        }, logical(1)))
        if (concordant && fullset) break
      }
      if (attempt >= 500L) {
        stop("architecture_nonconforming after 500 draws")
      }
    }

    ## identical orchestration to the runner's signal_trait_oracle path
    G <- sim$G
    marker_ids <- colnames(G)
    position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
    fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
    if (!isTRUE(fit$status$ok)) stop("fit_mt_null() did not converge.")
    basis <- derive_conditional_contrasts(fit$Sigma_P_ref)
    ps <- CondPED:::.oracle_signals_from_truth(sim$truth)
    effects_obj <- CondPED:::.effects_from_predefined(
      fit, G, ps, marker_ids, sqrt(.Machine$double.eps))
    attr <- attribute_traits(
      effects_obj, candidate_mode = "predefined",
      predefined_sets = sim$truth$candidate_traits)
    dec <- decompose_conditional_effects(
      effects_obj, attr, basis, tolerance = row$tolerance)

    signals_df <- effects_obj$signals
    if (nrow(signals_df) > 0L) {
      signals_df$position <- unname(position[signals_df$representative_snp])
    }
    estimates <- list(
      signals = signals_df,
      beta = effects_obj$beta,
      candidate_sets = attr$candidate_sets,
      minimum_representative_sets = dec$minimum_representative_sets,
      irreducible_modules = dec$irreducible_modules,
      subset_table = dec$subset_table,
      omnibus_summary = NULL,
      positions = position,
      asset = NULL
    )
    eval_settings <- list(
      experiment = row$experiment, architecture = row$scenario,
      n = row$n, locus_pve = row$locus_pve, correlation = row$correlation,
      target_loss = if (is.na(row$target_loss)) NA else row$target_loss,
      tolerance = row$tolerance, alpha_omnibus = row$alpha_omnibus
    )
    evaluation <- CondPED:::.evaluate_one_replicate(
      sim$truth, estimates, eval_settings, runtime = NA_real_, G = G)

    list(
      scenario_id = scenario_id, rep_id = rep_id, seed = seed,
      settings = settings, truth = sim$truth,
      omnibus_summary = NULL, loci = NULL,
      signals = signals_df, effect_estimates = effects_obj$beta,
      candidate_trait_estimates = attr,
      subset_table = dec$subset_table,
      minimum_representative_sets = dec$minimum_representative_sets,
      irreducible_modules = dec$irreducible_modules,
      tolerance_path = dec$tolerance_path,
      asset_results = NULL,
      comparison_pipeline = row$comparison_pipeline,
      Sigma_P_ref = fit$Sigma_P_ref,
      evaluation = evaluation, estimates = estimates,
      status = CondPED:::.new_status(ok = TRUE, code = "ok"),
      run_status = "ok",
      diagnostics = list(n_signals = nrow(signals_df),
                         mode = row$analysis_mode,
                         pipeline = row$comparison_pipeline),
      n_sim_attempts = attempt,
      runtime = NA_real_
    )
  }, error = function(e) e)
  runtime <- proc.time()[["elapsed"]] - t0
  if (inherits(result, "error")) {
    obj <- list(
      scenario_id = scenario_id, rep_id = rep_id, seed = seed,
      settings = settings,
      status = CondPED:::.new_status(ok = FALSE, code = "unstable",
                                     message = conditionMessage(result)),
      run_status = "failed", failure_stage = "runtime",
      error_class = class(result)[1L],
      error_message = conditionMessage(result), runtime = runtime
    )
  } else {
    obj <- result
    obj$runtime <- runtime
  }
  CondPED:::.save_replicate_atomic(obj, file)
  invisible(file)
}

for (spec in list(list("highly_representable", 0.05))) {
  scen <- spec[[1]]
  row <- iib_row(scen, spec[[2]])
  ## replace any earlier non-conforming outputs for these scenarios
  settings <- as.list(row)
  settings$master_seed <- master_seed
  key <- CondPED:::.scenario_dir_key(CondPED:::.canonical_scenario_id(settings))
  scen_dir <- file.path(out_dir, key)
  dir.create(scen_dir, recursive = TRUE, showWarnings = FALSE)
  old_files <- list.files(scen_dir, pattern = "^rep_.*\\.rds$",
                          full.names = TRUE)
  ## delete only first-pass (unconstrained) outputs: conforming replicates
  ## carry n_sim_attempts; wiping those would break resume
  if (length(old_files) > 0L) {
    conforming <- vapply(old_files, function(f) {
      o <- tryCatch(readRDS(f), error = function(e) NULL)
      !is.null(o) && !is.null(o$n_sim_attempts)
    }, logical(1))
    if (any(!conforming)) unlink(old_files[!conforming])
  }
  done_ids <- as.integer(sub("^rep_|\\.rds$", "", basename(
    list.files(scen_dir, pattern = "^rep_.*\\.rds$"))))
  for (cum in seq(batch, reps_total, by = batch)) {
    ids <- setdiff((cum - batch + 1L):cum, done_ids)
    if (length(ids) > 0L) {
      if (.Platform$OS.type == "unix") {
        parallel::mclapply(ids, iib_rep, row = row, mc.cores = workers)
      } else {
        lapply(ids, iib_rep, row = row)
      }
    }
    files <- list.files(scen_dir, pattern = "^rep_.*\\.rds$",
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
  iib_keys <- vapply(
    list(iib_row("highly_representable", 0.05)),
    function(r) {
      s <- as.list(r)
      s$master_seed <- master_seed
      CondPED:::.scenario_dir_key(CondPED:::.canonical_scenario_id(s))
    }, character(1))
  is_iib <- basename(dirname(bad_files)) %in% iib_keys
  ## II-A failures: retry via the runner
  if (any(!is_iib)) {
    unlink(bad_files[!is_iib])
    run_condped_simulation(
      grid, reps = reps_total, out_dir = out_dir,
      master_seed = master_seed, workers = workers, backend = "parallel",
      resume = TRUE, fail_policy = "save_unstable", progress = FALSE
    )
  }
  ## II-B failures: retry via the conforming loop (one extra attempt)
  if (any(is_iib)) {
    for (spec in list(list("highly_representable", 0.05))) {
      row <- iib_row(spec[[1]], spec[[2]])
      settings <- as.list(row)
      settings$master_seed <- master_seed
      key <- CondPED:::.scenario_dir_key(
        CondPED:::.canonical_scenario_id(settings))
      scen_dir <- file.path(out_dir, key)
      have <- as.integer(sub("^rep_|\\.rds$", "", basename(
        list.files(scen_dir, pattern = "^rep_.*\\.rds$"))))
      st_dir <- vapply(list.files(scen_dir, pattern = "^rep_.*\\.rds$",
                                  full.names = TRUE), function(f) {
        o <- tryCatch(readRDS(f), error = function(e) NULL)
        if (is.null(o)) "failed" else if (isTRUE(o$status$ok)) "ok" else "failed"
      }, character(1))
      missing_or_failed <- setdiff(seq_len(reps_total), have[st_dir == "ok"])
      if (length(missing_or_failed) > 0L) {
        unlink(file.path(scen_dir, sprintf("rep_%04d.rds",
                                           missing_or_failed)))
        if (.Platform$OS.type == "unix") {
          parallel::mclapply(missing_or_failed, iib_rep, row = row,
                             mc.cores = workers)
        } else {
          lapply(missing_or_failed, iib_rep, row = row)
        }
      }
    }
  }
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
cat("STAGE83 RUN DONE  wall:", round(wall, 1), "s  ok:",
    sum(final$status == "ok"), "/", nrow(final), "\n")
