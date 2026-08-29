## Stage 8.4: Simulation III formal 500 reps.
## E1 single_highly_representable (concordant ++++; conforming loop like
## Stage 8.3 R1), E2 single_nonredundant (natural feasible cone, standard
## runner), E3 linked_pseudo_multitrait (r2 = 0.3, spve = 0.020).
## Freeze doc 05 模拟20260825.md v2 sections 1.36-1.46; E1/E2 sign rules
## per the 2026-08-27 constraint fix v2 (strength matching + honest sign
## reporting); E3 main evidence = conditional effect contamination
## (freeze 1.27), trait-breadth inflation descriptive only.
##
## Runner config identical to Stage 8.2/8.3: paired three pipelines,
## deterministic seeds (master_seed 20260826), atomic per-rep RDS
## checkpoints, resume, one automatic retry, workers = 4, batches of 50
## with partial RDS + progress.
##
## Launch ONLY via:  bash run_sim.sh stage84 inst/formal-simu/run-stage84.R

devtools::load_all(quiet = TRUE)

master_seed <- 20260826L
reps_total <- 500L
batch <- 50L
workers <- 4L
out_dir <- "inst/formal-simu/output/stage84"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mk <- function(scenario, target_r2 = NA_real_, target_loss = NA_real_) {
  data.frame(
    experiment = "end_to_end", scenario = scenario,
    analysis_mode = "full", comparison_pipeline = "paired",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = target_r2,
    correlation = "block", target_loss = target_loss,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}
grid <- rbind(
  mk("single_nonredundant"),                  # E2 (natural cone, runner)
  mk("linked_pseudo_multitrait", target_r2 = 0.3)  # E3
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

## ---- E1: concordant conforming loop (full + paired orchestration) -----------
## Mirrors the runner's full + paired path exactly (scan -> BH -> loci ->
## resolve -> .run_comparison_pipelines -> attribution/decompose ->
## evaluator), with rejection sampling of concordant ++++ effect vectors
## within each replicate's deterministic RNG stream (n_sim_attempts
## recorded). E1's feasible cone is dominated by ++++ (87% natural
## acceptance); the loop enforces it per freeze v2 §1.42.
trait_names <- paste0("Trait", 1:4)

e1_row <- function() {
  data.frame(
    experiment = "end_to_end", scenario = "single_highly_representable",
    analysis_mode = "full", comparison_pipeline = "paired",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = NA_real_,
    correlation = "block", target_loss = 0.05,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}

e1_rep <- function(rep_id, row) {
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
      correlation = row$correlation, tolerance = row$tolerance,
      target_loss = row$target_loss
    )
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
      if (attempt >= 500L) stop("architecture_nonconforming after 500 draws")
    }

    ## identical orchestration to the runner's full + paired path
    G <- sim$G
    marker_ids <- colnames(G)
    position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
    fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
    if (!isTRUE(fit$status$ok)) stop("fit_mt_null() did not converge.")
    basis <- derive_conditional_contrasts(fit$Sigma_P_ref)

    scan <- scan_mt_omnibus(fit, G)
    om <- scan$omnibus
    valid <- !is.na(om$p_value)
    p_adj <- rep(NA_real_, nrow(om))
    p_adj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
    selected <- om$marker_id[valid][p_adj[valid] <= row$alpha_omnibus]

    loci_obj <- NULL
    resolution <- NULL
    asset_results <- NULL
    empty_long <- data.frame(
      locus_id = character(), signal_id = character(),
      representative_snp = character(), trait = character(),
      beta = numeric(), se = numeric(), p_value = numeric(),
      stringsAsFactors = FALSE
    )
    if (length(selected) > 0L) {
      loci_obj <- define_associated_loci(
        scan, G, chromosome = rep("chr1", ncol(G)),
        position = unname(position), selected_markers = selected,
        marker_ids = marker_ids, method = "physical", window_bp = 5000
      )
      resolution <- resolve_locus_signals(fit, G, loci_obj,
                                          marker_ids = marker_ids)
      pipes <- CondPED:::.run_comparison_pipelines(
        fit, G, loci_obj, scan, trait_names = fit$trait_names,
        tolerance = row$tolerance
      )
      asset_results <- list(lead_asset = pipes$lead_asset,
                            resolved_asset = pipes$resolved_asset)
      condped_out <- pipes$condped_full
      attr <- list(candidate_sets = condped_out$candidate_sets,
                   signal_table = condped_out$signal_table,
                   trait_table = condped_out$trait_table)
      dec <- condped_out$subset_analysis
      effects_obj <- resolution
    } else {
      attr <- attribute_traits(empty_long, candidate_mode = "all_traits")
      dec <- decompose_conditional_effects(empty_long, attr, basis,
                                           tolerance = row$tolerance)
      effects_obj <- NULL
      condped_out <- NULL
    }

    signals_df <- if (!is.null(effects_obj) &&
                      !is.null(effects_obj$signals)) {
      s <- effects_obj$signals
      if (nrow(s) > 0L) {
        s$position <- unname(position[s$representative_snp])
      }
      s
    } else {
      data.frame()
    }
    estimates <- list(
      signals = signals_df,
      beta = if (!is.null(effects_obj)) effects_obj$beta else NULL,
      candidate_sets = attr$candidate_sets,
      minimum_representative_sets = dec$minimum_representative_sets,
      irreducible_modules = dec$irreducible_modules,
      subset_table = dec$subset_table,
      omnibus_summary = scan$omnibus,
      positions = position,
      asset = CondPED:::.asset_estimates_for_eval(asset_results, condped_out)
    )
    eval_settings <- list(
      experiment = row$experiment, architecture = row$scenario,
      n = row$n, locus_pve = row$locus_pve, correlation = row$correlation,
      target_loss = row$target_loss, tolerance = row$tolerance,
      alpha_omnibus = row$alpha_omnibus
    )
    evaluation <- tryCatch(
      CondPED:::.evaluate_one_replicate(sim$truth, estimates, eval_settings,
                                        runtime = NA_real_, G = G),
      error = function(e) e
    )
    if (inherits(evaluation, "error")) {
      evaluation <- list(error = conditionMessage(evaluation))
    }

    list(
      scenario_id = scenario_id, rep_id = rep_id, seed = seed,
      settings = settings, truth = sim$truth,
      omnibus_summary = scan$omnibus,
      loci = if (!is.null(loci_obj)) loci_obj$loci else NULL,
      signals = signals_df,
      effect_estimates = if (!is.null(effects_obj)) effects_obj$beta else NULL,
      candidate_trait_estimates = attr,
      subset_table = dec$subset_table,
      minimum_representative_sets = dec$minimum_representative_sets,
      irreducible_modules = dec$irreducible_modules,
      tolerance_path = dec$tolerance_path,
      asset_results = asset_results,
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

row <- e1_row()
settings <- as.list(row)
settings$master_seed <- master_seed
key <- CondPED:::.scenario_dir_key(CondPED:::.canonical_scenario_id(settings))
scen_dir <- file.path(out_dir, key)
dir.create(scen_dir, recursive = TRUE, showWarnings = FALSE)
done_ids <- as.integer(sub("^rep_|\\.rds$", "", basename(
  list.files(scen_dir, pattern = "^rep_.*\\.rds$"))))
for (cum in seq(batch, reps_total, by = batch)) {
  ids <- setdiff((cum - batch + 1L):cum, done_ids)
  if (length(ids) > 0L) {
    if (.Platform$OS.type == "unix") {
      parallel::mclapply(ids, e1_rep, row = row, mc.cores = workers)
    } else {
      lapply(ids, e1_rep, row = row)
    }
  }
  files <- list.files(scen_dir, pattern = "^rep_.*\\.rds$", full.names = TRUE)
  st <- vapply(files, function(f) {
    o <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(o)) "corrupt" else if (isTRUE(o$status$ok)) "ok" else "failed"
  }, character(1))
  partial <- data.frame(scenario = "single_highly_representable",
                        reps_done = length(files), n_ok = sum(st == "ok"),
                        n_failed = sum(st == "failed"), batch_target = cum,
                        elapsed_s = difftime(Sys.time(), t0, units = "secs"),
                        stringsAsFactors = FALSE)
  saveRDS(partial, file.path(out_dir, sprintf(
    "partial_%s_%03d.rds", "single_highly_representable", cum)))
  cat(sprintf("PROGRESS %s %d/%d reps (ok %d, failed %d) elapsed %.0fs\n",
              "single_highly_representable", length(files), reps_total,
              sum(st == "ok"), sum(st == "failed"), partial$elapsed_s))
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
  e1_key <- CondPED:::.scenario_dir_key(
    CondPED:::.canonical_scenario_id(settings))
  is_e1 <- basename(dirname(bad_files)) == e1_key
  if (any(!is_e1)) {
    unlink(bad_files[!is_e1])
    run_condped_simulation(
      grid, reps = reps_total, out_dir = out_dir,
      master_seed = master_seed, workers = workers, backend = "parallel",
      resume = TRUE, fail_policy = "save_unstable", progress = FALSE
    )
  }
  if (any(is_e1)) {
    unlink(bad_files[is_e1])
    row <- e1_row()
    ids <- as.integer(sub("^rep_|\\.rds$", "", basename(bad_files[is_e1])))
    if (.Platform$OS.type == "unix") {
      parallel::mclapply(ids, e1_rep, row = row, mc.cores = workers)
    } else {
      lapply(ids, e1_rep, row = row)
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
cat("STAGE84 RUN DONE  wall:", round(wall, 1), "s  ok:",
    sum(final$status == "ok"), "/", nrow(final), "\n")
