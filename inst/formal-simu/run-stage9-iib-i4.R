## Stage 9 batch 12: II-B-I4 — strict same-sign (+ + + +) R1 vs R3 under
## Sigma_P = I4 (freeze §1.58 extension; approved after the feasibility
## probe showed 13% acceptance for both architectures under I4).
## Rejection sampling of concordant vectors within each replicate's
## deterministic RNG stream; acceptance rate recorded and must be reported
## (probe value ~0.13; a material deviation must be flagged).
## R = 500 per architecture, workers = 4.
## Launch: bash run_sim.sh stage9-iib-i4 inst/formal-simu/run-stage9-iib-i4.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-iib-i4"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
reps_total <- 500L
batch <- 50L
workers <- 4L
master_seed <- stage9_master_seed
trait_names <- paste0("Trait", 1:4)

i4_row <- function(scenario, target_loss = NA_real_) {
  data.frame(
    experiment = "trait_representation", scenario = scenario,
    analysis_mode = "signal_trait_oracle",
    comparison_pipeline = "condped_full",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = NA_real_,
    correlation = "independent", target_loss = target_loss,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}

canon_sign <- function(beta) {
  s <- sign(beta)
  nz <- which(s != 0)
  if (length(nz) == 0L) return(NA_character_)
  if (s[nz[1]] < 0) s <- -s
  paste(ifelse(s > 0, "+", "-"), collapse = "")
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
    attempt <- 0L
    sim <- NULL
    repeat {
      attempt <- attempt + 1L
      sim <- do.call(simulate_condped_data, sim_args)
      if (isTRUE(sim$status$ok)) {
        pat <- canon_sign(sim$truth$B_Q[1L, ])
        fullset <- all(vapply(sim$truth$candidate_traits, function(s) {
          setequal(s, trait_names)
        }, logical(1)))
        if (identical(pat, "++++") && fullset) break
      }
      if (attempt >= 500L) stop("architecture_nonconforming after 500 draws")
    }

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
      signals = signals_df, beta = effects_obj$beta,
      candidate_sets = attr$candidate_sets,
      minimum_representative_sets = dec$minimum_representative_sets,
      irreducible_modules = dec$irreducible_modules,
      subset_table = dec$subset_table,
      omnibus_summary = NULL, positions = position, asset = NULL
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

t0 <- Sys.time()
attempt_log <- list()
for (spec in list(list("highly_representable", 0.05),
                  list("strongly_nonredundant", NA_real_))) {
  scen <- spec[[1]]
  row <- i4_row(scen, spec[[2]])
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
        parallel::mclapply(ids, iib_rep, row = row, mc.cores = workers)
      } else {
        lapply(ids, iib_rep, row = row)
      }
    }
    files <- list.files(scen_dir, pattern = "^rep_.*\\.rds$", full.names = TRUE)
    st <- vapply(files, function(f) {
      o <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(o)) "corrupt" else if (isTRUE(o$status$ok)) "ok" else "failed"
    }, character(1))
    cat(sprintf("PROGRESS %s %d/%d reps (ok %d, failed %d) elapsed %.0fs\n",
                scen, length(files), reps_total, sum(st == "ok"),
                sum(st == "failed"), difftime(Sys.time(), t0, units = "secs")))
  }
  files <- list.files(scen_dir, pattern = "^rep_.*\\.rds$", full.names = TRUE)
  att <- vapply(files, function(f) {
    o <- readRDS(f)
    if (is.null(o$n_sim_attempts)) NA_real_ else o$n_sim_attempts
  }, numeric(1))
  attempt_log[[scen]] <- data.frame(
    scenario = scen, n = sum(!is.na(att)),
    mean_attempts = mean(att, na.rm = TRUE),
    median_attempts = stats::median(att, na.rm = TRUE),
    empirical_accept_rate = 1 / mean(att, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

attempts <- do.call(rbind, attempt_log)
write.csv(attempts, file.path(out_dir, "acceptance.csv"), row.names = FALSE)

## metrics
rows <- do.call(rbind, lapply(list.files(
  out_dir, pattern = "^rep_.*\\.rds$", recursive = TRUE, full.names = TRUE),
  function(f) {
    o <- readRDS(f)
    if (is.null(o) || !isTRUE(o$status$ok) ||
        !is.data.frame(o$evaluation)) return(NULL)
    ev <- o$evaluation
    data.frame(scenario = o$settings$scenario,
               rep_exact = ev$rep_exact_family_recovery,
               rep_card = ev$rep_min_cardinality_recovery,
               rho_mae = ev$representation_loss_mae,
               threshold_acc = ev$representation_threshold_accuracy,
               irr_recovery = ev$irr_family_recovery,
               cand_exact = ev$candidate_exact_recovery,
               direction = ev$direction_recovery,
               eta_bias = ev$conditional_effect_bias,
               eta_rmse = ev$conditional_effect_rmse,
               attempts = o$n_sim_attempts,
               stringsAsFactors = FALSE)
  }))
summ <- do.call(rbind, lapply(split(rows, rows$scenario), function(d) {
  data.frame(
    scenario = d$scenario[1], n = nrow(d),
    rep_exact = stage9_prop(d$rep_exact)["estimate"],
    rep_card = stage9_prop(d$rep_card)["estimate"],
    rho_mae = stage9_cont(d$rho_mae)["estimate"],
    threshold_acc = stage9_prop(d$threshold_acc)["estimate"],
    irr_recovery = stage9_prop(d$irr_recovery)["estimate"],
    cand_exact = stage9_prop(d$cand_exact)["estimate"],
    direction = stage9_prop(d$direction)["estimate"],
    eta_bias = stage9_cont(d$eta_bias)["estimate"],
    eta_rmse = stage9_cont(d$eta_rmse)["estimate"],
    mean_attempts = mean(d$attempts),
    stringsAsFactors = FALSE)
}))
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(rows, file.path(out_dir, "per_rep.csv"), row.names = FALSE)

wall <- difftime(Sys.time(), t0, units = "secs")
lines <- c(
  "# Stage 9 batch 12: II-B-I4 strict same-sign comparison (freeze §1.58 ext.)",
  "",
  sprintf("- Sigma_P = I4; R1/R3 with strict ++++ rejection, R=500 each, %d/%d ok",
          nrow(rows), reps_total * 2L),
  "",
  "| scenario | n | rep_exact | rep_card | rho MAE | thr acc | Irr | cand exact | direction | eta bias | eta RMSE | mean attempts |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|",
  apply(summ, 1, function(r) sprintf("| %s | %d | %.3f | %.3f | %.4f | %.3f | %.3f | %.3f | %.3f | %.4f | %.4f | %.1f |",
    r[["scenario"]], as.integer(r[["n"]]), as.numeric(r[["rep_exact"]]),
    as.numeric(r[["rep_card"]]), as.numeric(r[["rho_mae"]]),
    as.numeric(r[["threshold_acc"]]), as.numeric(r[["irr_recovery"]]),
    as.numeric(r[["cand_exact"]]), as.numeric(r[["direction"]]),
    as.numeric(r[["eta_bias"]]), as.numeric(r[["eta_rmse"]]),
    as.numeric(r[["mean_attempts"]]))),
  "",
  sprintf("- empirical acceptance rate: %s (probe expected ~0.13; material deviation must be flagged)",
          paste(sprintf("%s %.3f", attempts$scenario,
                        attempts$empirical_accept_rate), collapse = "; ")),
  sprintf("- wall %.0f s", wall)
)
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
print(attempts, row.names = FALSE)
cat("BATCH12 DONE\n")
