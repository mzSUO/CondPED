## Stage 8.4 analysis: Simulation III formal 500 reps.
## E1/E2: representation + attribution metrics from the frozen evaluator,
## pipeline breadth comparison (lead vs resolved vs CondPED), Q-form
## strength matching check, sign-pattern reporting (freeze v2).
## E3: MAIN evidence per freeze 1.27 = conditional effect contamination
## (marginal vs resolved vs oracle-conditional effect error at the true
## causals); trait-breadth inflation descriptive only.
devtools::load_all(quiet = TRUE)

workers <- as.integer(Sys.getenv("STAGE84_WORKERS", "4"))
out_dir <- "inst/formal-simu/output/stage84"
final <- readRDS(file.path(out_dir, "manifest.rds"))
runtime <- read.csv(file.path(out_dir, "runtime.csv"))
retry_log <- tryCatch(read.csv(file.path(out_dir, "retry_log.csv")),
                      error = function(e) data.frame())
peak_kb <- suppressWarnings(as.numeric(readLines(
  file.path(out_dir, "peak_rss_kb.txt"))))
trait_names <- paste0("Trait", 1:4)

arch_of <- function(o) {
  a <- o$settings$architecture
  if (is.null(a) || length(a) == 0L || !nzchar(a)) a <- o$settings$scenario
  a
}
canon_sign <- function(beta) {
  s <- sign(beta)
  nz <- which(s != 0)
  if (s[nz[1]] < 0) s <- -s
  paste(ifelse(s > 0, "+", "-"), collapse = "")
}

## ---- per-rep extraction (all scenarios, no regen) -----------------------------
rep_extract <- function(f) {
  o <- readRDS(f)
  if (is.null(o) || !isTRUE(o$status$ok)) return(NULL)
  ev <- o$evaluation
  if (!is.data.frame(ev)) return(NULL)
  data.frame(
    scenario = arch_of(o), rep_id = o$rep_id,
    locus_detection = ev$locus_detection_rate,
    cand_exact = ev$candidate_exact_recovery,
    direction = ev$direction_recovery,
    eta_bias = ev$conditional_effect_bias,
    eta_rmse = ev$conditional_effect_rmse,
    rho_mae = ev$representation_loss_mae,
    beta_coverage = ev$beta_coverage,
    threshold_acc = ev$representation_threshold_accuracy,
    rep_card = ev$rep_min_cardinality_recovery,
    rep_exact = ev$rep_exact_family_recovery,
    irr_recovery = ev$irr_family_recovery,
    beta_rmse = ev$beta_rmse,
    lead_breadth_infl = ev$lead_trait_breadth_inflation,
    resolved_breadth_infl = ev$resolved_trait_breadth_inflation,
    pseudo_rate = ev$pseudo_multitrait_interpretation_rate,
    asset_jaccard = ev$asset_trait_set_jaccard,
    asset_dir = ev$asset_direction_recovery,
    overall_recovery = ev$overall_recovery,
    secondary_power = ev$secondary_signal_power,
    stringsAsFactors = FALSE
  )
}
files <- final$file[final$status == "ok"]
t0 <- Sys.time()
rows <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(files, rep_extract, mc.cores = workers)
} else {
  lapply(files, rep_extract)
}
rows <- rows[!vapply(rows, is.null, logical(1))]
d <- do.call(rbind, rows)
write.csv(d, file.path(out_dir, "per_rep_metrics.csv"), row.names = FALSE)
cat(sprintf("extracted %d reps in %.0fs\n", nrow(d),
            difftime(Sys.time(), t0, units = "secs")))

## ---- E1/E2 strength matching (Q-form) + sign distribution ---------------------
qform_of <- function(f) {
  o <- readRDS(f)
  b <- o$truth$B_Q[1L, ]
  SPinv <- CondPED:::.safe_inverse(o$truth$Sigma_P_ref)$inverse
  as.numeric(diag(o$truth$Sigma_X)[1] * sum(b * (SPinv %*% b)))
}
pat_of <- function(f) canon_sign(readRDS(f)$truth$B_Q[1L, ])
e1_files <- files[d$scenario == "single_highly_representable"]
e2_files <- files[d$scenario == "single_nonredundant"]
q_e1 <- mean(vapply(e1_files[seq_len(100L)], qform_of, numeric(1)))
q_e2 <- mean(vapply(e2_files[seq_len(100L)], qform_of, numeric(1)))
e1_pat <- table(vapply(e1_files, pat_of, character(1)))
e2_pat <- table(vapply(e2_files, pat_of, character(1)))

## ---- E3 contamination (freeze 1.27 main evidence; regen + refit) -------------
e3_files <- files[d$scenario == "linked_pseudo_multitrait"]

regen <- function(o) {
  s <- o$settings
  args <- list(n = s$n, m = s$m, p = s$p, experiment = s$experiment,
               scenario = s$architecture, locus_pve = s$locus_pve,
               secondary_signal_pve = s$secondary_signal_pve,
               correlation = s$correlation, tolerance = s$tolerance)
  for (fld in c("local_ld", "target_r2", "target_loss")) {
    if (!is.null(s[[fld]]) && !is.na(s[[fld]])) args[[fld]] <- s[[fld]]
  }
  set.seed(o$seed)
  do.call(simulate_condped_data, args)
}

e3_rep <- function(f) {
  o <- readRDS(f)
  if (is.null(o) || !isTRUE(o$status$ok)) return(NULL)
  sim <- regen(o)
  G <- sim$G
  truth <- sim$truth
  causals <- truth$causal_markers
  marker_ids <- colnames(G)
  position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
  tn <- trait_names

  fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
  if (!isTRUE(fit$status$ok)) {
    return(data.frame(rep_id = o$rep_id, ok = FALSE))
  }
  ## marginal effect of causal1 on Trait2 (true effect = 0)
  eff1 <- estimate_mt_effects(fit, G, targets = causals[1])
  b1 <- stats::setNames(eff1$effects_long$beta,
                        eff1$effects_long$trait)[tn]
  contam_marginal <- unname(b1["Trait2"])
  ## oracle-conditional effect of causal1 given causal2
  X_C <- G[, causals[2], drop = FALSE]
  proj <- CondPED:::.build_conditional_projection(fit, X_C)
  cs <- CondPED:::.conditional_mt_scan(
    proj, G[, causals[1], drop = FALSE], marker_ids = causals[1],
    return_effects = TRUE)
  contam_conditional <- if (!is.null(cs$effects)) {
    unname(cs$effects$beta[1, "Trait2"])
  } else {
    NA_real_
  }
  ## resolved-signal effect on Trait2 (from the saved pipeline output)
  est <- o$estimates$signals
  mm <- CondPED:::.match_signals(
    truth, list(signals = est, positions = position), G = G)
  contam_resolved <- NA_real_
  matched <- mm$matches[mm$matches$truth_signal_id ==
                          truth$signals$signal_id[1], ]
  eff <- o$effect_estimates
  if (nrow(matched) > 0L && !is.null(eff)) {
    es <- matched$est_signal_id[1]
    rw <- eff[eff$signal_id == es & eff$trait == "Trait2", ]
    if (nrow(rw) > 0L) contam_resolved <- rw$beta[1]
  }
  data.frame(
    rep_id = o$rep_id, ok = TRUE,
    r2_causal = as.numeric(stats::cor(G[, causals[1]], G[, causals[2]])^2),
    truth_beta_T2 = truth$beta$beta[truth$beta$signal_id ==
                                      truth$signals$signal_id[1] &
                                      truth$beta$trait == "Trait2"],
    contam_marginal = contam_marginal,
    contam_resolved = contam_resolved,
    contam_conditional = contam_conditional,
    n_matched = nrow(mm$matches), n_extra = length(mm$extra),
    stringsAsFactors = FALSE
  )
}
t0 <- Sys.time()
e3 <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(e3_files, e3_rep, mc.cores = workers)
} else {
  lapply(e3_files, e3_rep)
}
e3 <- e3[!vapply(e3, is.null, logical(1))]
e3 <- do.call(rbind, e3)
write.csv(e3, file.path(out_dir, "e3_contamination.csv"), row.names = FALSE)
cat(sprintf("E3 contamination: %d reps in %.0fs\n", nrow(e3),
            difftime(Sys.time(), t0, units = "secs")))

## ---- aggregates ------------------------------------------------------------------
agg <- function(x) {
  data.frame(
    scenario = x$scenario[1], n = nrow(x),
    locus_detection = mean(x$locus_detection, na.rm = TRUE),
    overall_recovery = mean(x$overall_recovery, na.rm = TRUE),
    cand_exact = mean(x$cand_exact, na.rm = TRUE),
    direction = mean(x$direction, na.rm = TRUE),
    eta_bias = mean(x$eta_bias, na.rm = TRUE),
    eta_rmse = mean(x$eta_rmse, na.rm = TRUE),
    rho_mae = mean(x$rho_mae, na.rm = TRUE),
    beta_coverage = mean(x$beta_coverage, na.rm = TRUE),
    threshold_acc = mean(x$threshold_acc, na.rm = TRUE),
    rep_card = mean(x$rep_card, na.rm = TRUE),
    rep_exact = mean(x$rep_exact, na.rm = TRUE),
    irr_recovery = mean(x$irr_recovery, na.rm = TRUE),
    beta_rmse = mean(x$beta_rmse, na.rm = TRUE),
    lead_breadth_infl = mean(x$lead_breadth_infl, na.rm = TRUE),
    resolved_breadth_infl = mean(x$resolved_breadth_infl, na.rm = TRUE),
    pseudo_rate = mean(x$pseudo_rate, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summ <- do.call(rbind, lapply(split(d, d$scenario), agg))
rownames(summ) <- NULL
write.csv(summ, file.path(out_dir, "summary_scenarios.csv"),
          row.names = FALSE)

e3s <- e3[e3$ok, ]
e3_agg <- data.frame(
  n = nrow(e3s),
  mean_r2 = mean(e3s$r2_causal),
  contam_marginal_abs = mean(abs(e3s$contam_marginal)),
  contam_resolved_abs = mean(abs(e3s$contam_resolved), na.rm = TRUE),
  contam_conditional_abs = mean(abs(e3s$contam_conditional), na.rm = TRUE),
  reduction_marg_to_cond = mean(abs(e3s$contam_marginal) -
                                  abs(e3s$contam_conditional), na.rm = TRUE),
  reduction_marg_to_resolved = mean(abs(e3s$contam_marginal) -
                                      abs(e3s$contam_resolved), na.rm = TRUE),
  n_resolved_na = sum(is.na(e3s$contam_resolved))
)
write.csv(e3_agg, file.path(out_dir, "e3_contamination_summary.csv"),
          row.names = FALSE)

## ---- summary.md --------------------------------------------------------------------
n_ok <- sum(final$status == "ok")
wall <- runtime$value[runtime$metric == "wall_seconds"]
peak_mb <- if (length(peak_kb) && is.finite(peak_kb)) peak_kb / 1024 else NA
s1 <- summ[summ$scenario == "single_highly_representable", ]
s2 <- summ[summ$scenario == "single_nonredundant", ]
s3 <- summ[summ$scenario == "linked_pseudo_multitrait", ]
pat_txt <- function(p) paste(sprintf("%s n=%d (%.1f%%)", names(p), p,
                                     100 * p / sum(p)), collapse = "; ")

lines <- c(
  "# Stage 8.4 Summary — Simulation III formal 500 reps",
  "",
  sprintf("- E1 single_highly_representable (++++, rho=0.05) / E2 single_nonredundant / E3 linked_pseudo_multitrait (r2=0.3, spve=0.020)"),
  sprintf("- R = 500 per scenario; %d/%d replicates ok; retries: %d",
          n_ok, nrow(final), nrow(retry_log)),
  sprintf("- wall time %.0f s; peak RSS %.0f MB; workers = 4", wall, peak_mb),
  "",
  "## E1 vs E2: representation architecture (three paired pipelines)",
  "",
  "| scenario | locus det | overall recovery | cand exact | direction | rep card | Rep exact | Irr | rho MAE | beta cov95 |",
  "|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| E1 | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.4f | %.3f |",
          s1$locus_detection, s1$overall_recovery, s1$cand_exact,
          s1$direction, s1$rep_card, s1$rep_exact, s1$irr_recovery,
          s1$rho_mae, s1$beta_coverage),
  sprintf("| E2 | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.4f | %.3f |",
          s2$locus_detection, s2$overall_recovery, s2$cand_exact,
          s2$direction, s2$rep_card, s2$rep_exact, s2$irr_recovery,
          s2$rho_mae, s2$beta_coverage),
  "",
  sprintf("- eta bias E1 %.4f / E2 %.4f; eta RMSE E1 %.4f / E2 %.4f; threshold-side acc %.3f / %.3f",
          s1$eta_bias, s2$eta_bias, s1$eta_rmse, s2$eta_rmse,
          s1$threshold_acc, s2$threshold_acc),
  sprintf("- Q-form strength: E1 %.4f vs E2 %.4f (freeze 1.17 matched); beta RMSE %.4f / %.4f",
          q_e1, q_e2, s1$beta_rmse, s2$beta_rmse),
  sprintf("- sign patterns: E1: %s; E2: %s", pat_txt(e1_pat), pat_txt(e2_pat)),
  sprintf("- pipeline breadth (descriptive): lead infl E1 %.2f / E2 %.2f; resolved infl E1 %.2f / E2 %.2f",
          s1$lead_breadth_infl, s2$lead_breadth_infl,
          s1$resolved_breadth_infl, s2$resolved_breadth_infl),
  "",
  "## E3: conditional effect contamination (freeze 1.27, MAIN evidence)",
  "",
  sprintf("- empirical r2 between causals: mean %.3f", e3_agg$mean_r2),
  sprintf("- |contamination| marginal: %.4f", e3_agg$contam_marginal_abs),
  sprintf("- |contamination| resolved-signal: %.4f", e3_agg$contam_resolved_abs),
  sprintf("- |contamination| oracle-conditional: %.4f", e3_agg$contam_conditional_abs),
  sprintf("- reduction marginal -> conditional: %.4f (%.0f%%)",
          e3_agg$reduction_marg_to_cond,
          100 * e3_agg$reduction_marg_to_cond / e3_agg$contam_marginal_abs),
  sprintf("- reduction marginal -> resolved: %.4f (%.0f%%)",
          e3_agg$reduction_marg_to_resolved,
          100 * e3_agg$reduction_marg_to_resolved / e3_agg$contam_marginal_abs),
  "",
  "Trait-breadth inflation is descriptive only (freeze 1.27 downgrade):",
  sprintf("- lead breadth inflation %.3f; resolved breadth inflation %.3f; pseudo-multitrait rate %.3f",
          s3$lead_breadth_infl, s3$resolved_breadth_infl, s3$pseudo_rate),
  "",
  "Note: causal1 acts on Trait1 only (true Trait2 effect = 0); contamination",
  "= |estimated Trait2 effect of causal1| under each analysis mode."
)
writeLines(lines, file.path(out_dir, "summary.md"))
print(summ[, c("scenario", "n", "cand_exact", "rep_card", "rep_exact",
               "rho_mae", "beta_rmse")], row.names = FALSE)
print(e3_agg, row.names = FALSE)
cat("STAGE84 ANALYSIS DONE\n")
