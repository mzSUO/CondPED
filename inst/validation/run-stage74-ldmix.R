## Stage 7.4 E+F: controlled mechanistic benchmark for
## linked_pseudo_multitrait after the same-locus causal geometry fix
## (Stage 7.4 A), now with the REAL ASSET backend (Stage 7.4 F; Stage 7.3
## ran without ASSET installed, so all ASSET subset metrics were void).
## Truth: signal1 -> Trait1 only, signal2 -> Trait2 only, two LD-linked
## signals inside ONE truth locus. The same dataset feeds:
##   * marginal lead + ASSET
##   * conditional stepwise resolver + ASSET
##   * oracle-primary conditional contamination check
##   * COJO-like joint GLS (truth causals only)
## No CondPED statistic is modified.
stopifnot(requireNamespace("ASSET", quietly = TRUE))
devtools::load_all(quiet = TRUE)

master_seed <- 20260815L
reps <- 100L
workers <- as.integer(Sys.getenv("STAGE74_WORKERS", "16"))
secondary_signal_pve <- as.numeric(Sys.getenv("STAGE74_SPVE", "0.0075"))
out_dir <- "inst/validation/output/stage74/ldmix"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenario <- "linked_pseudo_multitrait"
experiment <- "end_to_end"
locus_pve <- 0.02
alpha_signal <- 0.05

one_rep <- function(rep_id) {
  tryCatch({
    psid <- sprintf("%s|lp%.3f_sp%.4f", scenario, locus_pve,
                    secondary_signal_pve)
    seed <- CondPED:::.seed_for_rep(master_seed, psid, rep_id)
    set.seed(seed)

  sim <- simulate_condped_data(
    n = 1000L, m = 4L, p = 1000L,
    experiment = experiment, scenario = scenario,
    locus_pve = locus_pve,
    secondary_signal_pve = secondary_signal_pve,
    correlation = "block", tolerance = 0.10
  )
  if (!isTRUE(sim$status$ok)) {
    return(list(rep_id = rep_id, status = "sim_failed"))
  }
  G <- sim$G
  truth <- sim$truth
  marker_ids <- colnames(G)
  position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
  causals <- truth$causal_markers
  truth_sets <- truth$candidate_traits
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
  if (!isTRUE(fit$status$ok)) {
    return(list(rep_id = rep_id, status = "fit_failed"))
  }
  trait_names <- fit$trait_names
  scan <- scan_mt_omnibus(fit, G)
  om <- scan$omnibus

  ## truth locus membership
  tl <- as.data.frame(truth$loci)[1, ]
  members <- marker_ids[position >= tl$start & position <= tl$end]
  lead <- members[which.min(om$p_value[match(members, om$marker_id)])]

  ## empirical LD between causals
  r_causal <- stats::cor(G[, causals[1]], G[, causals[2]])
  r2_causal <- r_causal^2

  ## marginal effects of each causal on all traits
  eff1 <- estimate_mt_effects(fit, G, targets = causals[1])
  eff2 <- estimate_mt_effects(fit, G, targets = causals[2])
  beta1 <- stats::setNames(eff1$effects_long$beta,
                           eff1$effects_long$trait)[trait_names]
  beta2 <- stats::setNames(eff2$effects_long$beta,
                           eff2$effects_long$trait)[trait_names]

  ## marginal cross-trait contamination:
  ## SNP1 (true Trait1) effect on Trait2; SNP2 (true Trait2) on Trait1
  contam_marginal_snp1 <- abs(beta1["Trait2"])
  contam_marginal_snp2 <- abs(beta2["Trait1"])

  ## lead ASSET (marginal lead effects through adapter)
  lead_eff <- estimate_mt_effects(fit, G, targets = lead)
  lead_beta <- stats::setNames(lead_eff$effects_long$beta,
                               lead_eff$effects_long$trait)[trait_names]
  lead_se <- stats::setNames(lead_eff$effects_long$se,
                             lead_eff$effects_long$trait)[trait_names]
  Sigma_Z_lead <- CondPED:::.sigma_z_from_cov(lead_eff$covariance[, , 1L])
  lead_asset <- if (is.null(Sigma_Z_lead)) {
    CondPED:::.asset_empty_result("rank_deficient")
  } else {
    CondPED:::.run_asset_comparison(
      lead_beta, se = lead_se, Sigma_Z = Sigma_Z_lead,
      trait_names = trait_names,
      sample_size = nrow(fit$rotation$Y_tilde), backend = NULL
    )
  }
  lead_subset <- if (is.null(lead_asset$asset_best_subset)) character() else
    lead_asset$asset_best_subset

  ## valid_mixing: lead subset contains Trait1 AND Trait2, and the extra
  ## trait must come from the other causal signal (not Trait3/Trait4 noise)
  lead_has_T1 <- "Trait1" %in% lead_subset
  lead_has_T2 <- "Trait2" %in% lead_subset
  lead_has_T3 <- "Trait3" %in% lead_subset
  lead_has_T4 <- "Trait4" %in% lead_subset
  valid_mixing <- lead_has_T1 && lead_has_T2 && !lead_has_T3 && !lead_has_T4

  ## oracle-locus pipeline (same resolver as full)
  oracle_locus <- list(
    loci = data.frame(
      locus_id = tl$locus_id, chromosome = tl$chromosome,
      start = tl$start, end = tl$end, lead_snp = lead,
      lead_p = om$p_value[match(lead, om$marker_id)],
      n_significant_markers = length(members),
      n_region_markers = length(members), status = "ok",
      stringsAsFactors = FALSE
    ),
    membership = data.frame(
      locus_id = tl$locus_id, marker_id = members,
      significant_in_marginal_scan = NA,
      r2_to_lead = as.numeric(stats::cor(G[, members], G[, lead])^2),
      position = unname(position[members]),
      stringsAsFactors = FALSE
    )
  )
  res <- resolve_locus_signals(fit, G, oracle_locus, marker_ids = marker_ids,
                               signal_adjust = "within_locus_bonferroni",
                               alpha_signal = alpha_signal)

  ## resolved ASSET subsets
  res_signals <- res$signals
  res_subsets <- list()
  res_betas <- list()
  for (k in seq_len(nrow(res_signals))) {
    sid <- res_signals$signal_id[k]
    rows <- res$beta[res$beta$signal_id == sid, , drop = FALSE]
    b <- stats::setNames(rows$beta, rows$trait)[trait_names]
    s <- stats::setNames(rows$se, rows$trait)[trait_names]
    blk_idx <- (res_signals$signal_order[k] - 1L) * length(trait_names) +
      seq_along(trait_names)
    cov_block <- res$covariance[[res_signals$locus_id[k]]][blk_idx, blk_idx,
                                                           drop = FALSE]
    Sigma_Z <- CondPED:::.sigma_z_from_cov(cov_block)
    adapter <- if (is.null(Sigma_Z)) {
      CondPED:::.asset_empty_result("rank_deficient")
    } else {
      CondPED:::.run_asset_comparison(
        b, se = s, Sigma_Z = Sigma_Z, trait_names = trait_names,
        sample_size = nrow(fit$rotation$Y_tilde), backend = NULL
      )
    }
    res_subsets[[sid]] <- if (is.null(adapter$asset_best_subset)) character() else
      adapter$asset_best_subset
    res_betas[[sid]] <- b
  }

  ## match resolved signals to truth signals by representative-to-causal r2
  res_matched_truth <- character(length(res_subsets))
  res_r2_causal <- numeric(length(res_subsets))
  for (k in seq_along(res_subsets)) {
    rep_snp <- res_signals$representative_snp[k]
    r2c <- stats::cor(G[, rep_snp], G[, causals])^2
    best <- which.max(r2c)
    res_matched_truth[k] <- names(truth_sets)[best]
    res_r2_causal[k] <- as.numeric(r2c[best])
  }

  ## resolved breadth errors
  lead_breadth <- length(lead_subset)
  lead_breadth_err <- lead_breadth - length(truth_sets[[1]])
  res_breadth_errs <- vapply(seq_along(res_subsets), function(k) {
    length(res_subsets[[k]]) - length(truth_sets[[res_matched_truth[k]]])
  }, numeric(1))

  ## oracle-primary conditional contamination:
  ## condition on true SNP1 -> does SNP2's Trait1 contamination shrink?
  X_C1 <- G[, match(causals[1], marker_ids), drop = FALSE]
  proj1 <- CondPED:::.build_conditional_projection(fit, X_C1)
  cond2 <- CondPED:::.conditional_mt_scan(
    proj1, G[, match(causals[2], marker_ids), drop = FALSE],
    marker_ids = causals[2], return_effects = TRUE
  )
  beta2_cond <- if (!is.null(cond2$effects)) {
    stats::setNames(cond2$effects$beta[1, ], trait_names)
  } else {
    stats::setNames(rep(NA_real_, length(trait_names)), trait_names)
  }

  X_C2 <- G[, match(causals[2], marker_ids), drop = FALSE]
  proj2 <- CondPED:::.build_conditional_projection(fit, X_C2)
  cond1 <- CondPED:::.conditional_mt_scan(
    proj2, G[, match(causals[1], marker_ids), drop = FALSE],
    marker_ids = causals[1], return_effects = TRUE
  )
  beta1_cond <- if (!is.null(cond1$effects)) {
    stats::setNames(cond1$effects$beta[1, ], trait_names)
  } else {
    stats::setNames(rep(NA_real_, length(trait_names)), trait_names)
  }
  contam_cond_snp1 <- abs(beta1_cond["Trait2"])   # SNP1 cond on SNP2
  contam_cond_snp2 <- abs(beta2_cond["Trait1"])   # SNP2 cond on SNP1

  ## COJO-like joint GLS: truth causals only
  joint <- CondPED:::.fit_joint_signal_effects(fit, G, causals,
                                               marker_ids = marker_ids)
  joint_beta1 <- stats::setNames(joint$beta[1, ], trait_names)
  joint_beta2 <- stats::setNames(joint$beta[2, ], trait_names)
  contam_joint_snp1 <- abs(joint_beta1["Trait2"])
  contam_joint_snp2 <- abs(joint_beta2["Trait1"])

  ## truth beta for causals (from truth$beta)
  truth_beta <- truth$beta
  tb1 <- stats::setNames(truth_beta$beta[truth_beta$signal_id ==
                                           names(truth_sets)[1]],
                         truth_beta$trait[truth_beta$signal_id ==
                                            names(truth_sets)[1]])[trait_names]
  tb2 <- stats::setNames(truth_beta$beta[truth_beta$signal_id ==
                                           names(truth_sets)[2]],
                         truth_beta$trait[truth_beta$signal_id ==
                                            names(truth_sets)[2]])[trait_names]

  list(
    rep_id = rep_id, status = "ok",
    r_causal = r_causal, r2_causal = r2_causal,
    lead_snp = lead, lead_p = om$p_value[match(lead, om$marker_id)],
    lead_subset = lead_subset,
    lead_has_T1 = lead_has_T1, lead_has_T2 = lead_has_T2,
    lead_has_T3 = lead_has_T3, lead_has_T4 = lead_has_T4,
    valid_mixing = valid_mixing,
    lead_breadth = lead_breadth, lead_breadth_err = lead_breadth_err,
    n_resolved = nrow(res_signals),
    signal_count_recovered = nrow(res_signals) == length(truth_sets),
    res_subsets = res_subsets,
    res_matched_truth = res_matched_truth,
    res_r2_causal = res_r2_causal,
    res_breadth_errs = res_breadth_errs,
    marginal_beta1 = beta1, marginal_beta2 = beta2,
    cond_beta1 = beta1_cond, cond_beta2 = beta2_cond,
    joint_beta1 = joint_beta1, joint_beta2 = joint_beta2,
    truth_beta1 = tb1, truth_beta2 = tb2,
    contam_marginal_snp1 = contam_marginal_snp1,
    contam_marginal_snp2 = contam_marginal_snp2,
    contam_cond_snp1 = contam_cond_snp1,
    contam_cond_snp2 = contam_cond_snp2,
    contam_joint_snp1 = contam_joint_snp1,
    contam_joint_snp2 = contam_joint_snp2,
    paired_improvement = abs(lead_breadth_err) -
      mean(abs(res_breadth_errs))
    )
  }, error = function(e) {
    list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
  })
}

t0 <- Sys.time()
if (.Platform$OS.type == "unix") {
  out <- parallel::mclapply(seq_len(reps), one_rep, mc.cores = workers)
} else {
  out <- lapply(seq_len(reps), one_rep)
}
wall <- difftime(Sys.time(), t0, units = "secs")
cat("WALL_SECONDS:", round(wall, 1), "\n")

is_ok <- function(x) {
  is.list(x) && !is.null(x$status) && is.character(x$status) &&
    length(x$status) == 1L && x$status == "ok"
}
ok <- vapply(out, is_ok, logical(1))
cat("OK replicates:", sum(ok), "/", length(out), "\n")
saveRDS(out, file.path(out_dir, "ldmix_rows.rds"))

## ---------- summary -------------------------------------------------------
rows <- out[ok]
getb <- function(f) vapply(rows, f, logical(1))
getn <- function(f) vapply(rows, f, numeric(1))

summ <- data.frame(
  metric = c(
    "P_valid_mixing",
    "median_empirical_r2",
    "median_abs_r",
    "mean_marginal_contam_snp1_Trait2",
    "mean_marginal_contam_snp2_Trait1",
    "mean_cond_contam_snp1_Trait2",
    "mean_cond_contam_snp2_Trait1",
    "mean_joint_contam_snp1_Trait2",
    "mean_joint_contam_snp2_Trait1",
    "contam_reduction_marginal_to_joint",
    "lead_breadth_error",
    "resolved_breadth_error",
    "P_resolution_improves_breadth",
    "P_signal_count_recovered",
    "median_rep_to_causal_r2"
  ),
  value = c(
    mean(getb(function(x) x$valid_mixing)),
    stats::median(getn(function(x) x$r2_causal)),
    stats::median(abs(getn(function(x) x$r_causal))),
    mean(getn(function(x) x$contam_marginal_snp1)),
    mean(getn(function(x) x$contam_marginal_snp2)),
    mean(getn(function(x) x$contam_cond_snp1)),
    mean(getn(function(x) x$contam_cond_snp2)),
    mean(getn(function(x) x$contam_joint_snp1)),
    mean(getn(function(x) x$contam_joint_snp2)),
    mean(getn(function(x) x$contam_marginal_snp1) -
           getn(function(x) x$contam_joint_snp1)),
    mean(getn(function(x) x$lead_breadth_err)),
    mean(unlist(lapply(rows, function(x) x$res_breadth_errs))),
    mean(getn(function(x) x$paired_improvement) > 0),
    mean(getb(function(x) x$signal_count_recovered)),
    stats::median(unlist(lapply(rows, function(x) x$res_r2_causal)))
  )
)
print(summ, row.names = FALSE)
write.csv(summ, file.path(out_dir, "ldmix_summary.csv"), row.names = FALSE)
cat("LDMIX CONTROLLED BENCHMARK DONE\n")
