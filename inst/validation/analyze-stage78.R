## Stage 7.8: simulation operating-regime audit.
## No production/statistical/evaluator changes. Candidate-regime table is
## reused from Stage 7.7 P1 (identical seeds/settings); only the two
## missing mechanism diagnostics are computed as small paired pilots:
##   E3: LD-induced pseudo-multitrait (linked_pseudo_multitrait @ r2=0.3,
##       lead-ASSET vs resolved-ASSET vs CondPED), 50 reps
##   R1/R3: conditional representability (highly_representable @ rho=0.05
##       vs strongly_nonredundant), 50 reps each
devtools::load_all(quiet = TRUE)

workers <- as.integer(Sys.getenv("STAGE78_WORKERS", "4"))
reps <- 50L
master_seed <- 20260819L
out_dir <- "inst/validation/output/stage78"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

alpha_omnibus <- 0.05
alpha_signal <- 0.05
alpha_trait <- 0.05
window_bp <- 5000L
tolerance <- 0.10

is_ok <- function(x) {
  is.list(x) && !is.null(x$status) && is.character(x$status) &&
    length(x$status) == 1L && x$status == "ok"
}
plapply <- function(X, FUN, ...) {
  if (.Platform$OS.type == "unix") {
    parallel::mclapply(X, FUN, ..., mc.cores = workers)
  } else {
    lapply(X, FUN, ...)
  }
}

## ==========================================================================
## E3: LD-induced pseudo-multitrait @ r2 = 0.3
## ==========================================================================
e3_rds <- file.path(out_dir, "e3_rows.rds")
if (!file.exists(e3_rds)) {
  e3_rep <- function(rep_id) {
    tryCatch({
      seed <- CondPED:::.seed_for_rep(master_seed, "e3_ldmix_r2.3", rep_id)
      set.seed(seed)
      sim <- simulate_condped_data(
        n = 1000L, m = 4L, p = 1000L,
        experiment = "end_to_end", scenario = "linked_pseudo_multitrait",
        locus_pve = 0.02, secondary_signal_pve = 0.0175, target_r2 = 0.3,
        correlation = "block", tolerance = tolerance
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
      tn <- fit$trait_names
      scan <- scan_mt_omnibus(fit, G)
      om <- scan$omnibus
      tl <- as.data.frame(truth$loci)[1, ]
      members <- marker_ids[position >= tl$start & position <= tl$end]
      lead <- members[which.min(om$p_value[match(members, om$marker_id)])]

      eff1 <- estimate_mt_effects(fit, G, targets = causals[1])
      beta1 <- stats::setNames(eff1$effects_long$beta,
                               eff1$effects_long$trait)[tn]
      X_C2 <- G[, causals[2], drop = FALSE]
      proj2 <- CondPED:::.build_conditional_projection(fit, X_C2)
      cond1 <- CondPED:::.conditional_mt_scan(
        proj2, G[, causals[1], drop = FALSE],
        marker_ids = causals[1], return_effects = TRUE
      )
      contam_cond <- if (!is.null(cond1$effects)) {
        abs(unname(cond1$effects$beta[1, "Trait2"]))
      } else {
        NA_real_
      }

      asset_of <- function(beta, se, covm) {
        SZ <- CondPED:::.sigma_z_from_cov(covm)
        if (is.null(SZ)) return(character())
        a <- CondPED:::.run_asset_comparison(
          beta, se = se, Sigma_Z = SZ, trait_names = tn,
          sample_size = nrow(fit$rotation$Y_tilde), backend = NULL)
        a$asset_best_subset
      }

      ## lead-ASSET (marginal lead marker)
      lead_eff <- estimate_mt_effects(fit, G, targets = lead)
      lead_b <- stats::setNames(lead_eff$effects_long$beta,
                                lead_eff$effects_long$trait)[tn]
      lead_s <- stats::setNames(lead_eff$effects_long$se,
                                lead_eff$effects_long$trait)[tn]
      lead_subset <- asset_of(lead_b, lead_s, lead_eff$covariance[, , 1L])

      ## oracle-locus resolver -> resolved ASSET
      pos_v <- unname(position[members])
      lid <- sprintf("chr1:%d-%d", min(pos_v), max(pos_v))
      oracle_locus <- list(
        loci = data.frame(
          locus_id = lid, chromosome = "chr1", start = min(pos_v),
          end = max(pos_v), lead_snp = lead,
          lead_p = om$p_value[match(lead, om$marker_id)],
          n_significant_markers = length(members),
          n_region_markers = length(members), status = "ok",
          stringsAsFactors = FALSE),
        membership = data.frame(
          locus_id = lid, marker_id = members,
          significant_in_marginal_scan = FALSE,
          r2_to_lead = as.numeric(stats::cor(G[, members], G[, lead])^2),
          position = pos_v, stringsAsFactors = FALSE)
      )
      res <- resolve_locus_signals(fit, G, oracle_locus,
                                   marker_ids = marker_ids,
                                   signal_adjust = "within_locus_bonferroni",
                                   alpha_signal = alpha_signal)
      res_signals <- res$signals
      res_breadth <- numeric(0)
      res_pseudo <- logical(0)
      res_breadth_err_v <- numeric(0)
      for (k in seq_len(nrow(res_signals))) {
        sid <- res_signals$signal_id[k]
        rw <- res$beta[res$beta$signal_id == sid, ]
        b <- stats::setNames(rw$beta, rw$trait)[tn]
        s <- stats::setNames(rw$se, rw$trait)[tn]
        blk <- (res_signals$signal_order[k] - 1L) * length(tn) +
          seq_along(tn)
        covb <- res$covariance[[res_signals$locus_id[k]]][blk, blk,
                                                          drop = FALSE]
        sub <- asset_of(b, s, covb)
        ## pair resolved signal with truth signal by representative-to-
        ## causal r2 (same rule as Stage 7.3/7.4 ldmix)
        rep_snp <- res_signals$representative_snp[k]
        r2c <- stats::cor(G[, rep_snp], G[, causals])^2
        tset <- truth_sets[[which.max(r2c)]]
        res_breadth <- c(res_breadth, length(sub))
        res_breadth_err_v <- c(res_breadth_err_v,
                               length(sub) - length(tset))
        res_pseudo <- c(res_pseudo,
                        length(sub) > 0 &&
                          (length(sub) > length(tset) ||
                             !all(sub %in% tset)))
      }
      truth_breadth1 <- length(truth_sets[[1]])

      list(
        rep_id = rep_id, status = "ok",
        r2_empirical = as.numeric(stats::cor(G[, causals[1]],
                                             G[, causals[2]])^2),
        contam_marginal = abs(unname(beta1["Trait2"])),
        contam_conditional = contam_cond,
        lead_breadth = length(lead_subset),
        lead_breadth_err = length(lead_subset) - truth_breadth1,
        lead_pseudo = length(lead_subset) > truth_breadth1,
        n_resolved = nrow(res_signals),
        res_breadth_mean = if (length(res_breadth)) mean(res_breadth) else NA,
        res_breadth_err = if (length(res_breadth_err_v)) {
          mean(res_breadth_err_v)
        } else {
          NA
        },
        res_pseudo_rate = if (length(res_pseudo)) mean(res_pseudo) else NA,
        resolved_count_recovered = nrow(res_signals) == 2L
      )
    }, error = function(e) {
      list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
    })
  }
  t0 <- Sys.time()
  e3_rows <- plapply(seq_len(reps), e3_rep)
  ok <- vapply(e3_rows, is_ok, logical(1))
  cat(sprintf("E3: ok %d/%d wall %.0fs\n", sum(ok), length(e3_rows),
              difftime(Sys.time(), t0, units = "secs")))
  saveRDS(e3_rows, e3_rds)
} else {
  e3_rows <- readRDS(e3_rds)
  cat("E3: cached\n")
}
e3_rows <- e3_rows[vapply(e3_rows, is_ok, logical(1))]
e3 <- data.frame(
  n = length(e3_rows),
  mean_empirical_r2 = mean(vapply(e3_rows, function(x) x$r2_empirical,
                                  numeric(1))),
  contam_marginal = mean(vapply(e3_rows, function(x) x$contam_marginal,
                                numeric(1))),
  contam_conditional = mean(vapply(e3_rows, function(x) x$contam_conditional,
                                   numeric(1))),
  contamination_reduction = mean(vapply(e3_rows, function(x) {
    x$contam_marginal - x$contam_conditional
  }, numeric(1))),
  lead_breadth_mean = mean(vapply(e3_rows, function(x) x$lead_breadth,
                                  numeric(1))),
  lead_breadth_err = mean(vapply(e3_rows, function(x) x$lead_breadth_err,
                                 numeric(1))),
  lead_pseudo_rate = mean(vapply(e3_rows, function(x) x$lead_pseudo,
                                 logical(1))),
  resolved_breadth_mean = mean(vapply(e3_rows, function(x) x$res_breadth_mean,
                                      numeric(1)), na.rm = TRUE),
  resolved_breadth_err = mean(vapply(e3_rows, function(x) x$res_breadth_err,
                                     numeric(1)), na.rm = TRUE),
  resolved_pseudo_rate = mean(vapply(e3_rows, function(x) x$res_pseudo_rate,
                                     numeric(1)), na.rm = TRUE),
  Delta_LD_breadth = mean(vapply(e3_rows, function(x) x$lead_breadth_err,
                                 numeric(1))) -
    mean(vapply(e3_rows, function(x) x$res_breadth_err, numeric(1)),
        na.rm = TRUE),
  P_count_recovered = mean(vapply(e3_rows, function(x) {
    x$resolved_count_recovered
  }, logical(1)))
)
write.csv(e3, file.path(out_dir, "stage78_e3_ld_mechanism.csv"),
          row.names = FALSE)
print(e3, row.names = FALSE)

## ==========================================================================
## R1 vs R3: conditional representability
## ==========================================================================
r13_rds <- file.path(out_dir, "r13_rows.rds")
if (!file.exists(r13_rds)) {
  r13_rep <- function(rep_id, arch) {
    tryCatch({
      seed <- CondPED:::.seed_for_rep(master_seed,
                                      sprintf("r13_%s", arch), rep_id)
      set.seed(seed)
      sim <- simulate_condped_data(
        n = 1000L, m = 4L, p = 1000L,
        experiment = "trait_representation", scenario = arch,
        locus_pve = 0.02, target_loss = if (arch == "highly_representable") 0.05 else NULL,
        correlation = "block", tolerance = tolerance
      )
      if (!isTRUE(sim$status$ok)) {
        return(list(rep_id = rep_id, status = "sim_failed"))
      }
      G <- sim$G
      truth <- sim$truth
      marker_ids <- colnames(G)
      position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
      fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
      if (!isTRUE(fit$status$ok)) {
        return(list(rep_id = rep_id, status = "fit_failed"))
      }
      tn <- fit$trait_names
      scan <- scan_mt_omnibus(fit, G)
      om <- scan$omnibus
      valid <- !is.na(om$p_value)
      padj <- rep(NA_real_, nrow(om))
      padj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
      selected <- om$marker_id[valid][padj[valid] <= alpha_omnibus]
      out <- list(rep_id = rep_id, status = "ok", arch = arch,
                  detected = length(selected) > 0L)
      if (length(selected) > 0L) {
        loci_obj <- define_associated_loci(
          scan, G, chromosome = rep("chr1", ncol(G)),
          position = unname(position), selected_markers = selected,
          marker_ids = marker_ids, method = "physical",
          window_bp = window_bp
        )
        if (nrow(loci_obj$loci) > 0L) {
          resolution <- resolve_locus_signals(
            fit, G, loci_obj, marker_ids = marker_ids,
            signal_adjust = "within_locus_bonferroni",
            alpha_signal = alpha_signal)
          pipe <- CondPED:::.pipeline_condped_full(
            fit, resolution, alpha_trait = alpha_trait,
            tolerance = tolerance)
          path <- pipe$subset_analysis$tolerance_path
          reps_df <- pipe$subset_analysis$minimum_representative_sets
          ## truth-side Rep sizes (minimum per signal)
          truth_reps <- truth$minimum_representative_sets
          truth_rep_size <- if (!is.null(truth_reps) &&
                                nrow(truth_reps) > 0L) {
            mean(vapply(split(truth_reps$set_size, truth_reps$signal_id),
                        min, numeric(1)))
          } else {
            NA_real_
          }
          est <- resolution$signals
          if (nrow(est) > 0L) {
            est$position <- unname(position[est$representative_snp])
            ## evaluator for exact Rep-family recovery
            estimates <- list(
              signals = est, candidate_sets = pipe$candidate_sets,
              minimum_representative_sets = reps_df,
              irreducible_modules =
                pipe$subset_analysis$irreducible_modules,
              subset_table = pipe$subset_analysis$subset_table,
              positions = position
            )
            ev <- tryCatch(
              CondPED:::.evaluate_one_replicate(
                truth, estimates,
                list(experiment = "trait_representation", architecture = arch,
                     n = 1000, locus_pve = 0.02, correlation = "block",
                     target_loss = if (arch == "highly_representable") 0.05 else NA,
                     tolerance = tolerance, alpha_omnibus = alpha_omnibus),
                runtime = NA_real_, G = G),
              error = function(e) NULL)
            out$rep_exact <- if (is.null(ev)) NA else ev$rep_exact_family_recovery
            out$rho_map_mae <- if (is.null(ev)) NA else ev$representation_map_mae
            ## per-signal min Rep size (estimated, at tolerance)
            out$min_rep_size <- if (nrow(path) > 0) {
              mean(path$minimum_set_size)
            } else {
              NA
            }
            out$mean_rho <- if (!is.null(pipe$subset_analysis$subset_table) &&
                                nrow(pipe$subset_analysis$subset_table) > 0) {
              st <- pipe$subset_analysis$subset_table
              mean(st$representation_loss[st$feasible_primary], na.rm = TRUE)
            } else {
              NA
            }
            ## resolved-signal ASSET breadth (first signal)
            sid1 <- est$signal_id[1]
            rw <- resolution$beta[resolution$beta$signal_id == sid1, ]
            b <- stats::setNames(rw$beta, rw$trait)[tn]
            s <- stats::setNames(rw$se, rw$trait)[tn]
            blk <- (est$signal_order[1] - 1L) * length(tn) + seq_along(tn)
            covb <- resolution$covariance[[est$locus_id[1]]][blk, blk,
                                                             drop = FALSE]
            SZ <- CondPED:::.sigma_z_from_cov(covb)
            out$asset_breadth <- if (is.null(SZ)) NA_integer_ else {
              length(CondPED:::.run_asset_comparison(
                b, se = s, Sigma_Z = SZ, trait_names = tn,
                sample_size = nrow(fit$rotation$Y_tilde),
                backend = NULL)$asset_best_subset)
            }
            out$truth_rep_size <- mean(truth_rep_size)
          }
        }
      }
      out
    }, error = function(e) {
      list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
    })
  }
  t0 <- Sys.time()
  r13_rows <- list()
  for (arch in c("highly_representable", "strongly_nonredundant")) {
    rows <- plapply(seq_len(reps), r13_rep, arch = arch)
    ok <- vapply(rows, is_ok, logical(1))
    cat(sprintf("R13 %s: ok %d/%d\n", arch, sum(ok), length(rows)))
    r13_rows <- c(r13_rows, rows[ok])
  }
  cat(sprintf("R13 wall %.0fs\n", difftime(Sys.time(), t0, units = "secs")))
  saveRDS(r13_rows, r13_rds)
} else {
  r13_rows <- readRDS(r13_rds)
  cat("R13: cached\n")
}

r13_tab <- do.call(rbind, lapply(
  split(r13_rows, vapply(r13_rows, function(x) x$arch, character(1))),
  function(rows) {
    g <- function(f) vapply(rows, function(x) {
      v <- f(x); if (is.null(v) || length(v) == 0L) NA_real_ else v
    }, numeric(1))
    data.frame(
      arch = rows[[1]]$arch, n = length(rows),
      detected = mean(g(function(x) x$detected)),
      rep_exact = mean(g(function(x) x$rep_exact) >= 1 - 1e-8, na.rm = TRUE),
      rho_map_mae = mean(g(function(x) x$rho_map_mae), na.rm = TRUE),
      min_rep_size = mean(g(function(x) x$min_rep_size), na.rm = TRUE),
      truth_rep_size = mean(g(function(x) x$truth_rep_size), na.rm = TRUE),
      mean_rho = mean(g(function(x) x$mean_rho), na.rm = TRUE),
      asset_breadth = mean(g(function(x) x$asset_breadth), na.rm = TRUE)
    )
  }))
rownames(r13_tab) <- NULL
write.csv(r13_tab, file.path(out_dir, "stage78_r1_r3_representation.csv"),
          row.names = FALSE)
print(r13_tab, row.names = FALSE)

## ==========================================================================
## Candidate regime scoring (reuse Stage 7.7 P1)
## ==========================================================================
p1 <- read.csv("inst/validation/output/stage77/stage77_p1_pve.csv")
cand <- p1[p1$secondary_pve %in% c(0.015, 0.0175, 0.02, 0.025), ]
audit <- data.frame(
  candidate = cand$secondary_pve,
  oracle_causal_recovery = cand$oracle_causal_set_recovery,
  oracle_primary_recovery = cand$oracle_primary_Mregion,
  full_recovery = cand$full_recovery,
  secondary_power = cand$oracle_primary_Mregion,
  attribution = cand$attribution_cand_exact,
  direction = cand$direction_exact,
  extra_rate = cand$extra_signal_rate,
  marker_FDP = cand$marker_fdp,
  LD_contamination_reduction = e3$contamination_reduction,
  LD_breadth_separation = e3$Delta_LD_breadth,
  runtime = cand$mean_runtime,
  stringsAsFactors = FALSE
)
r1 <- r13_tab[r13_tab$arch == "highly_representable", ]
r3 <- r13_tab[r13_tab$arch == "strongly_nonredundant", ]
audit$R1_R3_Rep_separation <- r3$min_rep_size - r1$min_rep_size
audit$R1_R3_rho_separation <- r3$mean_rho - r1$mean_rho
write.csv(audit, file.path(out_dir, "stage78_candidate_regimes.csv"),
          row.names = FALSE)
print(audit, row.names = FALSE)
cat("\nSTAGE78 AUDIT DONE\n")
