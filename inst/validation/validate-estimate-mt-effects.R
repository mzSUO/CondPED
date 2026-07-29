# Heavy statistical validation for estimate_mt_effects().
#
# NOT part of testthat. Heavy acceptance per contract section 4.4:
#   - bias, RMSE and 95% coverage of beta_hat at the focal QTL across
#     architectures (single_trait, shared_same, shared_opposite, dense);
#   - rank-deficient / failure proportions;
#   - all replicates (failures included) kept in the results table / RDS.
#
# Run from the package root:
#   Rscript inst/validation/validate-estimate-mt-effects.R

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(CondPED)
}

failures <- character()
report <- function(ok, label) {
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) failures <<- c(failures, label)
}

n_rep <- 30L
n_ind <- 300L
m_tr  <- 3L
p_snp <- 400L
out_dir <- file.path("inst", "validation", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

architectures <- c("single_trait", "shared_same", "shared_opposite", "dense")
all_out <- list()

for (ai in seq_along(architectures)) {
  arch <- architectures[ai]
  cat(sprintf("\n== %s ==\n", arch))
  beta_est <- matrix(NA_real_, n_rep, m_tr)
  beta_se  <- matrix(NA_real_, n_rep, m_tr)
  # True beta varies across replicates (per-unit effect scales with the
  # sampled QTL's MAF to hit the target PVE), so it must be stored per
  # replicate -- comparing all estimates against a single replicate's truth
  # would fake attenuation and destroy coverage.
  beta_true <- matrix(NA_real_, n_rep, m_tr)
  tab <- data.frame(replicate = seq_len(n_rep), ok = NA_integer_,
                    converged = NA_integer_, status = NA_character_,
                    stringsAsFactors = FALSE)
  for (r in seq_len(n_rep)) {
    out <- tryCatch({
      sim <- simulate_condped_data(n = n_ind, m = m_tr, p = p_snp,
                                   architecture = arch, n_qtl = 1L,
                                   locus_pve = 0.03, seed = 21000 + 1000 * ai + r)
      fit <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 2L,
                         control = list(maxit = 300))
      if (!isTRUE(fit$status$ok)) {
        return(list(ok = FALSE, status = "fit_failed"))
      }
      est <- estimate_mt_effects(fit, sim$G, loci = sim$qtl_index)
      list(ok = TRUE, beta = est$beta[1, ],
           se = est$effects_long$se,
           true = sim$truth$beta[1, ], converged = fit$convergence$code == 0)
    }, error = function(e) list(ok = FALSE, status = conditionMessage(e)))
    tab$ok[r] <- as.integer(isTRUE(out$ok))
    tab$status[r] <- if (isTRUE(out$ok)) "ok" else out$status
    if (isTRUE(out$ok)) {
      beta_est[r, ] <- out$beta
      beta_se[r, ] <- out$se
      beta_true[r, ] <- out$true
      tab$converged[r] <- as.integer(out$converged)
    }
  }
  all_out[[arch]] <- list(beta_est = beta_est, beta_se = beta_se,
                          beta_true = beta_true, tab = tab)
  keep <- which(tab$ok == 1L & !is.na(tab$ok))
  err <- beta_est[keep, , drop = FALSE] - beta_true[keep, , drop = FALSE]
  bias <- colMeans(err)
  rmse <- sqrt(colMeans(err^2))
  covered <- abs(err / beta_se[keep, , drop = FALSE]) <= 1.96
  coverage <- colMeans(covered)
  mean_true <- colMeans(beta_true[keep, , drop = FALSE])
  cat(sprintf("n_ok = %d/%d, convergence = %.0f%%\n",
              length(keep), n_rep,
              100 * mean(tab$converged == 1L, na.rm = TRUE)))
  cat("mean true beta:", paste(sprintf("%7.3f", mean_true)), "\n")
  cat("bias        :", paste(sprintf("%7.3f", bias)), "\n")
  cat("RMSE        :", paste(sprintf("%7.3f", rmse)), "\n")
  cat("95% coverage:", paste(sprintf("%7.3f", coverage)), "\n")
  report(length(keep) / n_rep >= 0.90,
         sprintf("%s: at least 90%% replicates completed", arch))
  # bias should be small relative to the scale of the effects
  report(all(abs(bias) < 0.5 * pmax(abs(mean_true), 0.2)),
         sprintf("%s: |bias| < 50%% of mean |beta| (or 0.1 floor)", arch))
  # coverage within Monte Carlo tolerance of 95% (se ~ 4%)
  report(all(coverage > 0.85),
         sprintf("%s: 95%% CI coverage > 85%%", arch))
}

rds_file <- file.path(out_dir, "estimate-mt-effects-validation.rds")
saveRDS(list(config = list(n_rep = n_rep, n_ind = n_ind, m = m_tr,
                           p = p_snp, locus_pve = 0.03),
             architectures = all_out), rds_file)
cat(sprintf("\nAll replicates (failures included) saved to %s\n", rds_file))

cat("\n")
if (length(failures) > 0L) {
  cat("VALIDATION FAILURES:\n")
  for (f in failures) cat(" -", f, "\n")
  quit(status = 1L)
} else {
  cat("All estimate_mt_effects validation checks passed.\n")
}
