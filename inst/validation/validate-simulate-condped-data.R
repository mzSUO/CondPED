# Heavy statistical validation for simulate_condped_data().
#
# This script is NOT part of testthat. It runs the slower Monte Carlo
# checks required by the interface contract (sections 2.1 and 4.1):
#   1. empirical vs theoretical phenotypic variances in a large sample;
#   2. eigenvalue diagnostics of K_bg (structured vs unstructured);
#   3. per-architecture ground-truth tables (A, D, eta identities, PSD);
#   4. PSD preservation across repeated random replicates.
#
# Run from the package root:
#   Rscript inst/validation/validate-simulate-condped-data.R

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

## ---------------------------------------------------------------------------
## 1. Large-sample empirical variances vs theoretical values
## ---------------------------------------------------------------------------
cat("\n== 1. Large-sample empirical variance check ==\n")
# n = 4000 keeps the K eigen decomposition tractable; the tolerance is a
# Monte Carlo tolerance, deliberately looser than any single-replicate
# small-sample bound (contract section 4.1, heavy acceptance).
big <- simulate_condped_data(
  n = 4000, m = 4, p = 600, architecture = "dense",
  n_qtl = 2, locus_pve = 0.02, seed = 20260727
)
emp <- diag(stats::cov(big$Y))
theo <- diag(big$truth$Sigma_P_total)
rel_err <- abs(emp - theo) / theo
cat("theoretical diag:", round(theo, 4), "\n")
cat("empirical  diag:", round(emp, 4), "\n")
cat("relative error :", round(rel_err, 4), "\n")
report(all(rel_err < 0.10), "empirical diag(cov(Y)) within 10% of theory")

## ---------------------------------------------------------------------------
## 2. K_bg eigenvalue diagnostics
## ---------------------------------------------------------------------------
cat("\n== 2. K eigenvalue diagnostics ==\n")
sim_s <- simulate_condped_data(n = 500, m = 4, p = 2000, structured = TRUE,
                               n_groups = 8, fst = 0.05, seed = 1)
sim_u <- simulate_condped_data(n = 500, m = 4, p = 2000, structured = FALSE,
                               seed = 1)
cat(sprintf(
  "structured   : sd(eigen K) = %.4f, effective rank = %d, mean(diag K) = %.12f\n",
  sim_s$diagnostics$K_eigen_sd, sim_s$diagnostics$K_rank,
  sim_s$diagnostics$mean_diag_K
))
cat(sprintf(
  "unstructured : sd(eigen K) = %.4f, effective rank = %d, mean(diag K) = %.12f\n",
  sim_u$diagnostics$K_eigen_sd, sim_u$diagnostics$K_rank,
  sim_u$diagnostics$mean_diag_K
))
report(abs(sim_s$diagnostics$mean_diag_K - 1) < 1e-10,
       "structured K: mean(diag K) = 1")
report(abs(sim_u$diagnostics$mean_diag_K - 1) < 1e-10,
       "unstructured K: mean(diag K) = 1")
report(sim_s$diagnostics$K_eigen_sd > 1.5 * sim_u$diagnostics$K_eigen_sd,
       "structured K has clearly larger eigenvalue dispersion")

## ---------------------------------------------------------------------------
## 3. Per-architecture ground-truth tables
## ---------------------------------------------------------------------------
cat("\n== 3. Per-architecture truth tables ==\n")
architectures <- c(
  "null", "single_trait", "shared_same", "shared_opposite", "dense",
  "covariance_aligned", "conditional_deviation", "projection_induced"
)
for (arch in architectures) {
  sim <- simulate_condped_data(
    n = 500, m = 4, p = 1000, architecture = arch,
    n_qtl = 1, delta = 0.5, seed = 42
  )
  tr <- sim$truth
  eta_target <- tr$eta[1, 1]
  beta_target <- tr$beta[1, 1]
  cat(sprintf(
    paste0(
      "%-22s A={%s} D={%s} beta_1=% .4f eta_1=% .4f ",
      "min_eigen(Sigma_G_bg)=% .2e"
    ),
    arch,
    paste(tr$A[[1]], collapse = ","), paste(tr$D[[1]], collapse = ","),
    beta_target, eta_target, sim$diagnostics$min_eigen_Sigma_G_bg
  ))
  ok_identity <- isTRUE(all.equal(
    tr$eta, tr$beta %*% tr$C, tolerance = 1e-10, check.attributes = FALSE
  ))
  ok_psd <- sim$diagnostics$min_eigen_Sigma_G_bg >= -1e-10
  ok_arch <- switch(
    arch,
    null = all(tr$beta == 0),
    covariance_aligned = abs(eta_target) < 1e-10,
    conditional_deviation = eta_target > 0,
    projection_induced = beta_target == 0 && abs(eta_target) > 1e-3,
    TRUE
  )
  report(ok_identity && ok_psd && ok_arch,
         paste0(arch, ": eta identity + PSD + architecture-specific truth"))
}

## ---------------------------------------------------------------------------
## 4. PSD preservation across repeated replicates (small, NOT 500 reps)
## ---------------------------------------------------------------------------
cat("\n== 4. PSD and scaling consistency across 20 replicates ==\n")
n_bad <- 0L
n_scaled <- 0L
for (rep in seq_len(20)) {
  sim <- simulate_condped_data(
    n = 300, m = 4, p = 800, architecture = "dense",
    n_qtl = 2, locus_pve = sample(c(0.005, 0.01, 0.05, 0.2), 1),
    max_attempts = 5, seed = 1000 + rep
  )
  tr <- sim$truth
  ok <- sim$diagnostics$min_eigen_Sigma_G_bg >= -1e-10 &&
    isTRUE(all.equal(tr$Sigma_G_bg + tr$Sigma_Q, tr$Sigma_G_total,
                     tolerance = 1e-8, check.attributes = FALSE)) &&
    isTRUE(all.equal(tr$eta, tr$beta %*% tr$C,
                     tolerance = 1e-8, check.attributes = FALSE))
  if (!ok) n_bad <- n_bad + 1L
  if (sim$generator$scaled) n_scaled <- n_scaled + 1L
}
cat(sprintf("replicates with consistent truth: %d/20 (scaled: %d)\n",
            20L - n_bad, n_scaled))
report(n_bad == 0L, "all 20 replicates: PSD holds and truth stays in sync")

## ---------------------------------------------------------------------------
cat("\n== Summary ==\n")
if (length(failures) == 0L) {
  cat("All heavy validation checks PASSED.\n")
} else {
  cat("FAILED checks:\n")
  for (f in failures) cat(" -", f, "\n")
  quit(status = 1L)
}
