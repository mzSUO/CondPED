# Heavy statistical validation for fit_mt_null().
#
# This script is NOT part of testthat and is never run by R CMD check.
# It implements the heavy acceptance of contract section 4.2:
#   - identifiable, structured K (contract section 2.1);
#   - scenarios covering n in {500, 1000} and m in {2, 4};
#   - a pilot of 20 replicates per scenario (the formal 500-replicate run
#     is configured below but NOT executed by default);
#   - element-wise bias, relative bias and RMSE of Sigma_G and Sigma_E;
#   - convergence rate, boundary-solution proportion,
#     identifiable_warning proportion, K diagnostics and run time;
#   - every replicate, converged or not, is kept in the results table and
#     saved to RDS (contract section 1.3: failures are never dropped).
#
# Run from the package root:
#   Rscript inst/validation/validate-fit-mt-null.R

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
## Configuration
## ---------------------------------------------------------------------------
pilot_reps  <- 20L          # pilot size per scenario (contract: 20-50)
formal_reps <- 500L         # formal run size; NOT executed by default
run_formal  <- FALSE
scenarios <- list(
  list(n = 500L,  m = 2L, p = 1000L),
  list(n = 500L,  m = 4L, p = 1000L),
  list(n = 1000L, m = 2L, p = 1000L)
)
n_starts <- 3L
maxit    <- 400L
boundary_tol <- 1e-6        # min eigenvalue below this => boundary solution
out_dir <- file.path("inst", "validation", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

run_scenario <- function(sc, reps, seed0) {
  n <- sc$n; m <- sc$m
  est_G <- array(NA_real_, dim = c(m, m, reps))
  est_E <- array(NA_real_, dim = c(m, m, reps))
  truth_G <- truth_E <- NULL
  res <- data.frame(
    replicate = seq_len(reps),
    converged = NA_integer_,
    ok = NA_integer_,
    status_code = NA_character_,
    boundary = NA_integer_,
    identifiable_warning = NA_integer_,
    K_eigen_sd = NA_real_,
    K_rank = NA_integer_,
    seconds = NA_real_,
    stringsAsFactors = FALSE
  )
  for (r in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch({
      sim <- simulate_condped_data(n = n, m = m, p = sc$p,
                                   architecture = "null",
                                   structured = TRUE, n_groups = 8L,
                                   fst = 0.05, seed = seed0 + r)
      # tryCatch evaluates its expression in this frame: plain <- suffices
      truth_G <- sim$truth$Sigma_G_bg
      truth_E <- sim$truth$Sigma_E
      f <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = n_starts,
                       control = list(maxit = maxit))
      attr(f, "K_eigen_sd") <- sim$diagnostics$K_eigen_sd
      attr(f, "K_rank") <- sim$diagnostics$K_rank
      f
    }, error = function(e) e)
    res$seconds[r] <- proc.time()[["elapsed"]] - t0
    if (inherits(fit, "error")) {
      # failed replicates stay in the table (contract section 1.3)
      res$converged[r] <- 0L
      res$ok[r] <- 0L
      res$status_code[r] <- paste("error:", conditionMessage(fit))
      next
    }
    est_G[, , r] <- unname(fit$Sigma_G)
    est_E[, , r] <- unname(fit$Sigma_E)
    res$converged[r] <- as.integer(fit$convergence$code == 0)
    res$ok[r] <- as.integer(isTRUE(fit$status$ok))
    res$status_code[r] <- fit$status$code
    res$boundary[r] <- as.integer(
      fit$diagnostics$min_eigen_Sigma_G < boundary_tol ||
        fit$diagnostics$min_eigen_Sigma_E < boundary_tol)
    res$identifiable_warning[r] <-
      as.integer(fit$diagnostics$identifiable_warning)
    res$K_eigen_sd[r] <- attr(fit, "K_eigen_sd")
    res$K_rank[r] <- attr(fit, "K_rank")
  }
  keep <- which(res$ok == 1L & !is.na(res$ok))
  list(est_G = est_G[, , keep, drop = FALSE],
       est_E = est_E[, , keep, drop = FALSE],
       truth_G = truth_G, truth_E = truth_E, res = res)
}

element_stats <- function(est_arr, true_mat) {
  # est_arr: m x m x reps; returns element-wise bias / relative bias / RMSE
  err <- sweep(est_arr, c(1L, 2L), true_mat)
  bias <- apply(err, c(1L, 2L), mean)
  rmse <- sqrt(apply(err^2, c(1L, 2L), mean))
  rel <- sweep(err, c(1L, 2L), true_mat, `/`)
  rel_bias <- apply(rel, c(1L, 2L), mean)
  rel_rmse <- sqrt(apply(rel^2, c(1L, 2L), mean))
  list(bias = bias, rmse = rmse, rel_bias = rel_bias, rel_rmse = rel_rmse)
}

all_results <- list()
for (si in seq_along(scenarios)) {
  sc <- scenarios[[si]]
  label <- sprintf("n = %d, m = %d (structured K)", sc$n, sc$m)
  cat(sprintf("\n== Scenario %d: %s, pilot = %d replicates ==\n",
              si, label, pilot_reps))
  out <- run_scenario(sc, pilot_reps, seed0 = 10000 + 1000 * si)
  all_results[[label]] <- out
  res <- out$res
  conv_rate <- mean(res$converged == 1L, na.rm = TRUE)
  cat(sprintf("convergence rate           : %.0f%%\n", 100 * conv_rate))
  cat(sprintf("status ok rate             : %.0f%%\n",
              100 * mean(res$ok == 1L, na.rm = TRUE)))
  cat(sprintf("boundary solution share    : %.0f%%\n",
              100 * mean(res$boundary == 1L, na.rm = TRUE)))
  cat(sprintf("identifiable_warning share : %.0f%%\n",
              100 * mean(res$identifiable_warning == 1L, na.rm = TRUE)))
  cat(sprintf("K eigen sd (mean)          : %.4f\n",
              mean(res$K_eigen_sd, na.rm = TRUE)))
  cat(sprintf("K rank (mean)              : %.1f\n",
              mean(res$K_rank, na.rm = TRUE)))
  cat(sprintf("time per fit (mean, sd)    : %.2f s, %.2f s\n",
              mean(res$seconds), stats::sd(res$seconds)))

  sG <- element_stats(out$est_G, out$truth_G)
  sE <- element_stats(out$est_E, out$truth_E)
  cat("diag relative bias Sigma_G :",
      paste(sprintf("%7.3f", diag(sG$rel_bias))), "\n")
  cat("diag relative RMSE Sigma_G :",
      paste(sprintf("%7.3f", diag(sG$rel_rmse))), "\n")
  cat("diag relative bias Sigma_E :",
      paste(sprintf("%7.3f", diag(sE$rel_bias))), "\n")
  cat("diag relative RMSE Sigma_E :",
      paste(sprintf("%7.3f", diag(sE$rel_rmse))), "\n")
  cat(sprintf("max |element bias|  Sigma_G: %.4f (truth diag ~ %.3f)\n",
              max(abs(sG$bias)), mean(diag(out$truth_G))))
  cat(sprintf("max |element bias|  Sigma_E: %.4f (truth diag ~ %.3f)\n",
              max(abs(sE$bias)), mean(diag(out$truth_E))))

  report(conv_rate >= 0.90,
         sprintf("%s: convergence rate >= 90%%", label))
  report(all(abs(diag(sG$rel_bias)) < 0.10),
         sprintf("%s: |relative bias| of diag(Sigma_G) < 10%%", label))
  report(all(abs(diag(sE$rel_bias)) < 0.10),
         sprintf("%s: |relative bias| of diag(Sigma_E) < 10%%", label))
}

## ---------------------------------------------------------------------------
## Quick add-ons: fixed effects and K = I diagnostics
## ---------------------------------------------------------------------------
cat("\n== Fixed-effect recovery (5 replicates) ==\n")
# Single-replicate GLS errors have SE ~ 0.05-0.10 at n = 500 under a
# structured K, so the meaningful gate is on the bias across replicates,
# not on one noisy replicate.
fe_reps <- 5L
fe_errs <- vector("list", fe_reps)
B_true <- matrix(c(0.5, -0.3, 1.0, 0.2, -0.8, 0.4), nrow = 2, ncol = 3)
for (r in seq_len(fe_reps)) {
  sim <- simulate_condped_data(n = 500, m = 3, p = 1000,
                               architecture = "null", structured = TRUE,
                               seed = 3001 + r)
  W <- cbind(1, stats::rnorm(500))
  fit_f <- fit_mt_null(sim$Y + W %*% B_true, W = W[, 2, drop = FALSE],
                       K = sim$K_bg, add_intercept = TRUE,
                       n_starts = n_starts, control = list(maxit = maxit))
  fe_errs[[r]] <- unname(fit_f$fixed_effects - B_true)
}
fe_arr <- simplify2array(fe_errs)
fe_bias <- apply(fe_arr, c(1L, 2L), mean)
fe_rmse <- sqrt(apply(fe_arr^2, c(1L, 2L), mean))
cat("per-replicate max |error|:",
    paste(round(apply(fe_arr, 3L, function(x) max(abs(x))), 4)), "\n")
cat("max |bias| across replicates:", max(abs(fe_bias)), "\n")
cat("max entry RMSE              :", max(fe_rmse), "\n")
report(max(abs(fe_bias)) < 0.10,
       "fixed-effect bias < 0.10 across 5 replicates (n = 500)")

cat("\n== K = I diagnostics ==\n")
set.seed(4001)
Y_i <- matrix(stats::rnorm(200), 100, 2)
fit_i <- fit_mt_null(Y_i, K = diag(100), n_starts = 2L,
                     control = list(maxit = 200))
report(isTRUE(fit_i$diagnostics$identifiable_warning),
       "K = I triggers identifiable_warning")
report(any(grepl("weakly identifiable", fit_i$status$warnings)),
       "K = I surfaces a status warning")

## ---------------------------------------------------------------------------
## Persist everything (failures included), then summarise
## ---------------------------------------------------------------------------
rds_file <- file.path(out_dir, "fit-mt-null-pilot-results.rds")
saveRDS(list(config = list(pilot_reps = pilot_reps,
                           formal_reps = formal_reps,
                           run_formal = run_formal,
                           scenarios = scenarios,
                           n_starts = n_starts, maxit = maxit,
                           boundary_tol = boundary_tol),
             scenarios = all_results),
        rds_file)
cat(sprintf("\nResults (all replicates, failures included) saved to %s\n",
            rds_file))
cat("Formal 500-replicate run is configured but was NOT executed",
    "(run_formal = FALSE).\n")

cat("\n")
if (length(failures) > 0L) {
  cat("VALIDATION FAILURES:\n")
  for (f in failures) cat(" -", f, "\n")
  quit(status = 1L)
} else {
  cat("All fit_mt_null validation checks passed.\n")
}
