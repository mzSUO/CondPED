# Heavy numerical validation for .compute_subset_loss() (v1.0 Stage 1).
#
# NOT part of testthat. Per Runbook 3.7: generate 100 random positive
# definite covariance matrices and non-zero effect vectors, iterate
# over random subsets, and record:
#   max_normal_equation_error
#   max_decomposition_error
#   min_loss / max_loss
#   n_rank_deficient
#   n_unstable
# Plus a near-singular sensitivity block.
#
# Run from the package root:
#   Rscript inst/validation/validate_subset_loss.R

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(CondPED)
}
csl <- CondPED:::.compute_subset_loss

failures <- character()
report <- function(ok, label) {
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) failures <<- c(failures, label)
}

set.seed(20260804)
n_rep <- 100L
m <- 5L
trait_names <- paste0("T", seq_len(m))

records <- data.frame(
  rep = integer(), size_S = integer(),
  normal_equation_error = numeric(),
  decomposition_error = numeric(),
  loss = numeric(),
  status = character(),
  stringsAsFactors = FALSE
)
row <- 0L

for (rep in seq_len(n_rep)) {
  # Random PD covariance: crossprod of a random matrix + ridge.
  X <- matrix(stats::rnorm(m * (m + 2)), m, m + 2)
  Sigma <- tcrossprod(X) + diag(0.2, m)
  dimnames(Sigma) <- list(trait_names, trait_names)
  beta <- stats::rnorm(m)
  names(beta) <- trait_names

  # Random subsets of every size, plus the two boundaries.
  subsets <- list(character())
  for (k in seq_len(m - 1L)) {
    subsets <- c(subsets, list(sample(trait_names, k)))
  }
  subsets <- c(subsets, list(trait_names))

  for (S in subsets) {
    out <- csl(beta, Sigma, S)
    ne_err <- NA_real_
    if (!is.null(out$Gamma)) {
      # index with the NORMALISED sets returned by the function
      Sn <- out$representing_set
      Rn <- out$complement_set
      ne_err <- max(abs(Sigma[Sn, Sn, drop = FALSE] %*% out$Gamma -
                          Sigma[Sn, Rn, drop = FALSE]))
    }
    row <- row + 1L
    records[row, ] <- list(
      rep, length(S), ne_err, out$decomposition_error,
      out$representation_loss, out$status$code
    )
  }
}

n_eval <- nrow(records)
general <- records[!is.na(records$normal_equation_error), ]
max_ne <- max(general$normal_equation_error)
max_de <- max(records$decomposition_error, na.rm = TRUE)
min_loss <- min(records$loss, na.rm = TRUE)
max_loss <- max(records$loss, na.rm = TRUE)
n_rd <- sum(records$status == "rank_deficient")
n_unst <- sum(records$status == "unstable")
n_na <- sum(records$status == "not_applicable")

cat(sprintf("evaluations: %d (%d reps x %d subsets)\n",
            n_eval, n_rep, m + 1L))
cat(sprintf("max_normal_equation_error : %.3e\n", max_ne))
cat(sprintf("max_decomposition_error   : %.3e\n", max_de))
cat(sprintf("loss range                : [%.6f, %.6f]\n", min_loss, max_loss))
cat(sprintf("n_rank_deficient          : %d\n", n_rd))
cat(sprintf("n_unstable                : %d\n", n_unst))
cat(sprintf("n_not_applicable          : %d\n", n_na))

report(max_ne < 1e-8, "normal equations hold across 100 random PD cases")
report(max_de < 1e-8, "quadratic-form decomposition holds (rel tol 1e-8)")
report(min_loss >= 0 && max_loss <= 1, "all losses lie in [0, 1]")
report(n_unst == 0L, "no unstable results on well-conditioned PD inputs")

# ---- near-singular sensitivity ---------------------------------------------
cat("\n== near-singular sensitivity ==\n")
eps_levels <- c(1e-6, 1e-9, 1e-13)
for (eps in eps_levels) {
  Sigma <- diag(c(1, 0.5, eps))
  dimnames(Sigma) <- list(c("A", "B", "C"), c("A", "B", "C"))
  beta <- c(A = 1, B = 0.5, C = 0.2)
  out <- csl(beta, Sigma, c("A", "B"))
  cat(sprintf(paste0(
    "eps = %6.0e : status = %-14s ok = %-5s used_pinv = %-5s ",
    "rank_Omega = %d loss = %s raw_loss = %.6g decomp_err = %.2e\n"),
    eps, out$status$code, out$status$ok, out$used_pseudoinverse,
    out$rank_Omega, format(out$representation_loss),
    out$raw_representation_loss, out$decomposition_error))
}
# well-conditioned but high-loss case stays usable (eps = 1e-6)
Sigma <- diag(c(1, 0.5, 1e-6))
dimnames(Sigma) <- list(c("A", "B", "C"), c("A", "B", "C"))
out_ok <- csl(c(A = 1, B = 0.5, C = 0.2), Sigma, c("A", "B"))
report(out_ok$status$ok && out_ok$status$code == "ok" &&
         !is.na(out_ok$representation_loss),
       "well-conditioned near-boundary case remains usable (eps = 1e-6)")
# at eps = 1e-9 / 1e-13 the covariance violates Sigma > 0 numerically:
# rank_deficient, ok = FALSE, formal loss NA, raw diagnostics kept
for (eps in c(1e-9, 1e-13)) {
  Sigma <- diag(c(1, 0.5, eps))
  dimnames(Sigma) <- list(c("A", "B", "C"), c("A", "B", "C"))
  out <- csl(c(A = 1, B = 0.5, C = 0.2), Sigma, c("A", "B"))
  report(out$status$code == "rank_deficient" &&
           !out$status$ok &&
           is.na(out$representation_loss) &&
           is.finite(out$raw_representation_loss),
         sprintf("eps = %.0e: rank_deficient, ok = FALSE, loss = NA, raw kept", eps))
}

# ---- monotonicity spot check on random cases --------------------------------
cat("\n== monotonicity (nested chains, 20 random cases) ==\n")
max_viol <- 0
for (rep in seq_len(20)) {
  X <- matrix(stats::rnorm(m * (m + 2)), m, m + 2)
  Sigma <- tcrossprod(X) + diag(0.2, m)
  dimnames(Sigma) <- list(trait_names, trait_names)
  beta <- stats::rnorm(m); names(beta) <- trait_names
  # a random nested chain: empty -> ... -> full set
  chain <- sample(trait_names)
  rho_prev <- 1  # rho(empty)
  for (k in seq_len(m - 1L)) {
    rho_k <- csl(beta, Sigma, chain[seq_len(k)])$representation_loss
    max_viol <- max(max_viol, rho_k - rho_prev)
    rho_prev <- rho_k
  }
  max_viol <- max(max_viol, 0 - rho_prev)  # rho(full) = 0
}
cat(sprintf("max monotonicity violation (nested chains): %.3e\n", max_viol))
report(max_viol < 1e-8, "loss is non-increasing along nested chains")

cat("\n")
if (length(failures) > 0L) {
  cat("VALIDATION FAILURES:\n")
  for (f in failures) cat(" -", f, "\n")
  quit(status = 1L)
} else {
  cat("All .compute_subset_loss validation checks passed.\n")
}
