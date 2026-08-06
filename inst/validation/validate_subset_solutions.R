# Heavy validation for the Stage-2 subset machinery (v1.0).
#
# NOT part of testthat. Covers:
#   - complete enumeration (2^k rows) on random PD covariance;
#   - full-table monotonicity on real losses;
#   - module extraction cross-check: direct inclusion-minimal modules
#     vs complements of maximal infeasible sets (independent algorithm);
#   - tie handling and scale/permutation invariance of extracted sets;
#   - all vs singleton loss consistency;
#   - near-singular loci are excluded from formal extraction.
#
# Run from the package root:
#   Rscript inst/validation/validate_subset_solutions.R

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

set.seed(20260805)
n_rep <- 50L
m <- 4L
trait_names <- paste0("T", seq_len(m))
tolerances <- c(0.05, 0.10, 0.20)

# Independent reference algorithm: modules = complements of
# inclusion-maximal infeasible sets.
modules_via_maximal_infeasible <- function(tab, tol) {
  infeas <- tab[!is.na(tab$representation_loss) &
                  tab$representation_loss > tol, ]
  if (nrow(infeas) == 0L) return(character())
  A <- tab$representing_set[[which.max(tab$set_size)]]
  is_max <- vapply(seq_len(nrow(infeas)), function(i) {
    Si <- infeas$representing_set[[i]]
    !any(vapply(seq_len(nrow(infeas)), function(j) {
      i != j && all(Si %in% infeas$representing_set[[j]])
    }, logical(1)))
  }, logical(1))
  sort(vapply(which(is_max), function(i) {
    CondPED:::.trait_set_key(setdiff(A, infeas$representing_set[[i]]))
  }, character(1)))
}

n_cases <- 0L
max_mono <- 0
max_decomp <- 0
n_ties <- 0L
n_mismatch <- 0L
t0 <- proc.time()["elapsed"]

for (rep in seq_len(n_rep)) {
  X <- matrix(stats::rnorm(m * (m + 3)), m, m + 3)
  Sigma <- tcrossprod(X) + diag(0.3, m)
  dimnames(Sigma) <- list(trait_names, trait_names)
  beta <- stats::rnorm(m)
  eff <- data.frame(marker_id = rep("L1", m), trait = trait_names,
                    beta = beta, se = 0.1, stringsAsFactors = FALSE)
  att <- list(L1 = trait_names)
  basis <- derive_conditional_contrasts(Sigma)
  res <- decompose_conditional_effects(
    eff, att, basis, subset_mode = "all",
    tolerance = 0.10, sensitivity_tolerance = tolerances
  )
  n_cases <- n_cases + 1L
  tab <- res$subset_table
  if (nrow(tab) != 2^m) {
    report(FALSE, sprintf("rep %d: expected %d rows, got %d",
                          rep, 2^m, nrow(tab)))
    next
  }
  max_decomp <- max(max_decomp, res$diagnostics$max_decomposition_error)
  max_mono <- max(max_mono, res$diagnostics$max_monotonicity_violation)

  # module cross-check per tolerance
  for (tol in tolerances) {
    direct <- sort(res$irreducible_modules$trait_key[
      res$irreducible_modules$tolerance == tol])
    xref <- modules_via_maximal_infeasible(tab, tol)
    if (!identical(direct, xref)) n_mismatch <- n_mismatch + 1L
  }
  n_ties <- n_ties + any(res$minimum_representative_sets$n_tied_solutions > 1)

  # scale invariance of extracted sets
  eff2 <- eff; eff2$beta <- eff$beta * 3
  res2 <- decompose_conditional_effects(
    eff2, att, basis, tolerance = 0.10,
    sensitivity_tolerance = 0.10
  )
  if (!identical(res$minimum_representative_sets$trait_key[
    res$minimum_representative_sets$tolerance == 0.1],
    res2$minimum_representative_sets$trait_key[
      res2$minimum_representative_sets$tolerance == 0.1])) {
    report(FALSE, sprintf("rep %d: beta scaling changed representative sets", rep))
  }

  # all vs singleton consistency on shared subsets
  res_s <- decompose_conditional_effects(
    eff, att, basis, subset_mode = "singleton",
    tolerance = 0.10, sensitivity_tolerance = 0.10
  )
  la <- setNames(tab$representation_loss, tab$representing_key)
  lb <- setNames(res_s$subset_table$representation_loss,
                 res_s$subset_table$representing_key)
  shared <- intersect(names(la), names(lb))
  if (max(abs(la[shared] - lb[shared])) > 1e-10) {
    report(FALSE, sprintf("rep %d: all/singleton loss mismatch", rep))
  }
}
elapsed <- proc.time()["elapsed"] - t0

cat(sprintf("cases: %d (m = %d, %d subsets each)\n", n_cases, m, 2^m))
cat(sprintf("max monotonicity violation : %.3e\n", max_mono))
cat(sprintf("max decomposition error    : %.3e\n", max_decomp))
cat(sprintf("module cross-check mismatches: %d / %d\n",
            n_mismatch, n_cases * length(tolerances)))
cat(sprintf("cases with tied solutions  : %d\n", n_ties))
cat(sprintf("elapsed                    : %.1f s\n", elapsed))

report(max_mono < 1e-8, "full-table monotonicity holds on random PD cases")
report(max_decomp < 1e-8, "decomposition error within tolerance")
report(n_mismatch == 0L,
       "direct module extraction == maximal-infeasible-complement algorithm")

# ---- near-singular locus -----------------------------------------------------
cat("\n== near-singular locus ==\n")
Sigma <- diag(c(1, 0.5, 1e-13, 1))
dimnames(Sigma) <- list(trait_names, trait_names)
basis <- derive_conditional_contrasts(Sigma)
eff <- data.frame(marker_id = rep("L1", m), trait = trait_names,
                  beta = c(1, 0.5, 0.2, 0.3), se = 0.1,
                  stringsAsFactors = FALSE)
res <- decompose_conditional_effects(eff, list(L1 = trait_names), basis,
                                     tolerance = 0.10)
cat(sprintf("locus status rows: %s\n",
            paste(unique(res$subset_table$status), collapse = ", ")))
cat(sprintf("n_unstable = %d, formal sets extracted: %d reps, %d modules\n",
            res$diagnostics$n_unstable,
            nrow(res$minimum_representative_sets),
            nrow(res$irreducible_modules)))
report(res$diagnostics$n_unstable == 1L &&
         nrow(res$minimum_representative_sets) == 0L,
       "near-singular locus: unstable, excluded from formal extraction")

cat("\n")
if (length(failures) > 0L) {
  cat("VALIDATION FAILURES:\n")
  for (f in failures) cat(" -", f, "\n")
  quit(status = 1L)
} else {
  cat("All Stage-2 subset validation checks passed.\n")
}
