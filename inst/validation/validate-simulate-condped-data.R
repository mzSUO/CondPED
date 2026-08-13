# Stage 6C validation: deterministic scenario checks.
#
# NOT part of testthat. Covers the frozen validation scenarios:
# null / one signal / two independent signals / two linked
# trait-specific signals / R1 / R2 / R3 / linked pseudo-multitrait /
# mixed multi-signal / all three analysis modes / ASSET adapter.
# Population truth and structural checks only -- no formal
# 500-replicate simulation here.
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

n_base <- 500L
m_base <- 4L
p_base <- 400L
tol <- 0.10

mk <- function(experiment, scenario, seed, ...) {
  simulate_condped_data(n = n_base, m = m_base, p = p_base,
                        experiment = experiment, scenario = scenario,
                        local_region_size = 30L, seed = seed, ...)
}

check_common <- function(s, label) {
  q <- unname(s$truth$signal_count)
  list(
    psd = CondPED:::.min_eigen_sym(s$truth$Sigma_G_bg) > -1e-8,
    sigma_q = if (q == 0L) {
      isTRUE(all.equal(unname(s$truth$Sigma_Q),
                       matrix(0, m_base, m_base), tolerance = 1e-10))
    } else {
      isTRUE(all.equal(
        s$truth$Sigma_Q,
        crossprod(s$truth$B_Q, s$truth$Sigma_X %*% s$truth$B_Q),
        tolerance = 1e-10))
    },
    sigma_bg = isTRUE(all.equal(
      s$truth$Sigma_G_bg, s$truth$Sigma_G_total - s$truth$Sigma_Q,
      tolerance = 1e-10)),
    trace = abs(s$diagnostics$trace_K_over_n - 1) < 1e-8,
    mono = q == 0L || CondPED:::.check_loss_monotonicity(
      s$truth$subset_table)$n_violations == 0L
  )
}

cat("== 1. null ==\n")
s <- mk("signal_resolution", "null", 1)
ck <- check_common(s, "null")
report(all(unlist(ck)), "null: covariance/GRM/monotonicity checks")
report(s$truth$signal_count == 0L &&
         nrow(s$truth$subset_table) == 0L,
       "null: zero signals, empty truth tables")

cat("== 2. one signal ==\n")
s <- mk("signal_resolution", "single_multi_trait", 2)
ck <- check_common(s)
report(all(unlist(ck)), "single_multi_trait: common checks")
report(s$truth$signal_count == 1L &&
         identical(unname(s$truth$candidate_traits[[1L]]),
                   c("Trait1", "Trait2")),
       "single_multi_trait: candidate set {T1, T2}")

cat("== 3. two independent signals ==\n")
s <- mk("signal_resolution", "two_heterogeneous", 3,
        local_ld = "none")
ck <- check_common(s)
report(all(unlist(ck)), "two_heterogeneous: common checks")
report(s$truth$signal_count == 2L &&
         s$truth$local_ld$realized_mean_r2 < 0.10,
       "two independent signals: LD ~ 0, both signals present")

cat("== 4. two linked trait-specific signals ==\n")
s <- mk("signal_resolution", "two_linked_trait_specific", 4,
        target_r2 = 0.3)
ck <- check_common(s)
report(all(unlist(ck)), "two_linked_trait_specific: common checks")
report(abs(s$truth$local_ld$realized_mean_r2 - 0.3) < 0.05,
       "linked signals: realized LD matches target 0.3")

cat("== 5-7. R1 / R2 / R3 ==\n")
r1 <- mk("trait_representation", "highly_representable", 5)
r2 <- mk("trait_representation", "partially_representable", 6)
r3 <- mk("trait_representation", "strongly_nonredundant", 7)
tab1 <- r1$truth$subset_table
reps1 <- r1$truth$minimum_representative_sets
reps1 <- reps1[reps1$tolerance == tol, ]
report(min(reps1$set_size) == 1L,
       sprintf("R1: min Rep cardinality = 1 (loss(S*) = %.3f)",
               tab1$representation_loss[
                 tab1$representing_key == "Trait1"]))
tab2 <- r2$truth$subset_table
reps2 <- r2$truth$minimum_representative_sets
reps2 <- reps2[reps2$tolerance == tol, ]
report(all(tab2$representation_loss[tab2$set_size == 1L] > tol) &&
         min(reps2$set_size) == 2L,
       sprintf(paste0("R2: all singletons > tau, min pair = %.3f, ",
                      "min Rep cardinality = 2"),
               min(tab2$representation_loss[tab2$set_size == 2L])))
tab3 <- r3$truth$subset_table
k3 <- max(tab3$set_size)
min_proper <- min(tab3$representation_loss[tab3$set_size < k3])
report(min_proper > tol,
       sprintf("R3: all proper subsets > tau (min proper rho = %.3f)",
               min_proper))
# scale matching does not change the maps
r1b <- mk("trait_representation", "highly_representable", 5,
          locus_pve = 0.02)
report(isTRUE(all.equal(r1$truth$subset_table$representation_loss,
                        r1b$truth$subset_table$representation_loss,
                        tolerance = 1e-10)),
       "R1: scale matching leaves the rho map unchanged")

cat("== 8. linked pseudo-multitrait ==\n")
s <- mk("end_to_end", "linked_pseudo_multitrait", 8)
ck <- check_common(s)
report(all(unlist(ck)) && s$truth$signal_count == 2L &&
         s$truth$local_ld$realized_mean_r2 > 0.4,
       "linked_pseudo_multitrait: 2 signals, high LD")

cat("== 9. mixed multi-signal ==\n")
s <- mk("end_to_end", "mixed_multisignal", 9)
ck <- check_common(s)
report(all(unlist(ck)) && s$truth$signal_count == 2L,
       "mixed_multisignal: common checks, 2 signals")

cat("== 10. analysis modes ==\n")
mm <- CondPED:::.analysis_mode_mapping
report(identical(mm("full")$signal_mode, "resolve") &&
         identical(mm("signal_oracle")$candidate_mode, "holm_fwer") &&
         identical(mm("signal_trait_oracle")$candidate_mode,
                   "predefined"),
       "analysis-mode mapping frozen")
# signal_trait_oracle end-to-end smoke on a small replicate
s <- mk("end_to_end", "single_highly_representable", 10)
ps <- data.frame(
  locus_id = s$truth$loci$locus_id,
  representative_snp = s$truth$signals$representative_snp,
  stringsAsFactors = FALSE
)
ps$conditioning_snps <- list(character())
pre_sets <- stats::setNames(list(s$truth$candidate_traits[[1L]]),
                            s$truth$signals$signal_id)
f <- condped(
  s$Y, G = s$G, K = s$K_bg,
  signal_mode = "predefined", candidate_mode = "predefined",
  control = list(null_control = list(maxit = 300L),
                 predefined_signals = ps,
                 predefined_candidate_sets = pre_sets)
)
report(isTRUE(f$status$ok) &&
         identical(f$candidate_traits$candidate_sets[[1L]],
                   s$truth$candidate_traits[[1L]]),
       "signal_trait_oracle: truth sets flow through condped")

cat("== 11. ASSET adapter ==\n")
if (exists(".run_asset_comparison", envir = asNamespace("CondPED"))) {
  out <- CondPED:::.run_asset_comparison(
    beta = c(T1 = 0.3, T2 = -0.2),
    se = c(T1 = 0.1, T2 = 0.1),
    Sigma_Z = matrix(c(1, 0.4, 0.4, 1), 2),
    trait_names = c("T1", "T2")
  )
  report(identical(out$status, "asset_not_available") ||
           identical(out$status, "ok"),
         sprintf("ASSET adapter status: %s", out$status))
} else {
  report(FALSE, "ASSET adapter .run_asset_comparison() missing")
}

cat("\n")
if (length(failures) > 0L) {
  cat("VALIDATION FAILURES:\n")
  for (f in failures) cat(" -", f, "\n")
  quit(status = 1L)
} else {
  cat("All Stage 6C validation checks passed.\n")
}
