# Tests for condped() (v1.0 Stage 5): thin end-to-end orchestration.
#
# The pipeline pieces are tested elsewhere; here we test orchestration:
# exactly one null fit, no LOCO, legal empty/boundary results, P2
# switches never calling P2 stubs, output structure and
# reproducibility.

run_toy <- function(architecture = "candidate_pair", seed = 1,
                    n = 400L, m = 3L, p = 200L, locus_pve = 0.08, ...) {
  sim <- simulate_condped_data(n = n, m = m, p = p,
                               architecture = architecture,
                               locus_pve = locus_pve,
                               correlation = "block", seed = seed)
  list(sim = sim,
       fit = condped(sim$Y, G = sim$G, K = sim$K_bg,
                     control = list(null_control = list(maxit = 300L)),
                     ...))
}

test_that("minimal end-to-end run returns the contracted structure", {
  f <- run_toy()$fit
  expect_s3_class(f, "condped_fit")
  expect_identical(
    names(f),
    c("null_fit", "omnibus", "effects", "candidate_traits",
      "subset_analysis", "pve", "crossfit", "bootstrap", "settings",
      "status", "diagnostics")
  )
  expect_identical(names(f$subset_analysis),
                   c("subset_table", "minimum_representative_sets",
                     "irreducible_modules", "tolerance_path"))
  expect_identical(f$status$code, "ok")
  expect_s3_class(f$null_fit, "condped_mt_null")
  expect_true(is.data.frame(f$omnibus$omnibus))
  expect_true(is.data.frame(f$candidate_traits$trait_table))
  expect_true(is.data.frame(f$subset_analysis$subset_table))
  expect_null(f$pve)
  expect_null(f$crossfit)
  expect_null(f$bootstrap)
})

test_that("fit_mt_null is called exactly once", {
  sim <- simulate_condped_data(n = 100, m = 3L, p = 150L,
                               architecture = "candidate_pair",
                               correlation = "block", seed = 2)
  calls <- 0L
  real_fit <- fit_mt_null   # capture before rebinding
  local_mocked_bindings(
    fit_mt_null = function(...) {
      calls <<- calls + 1L
      real_fit(...)
    },
    .package = "CondPED"
  )
  f <- condped(sim$Y, G = sim$G, K = sim$K_bg,
               control = list(null_control = list(maxit = 200L)))
  expect_identical(calls, 1L)
})

test_that("P2 stubs are never called when the switches are FALSE", {
  sim <- simulate_condped_data(n = 100, m = 3L, p = 150L,
                               architecture = "candidate_pair",
                               correlation = "block", seed = 3)
  local_mocked_bindings(
    estimate_locus_pve = function(...) stop("PVE called"),
    crossfit_effect_pve = function(...) stop("crossfit called"),
    bootstrap_condped = function(...) stop("bootstrap called"),
    .package = "CondPED"
  )
  expect_no_error(
    condped(sim$Y, G = sim$G, K = sim$K_bg,
            control = list(null_control = list(maxit = 200L)))
  )
  # switching them on warns and keeps NULL, still without calling stubs
  expect_warning(
    expect_warning(
      f <- condped(sim$Y, G = sim$G, K = sim$K_bg, pve = TRUE,
                   bootstrap = TRUE,
                   control = list(null_control = list(maxit = 200L))),
      "not implemented"
    ),
    "not implemented"
  )
  expect_null(f$pve)
  expect_null(f$bootstrap)
})

test_that("no significant locus returns a complete empty result", {
  sim <- simulate_condped_data(n = 100, m = 3L, p = 150L,
                               architecture = "null",
                               correlation = "block", seed = 4)
  f <- condped(sim$Y, G = sim$G, K = sim$K_bg,
               alpha_omnibus = 1e-12, omnibus_adjust = "bonferroni",
               control = list(null_control = list(maxit = 200L)))
  expect_identical(f$status$code, "empty_selection")
  expect_true(f$status$ok)
  expect_null(f$effects)
  expect_identical(nrow(f$candidate_traits$trait_table), 0L)
  expect_identical(nrow(f$subset_analysis$subset_table), 0L)
  expect_identical(nrow(f$subset_analysis$minimum_representative_sets), 0L)
  expect_identical(nrow(f$subset_analysis$irreducible_modules), 0L)
})

test_that("single-candidate locus is a legal trait_restricted result", {
  sim <- simulate_condped_data(n = 300, m = 3L, p = 200L,
                               architecture = "candidate_single",
                               locus_pve = 0.05,
                               correlation = "block", seed = 5)
  f <- condped(sim$Y, G = sim$G, K = sim$K_bg,
               control = list(null_control = list(maxit = 300L)))
  expect_identical(f$status$code, "ok")
  lt <- f$candidate_traits$locus_table
  expect_true(all(lt$locus_status %in%
                    c("omnibus_only", "trait_restricted", "multi_trait")))
  # single-candidate loci appear in the subset table as not_applicable
  # and never reach formal extraction
  if (any(lt$locus_status == "trait_restricted")) {
    expect_true(all(
      f$subset_analysis$subset_table$status[
        f$subset_analysis$subset_table$marker_id %in%
          lt$marker_id[lt$locus_status == "trait_restricted"]
      ] == "not_applicable"
    ))
  }
})

test_that("multi-candidate locus produces formal set results", {
  sim <- simulate_condped_data(n = 400, m = 4L, p = 300L,
                               architecture = "candidate_dense",
                               locus_pve = 0.05,
                               correlation = "block", seed = 6)
  f <- condped(sim$Y, G = sim$G, K = sim$K_bg,
               control = list(null_control = list(maxit = 300L)))
  expect_identical(f$status$code, "ok")
  expect_gt(nrow(f$subset_analysis$subset_table), 0L)
  # 2^k rows for the causal locus when it is multi_trait
  qtl <- colnames(sim$G)[sim$causal_index]
  if (qtl %in% f$candidate_traits$locus_table$marker_id[
    f$candidate_traits$locus_table$locus_status == "multi_trait"]) {
    k <- sum(f$candidate_traits$candidate_sets[[qtl]] ==
               f$candidate_traits$candidate_sets[[qtl]]) # length
    k <- length(f$candidate_traits$candidate_sets[[qtl]])
    expect_equal(
      sum(f$subset_analysis$subset_table$marker_id == qtl), 2^k
    )
  }
})

test_that("manual followup_loci are honoured and validated", {
  sim <- simulate_condped_data(n = 150, m = 3L, p = 200L,
                               architecture = "candidate_pair",
                               locus_pve = 0.05,
                               correlation = "block", seed = 7)
  qtl <- colnames(sim$G)[sim$causal_index]
  f <- condped(sim$Y, G = sim$G, K = sim$K_bg,
               control = list(followup_loci = qtl,
                              null_control = list(maxit = 300L)))
  expect_identical(f$diagnostics$n_followup, 1L)
  expect_true(all(f$subset_analysis$subset_table$marker_id == qtl))
  expect_error(
    condped(sim$Y, G = sim$G, K = sim$K_bg,
            control = list(followup_loci = "no_such_marker")),
    class = "condped_invalid_input"
  )
})

test_that("K = NULL builds one fixed GRM; supplying the same K matches", {
  sim <- simulate_condped_data(n = 120, m = 3L, p = 150L,
                               architecture = "candidate_pair",
                               correlation = "block", seed = 8)
  f_null <- condped(sim$Y, G = sim$G, K = NULL, seed = 11L,
                    control = list(null_control = list(maxit = 200L)))
  # rebuild K exactly as condped() does
  maf_hat <- colMeans(sim$G) / 2
  sd_hat <- sqrt(2 * maf_hat * (1 - maf_hat))
  keep <- sd_hat > 0
  Z <- sweep(sim$G[, keep], 2L, 2 * maf_hat[keep], `-`)
  Z <- sweep(Z, 2L, sd_hat[keep], `/`)
  K_manual <- CondPED:::.make_grm(Z)
  f_manual <- condped(sim$Y, G = sim$G, K = K_manual, seed = 11L,
                      control = list(null_control = list(maxit = 200L)))
  expect_equal(f_null$null_fit$Sigma_P_ref, f_manual$null_fit$Sigma_P_ref,
               tolerance = 1e-8)
  expect_equal(f_null$omnibus$omnibus$Q, f_manual$omnibus$omnibus$Q,
               tolerance = 1e-8)
})

test_that("same seed reproduces the full result", {
  sim <- simulate_condped_data(n = 120, m = 3L, p = 150L,
                               architecture = "candidate_pair",
                               correlation = "block", seed = 9)
  f1 <- condped(sim$Y, G = sim$G, K = sim$K_bg, seed = 21L,
                control = list(null_control = list(maxit = 200L)))
  f2 <- condped(sim$Y, G = sim$G, K = sim$K_bg, seed = 21L,
                control = list(null_control = list(maxit = 200L)))
  expect_identical(f1$subset_analysis$subset_table,
                   f2$subset_analysis$subset_table)
  expect_identical(f1$omnibus$omnibus, f2$omnibus$omnibus)
})

test_that("chromosome is annotation only; unknown control fields warn", {
  sim <- simulate_condped_data(n = 100, m = 3L, p = 150L,
                               architecture = "candidate_pair",
                               correlation = "block", seed = 10)
  chr <- rep(c("chr1", "chr2"), length.out = ncol(sim$G))
  expect_warning(
    f <- condped(sim$Y, G = sim$G, K = sim$K_bg, chromosome = chr,
                 control = list(null_control = list(maxit = 200L),
                                bogus_field = 1)),
    "Unrecognised control field"
  )
  expect_identical(f$settings$chromosome, chr)
  expect_error(
    condped(sim$Y, G = sim$G, K = sim$K_bg, chromosome = c("chr1")),
    class = "condped_invalid_input"
  )
})

test_that("input validation", {
  sim <- simulate_condped_data(n = 60, m = 2L, p = 50L,
                               architecture = "null", seed = 12)
  expect_error(condped(sim$Y[, 1, drop = FALSE], G = sim$G, K = sim$K_bg,
                       alpha_omnibus = 2),
               class = "condped_invalid_input")
  expect_error(condped(sim$Y, G = sim$G[1:10, ], K = sim$K_bg),
               class = "condped_invalid_input")
})
