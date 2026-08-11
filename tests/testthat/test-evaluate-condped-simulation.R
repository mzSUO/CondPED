# Tests for evaluate_condped_simulation() (v1.0 Stage 6).
#
# All eight hand-constructed cases of Runbook 8.4 are checked against
# hand-computed values. Set families are always compared as sets of
# normalised keys.

eval_one <- CondPED:::.evaluate_one_replicate

# helpers to build minimal truth/estimate objects ---------------------------

mk_sets_tab <- function(sets, tol = 0.10) {
  if (length(sets) == 0L) {
    return(data.frame(
      tolerance = numeric(), trait_set = I(list()),
      trait_key = character(), set_size = integer(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    tolerance = tol,
    trait_set = I(sets),
    trait_key = vapply(sets, paste, character(1), collapse = "|"),
    set_size = lengths(sets),
    stringsAsFactors = FALSE
  )
}

mk_truth <- function(cand, reps = list(), mods = list()) {
  list(
    beta = matrix(1, nrow = 1L, dimnames = list("M1", "T1")),
    candidate_traits = cand,
    minimum_representative_sets = mk_sets_tab(reps),
    irreducible_modules = mk_sets_tab(mods)
  )
}

mk_est <- function(cand, reps = list(), mods = list(),
                   om = NULL) {
  cs <- list()
  if (!is.null(cand)) cs <- list(M1 = cand)
  list(
    causal_marker = "M1",
    omnibus_summary = om,
    candidate_sets = cs,
    minimum_representative_sets = mk_sets_tab(reps),
    irreducible_modules = mk_sets_tab(mods)
  )
}

row_of <- function(truth, est, settings = list(tolerance = 0.10)) {
  eval_one(truth, est, settings = settings, runtime = 1)
}

# ---- cases 1-3: empty-set rules ---------------------------------------------

test_that("case 1: truth and estimate both empty", {
  r <- row_of(mk_truth(character()), mk_est(character()))
  expect_identical(r$candidate_exact_recovery, 1)
  expect_identical(r$candidate_jaccard, 1)
  expect_true(is.na(r$candidate_tpr))
  expect_identical(r$candidate_fdp, 0)
  expect_identical(r$representative_exact_recovery, 1)
  expect_identical(r$tie_coverage, 1)
})

test_that("case 2: truth non-empty, estimate empty", {
  r <- row_of(mk_truth(c("T1", "T2")), mk_est(character()))
  expect_identical(r$candidate_tpr, 0)
  expect_identical(r$candidate_fdp, 0)
  expect_identical(r$candidate_jaccard, 0)
  expect_identical(r$candidate_exact_recovery, 0)
  expect_identical(r$tie_coverage, 1)  # both rep families empty here
})

test_that("case 3: truth empty, estimate non-empty", {
  r <- row_of(mk_truth(character()), mk_est(c("T1")))
  expect_true(is.na(r$candidate_tpr))
  expect_identical(r$candidate_fdp, 1)
  expect_identical(r$candidate_jaccard, 0)
  expect_identical(r$candidate_exact_recovery, 0)
})

# ---- case 4: exact recovery ----------------------------------------------------

test_that("case 4: single set exactly recovered", {
  truth <- mk_truth(c("T1", "T2"),
                    reps = list(c("T1", "T2")),
                    mods = list("T1"))
  est <- mk_est(c("T1", "T2"),
                reps = list(c("T1", "T2")),
                mods = list("T1"))
  r <- row_of(truth, est)
  expect_identical(r$candidate_exact_recovery, 1)
  expect_identical(r$representative_exact_recovery, 1)
  expect_identical(r$irreducible_exact_recovery, 1)
  expect_equal(r$representative_size_bias, 0)
  expect_identical(r$tie_coverage, 1)
  expect_identical(r$overall_recovery, 1)
  expect_identical(r$conditional_recovery, 1)
})

# ---- case 5: tied representative sets, only one recovered -----------------------

test_that("case 5: two tied minimum sets, one recovered", {
  truth <- mk_truth(c("T1", "T2", "T3"),
                    reps = list(c("T1", "T2"), c("T1", "T3")))
  est <- mk_est(c("T1", "T2", "T3"),
                reps = list(c("T1", "T2")))
  r <- row_of(truth, est)
  expect_identical(r$representative_exact_recovery, 0)
  expect_equal(r$representative_jaccard, 1 / 2)   # union is {AB, AC}
  expect_identical(r$tie_coverage, 0.5)
  expect_equal(r$representative_size_bias, 0)  # both size 2
  expect_identical(r$overall_recovery, 0)
  # conditional on the (correct) candidate set the family was wrong
  expect_identical(r$conditional_recovery, 0)
})

# ---- case 6: joint module completely missed --------------------------------------

test_that("case 6: joint module completely missed", {
  truth <- mk_truth(c("T1", "T2", "T3"), mods = list(c("T1", "T2")))
  est <- mk_est(c("T1", "T2", "T3"), mods = list("T3"))
  r <- row_of(truth, est)
  expect_identical(r$joint_module_miss_rate, 1)
  expect_identical(r$over_fragmentation_rate, 0)  # T3 is not a subset of T1|T2
  expect_identical(r$irreducible_tpr, 0)
  expect_identical(r$irreducible_fdp, 1)          # est T3 not in truth
  expect_identical(r$irreducible_exact_recovery, 0)
})

# ---- case 7: joint module fragmented ----------------------------------------------

test_that("case 7: joint module split into sub-modules", {
  truth <- mk_truth(c("T1", "T2"), mods = list(c("T1", "T2")))
  est <- mk_est(c("T1", "T2"), mods = list("T1", "T2"))
  r <- row_of(truth, est)
  expect_identical(r$joint_module_miss_rate, 1)
  expect_identical(r$over_fragmentation_rate, 1)  # two strict sub-modules
  expect_identical(r$irreducible_tpr, 0)
  expect_identical(r$irreducible_fdp, 1)
})

# ---- case 8: multiple modules partially recovered ------------------------------------

test_that("case 8: multiple modules partially recovered", {
  truth <- mk_truth(c("T1", "T2", "T3"),
                    mods = list("T1", c("T2", "T3")))
  est <- mk_est(c("T1", "T2", "T3"),
                mods = list("T1"))
  r <- row_of(truth, est)
  expect_identical(r$irreducible_tpr, 0.5)
  expect_identical(r$irreducible_fdp, 0)
  expect_equal(r$irreducible_jaccard, 1 / 2)   # union is {T1, T2|T3}
  expect_identical(r$irreducible_exact_recovery, 0)
  expect_identical(r$joint_module_miss_rate, 1)   # T2|T3 missed
  expect_identical(r$over_fragmentation_rate, 0)  # no sub-modules of T2|T3
})

# ---- omnibus metrics ------------------------------------------------------------------

test_that("type1 and power use the omnibus summary", {
  om <- data.frame(
    marker_id = c("M1", "N1", "N2", "N3", "N4"),
    p_value = c(0.001, 0.01, 0.20, 0.60, NA),
    stringsAsFactors = FALSE
  )
  truth <- mk_truth(c("T1", "T2"))
  est <- mk_est(c("T1", "T2"), om = om)
  r <- row_of(truth, est, settings = list(tolerance = 0.10,
                                          alpha_omnibus = 0.05))
  expect_identical(r$power_omnibus, 1)
  expect_equal(r$type1_omnibus, 1 / 3)  # N1 of N1,N2,N3 (NA excluded)
  # null scene: no power
  r0 <- row_of(mk_truth(character()), mk_est(character(), om = om))
  expect_true(is.na(r0$power_omnibus))
  expect_equal(r0$type1_omnibus, 1 / 3)
})

# ---- public function: files, failures, grouping -------------------------------------------

test_that("public function reads files, summarises groups and lists failures", {
  dir <- tempfile(); dir.create(dir)
  good <- list(
    truth = mk_truth(c("T1", "T2"), reps = list(c("T1", "T2")),
                     mods = list("T1")),
    estimates = mk_est(c("T1", "T2"), reps = list(c("T1", "T2")),
                       mods = list("T1")),
    settings = list(experiment = "sim2", n = 500, architecture = "rep_pair",
                    locus_pve = 0.01, correlation = "block",
                    target_loss = 0.05, tolerance = 0.10,
                    alpha_omnibus = 0.05),
    status = list(ok = TRUE, code = "ok"),
    runtime = 12
  )
  unstable <- good
  unstable$status <- list(ok = FALSE, code = "unstable")
  unstable$estimates <- mk_est(c("T1"), reps = list(), mods = list())
  saveRDS(good, file.path(dir, "r1.rds"))
  saveRDS(good, file.path(dir, "r2.rds"))
  saveRDS(unstable, file.path(dir, "r3.rds"))
  saveRDS(list(junk = 1), file.path(dir, "broken.rds"))

  files <- file.path(dir, c("r1.rds", "r2.rds", "r3.rds", "broken.rds"))
  res <- evaluate_condped_simulation(files)
  expect_identical(res$status$code, "ok")
  expect_identical(nrow(res$replicate_metrics), 3L)
  expect_identical(nrow(res$failures), 2L)   # unstable + broken
  expect_true(res$failures$included[res$failures$reason ==
                                      "replicate status: unstable"])
  expect_equal(res$summary$overall_recovery, 2 / 3)
  expect_equal(res$mcse$overall_recovery,
               sqrt((2 / 3) * (1 / 3) / 3), tolerance = 1e-10)
  expect_true("architecture" %in% names(res$summary))

  # include_unstable = FALSE drops the replicate from metrics only
  res2 <- evaluate_condped_simulation(files, include_unstable = FALSE)
  expect_identical(nrow(res2$replicate_metrics), 2L)
  expect_identical(nrow(res2$failures), 2L)
  expect_equal(res2$summary$overall_recovery, 1)

  # metrics selection
  res3 <- evaluate_condped_simulation(files[1:2],
                                      metrics = c("candidate_tpr", "runtime"))
  expect_true(all(c("candidate_tpr", "runtime") %in% names(res3$summary)))
  expect_false("power_omnibus" %in% names(res3$summary))
})

test_that("edge behaviours: no evaluable files, group_by fallback", {
  dir <- tempfile(); dir.create(dir)
  saveRDS(list(junk = 1), file.path(dir, "broken.rds"))
  res <- evaluate_condped_simulation(file.path(dir, "broken.rds"))
  expect_identical(res$status$code, "empty_selection")
  expect_true(res$status$ok)
  expect_identical(nrow(res$failures), 1L)

  saveRDS(list(truth = mk_truth("T1"), estimates = mk_est("T1"),
               settings = list(tolerance = 0.10),
               status = list(ok = TRUE, code = "ok"), runtime = 1),
          file.path(dir, "r.rds"))
  expect_warning(
    res2 <- evaluate_condped_simulation(file.path(dir, "r.rds"),
                                        group_by = c("architecture", "n")),
    "dropped"
  )
  expect_identical(nrow(res2$summary), 1L)
})

# ---- Stage 6A quantitative metric anchors -----------------------------------

mk_subset_tab <- function(keys, comps, etas, losses) {
  if (length(keys) == 0L) {
    return(data.frame(
      marker_id = character(), representing_key = character(),
      complement_set = I(list()), conditional_effect = I(list()),
      representation_loss = numeric(), stringsAsFactors = FALSE
    ))
  }
  data.frame(
    marker_id = "M1",
    representing_key = keys,
    complement_set = I(comps),
    conditional_effect = I(etas),
    representation_loss = losses,
    stringsAsFactors = FALSE
  )
}

test_that("beta metrics match the hand anchor", {
  truth <- list(
    beta = matrix(c(1.0, 0.5), nrow = 1L,
                  dimnames = list("M1", c("T1", "T2"))),
    candidate_traits = c("T1", "T2"),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = mk_subset_tab(character(), list(), list(), numeric())
  )
  est <- list(
    causal_marker = "M1",
    omnibus_summary = NULL,
    candidate_sets = list(M1 = c("T1", "T2")),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = mk_subset_tab(character(), list(), list(), numeric()),
    beta = c(T1 = 1.1, T2 = 0.4),
    se = c(T1 = 0.06, T2 = 0.06)
  )
  r <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_equal(r$beta_bias, 0)                    # (0.1 + -0.1) / 2
  expect_equal(r$beta_rmse, 0.1)                  # sqrt((0.01 + 0.01) / 2)
  expect_equal(r$beta_coverage, 1)                # 0.1 <= 1.96 * 0.06 twice
  expect_equal(r$beta_sign_accuracy, 1)

  est$se <- c(T1 = 0.05, T2 = 0.05)
  r2 <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_equal(r2$beta_coverage, 0)               # 0.1 > 1.96 * 0.05 twice
  est$se <- NULL
  r3 <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_true(is.na(r3$beta_coverage))
})

test_that("eta metrics match the hand anchor", {
  st_true <- mk_subset_tab(
    "A", list(c("B", "C")), list(c(B = 0.20, C = -0.10)), 0.08
  )
  st_est <- mk_subset_tab(
    "A", list(c("B", "C")), list(c(B = 0.25, C = -0.08)), 0.11
  )
  truth <- list(
    beta = matrix(1, nrow = 1L, dimnames = list("M1", "A")),
    candidate_traits = c("A", "B", "C"),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = st_true
  )
  est <- list(
    causal_marker = "M1",
    omnibus_summary = NULL,
    candidate_sets = list(M1 = c("A", "B", "C")),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = st_est
  )
  r <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_equal(r$conditional_effect_bias, 0.035)          # (0.05 + 0.02)/2
  expect_equal(r$conditional_effect_rmse, sqrt(0.00145))
  expect_equal(r$conditional_effect_sign_accuracy, 1)
})

test_that("rho metrics match the hand anchor, including threshold side", {
  st_true <- mk_subset_tab(
    c("A", "A|B"), list(c("B", "C"), character()),
    list(c(B = 0.2, C = -0.1), NULL), c(0.08, 0.00)
  )
  st_est <- mk_subset_tab(
    c("A", "A|B"), list(c("B", "C"), character()),
    list(c(B = 0.25, C = -0.08), NULL), c(0.11, 0.00)
  )
  truth <- list(
    beta = matrix(1, nrow = 1L, dimnames = list("M1", "A")),
    candidate_traits = c("A", "B", "C"),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = st_true
  )
  est <- list(
    causal_marker = "M1",
    omnibus_summary = NULL,
    candidate_sets = list(M1 = c("A", "B", "C")),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = st_est
  )
  r <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_equal(r$representation_loss_bias, 0.015)   # (0.03 + 0)/2
  expect_equal(r$representation_loss_rmse, sqrt((0.03^2 + 0) / 2))
  expect_equal(r$representation_loss_mae, 0.015)
  expect_equal(r$representation_map_mae, 0.015)
  # truth side: feasible, feasible; est side: infeasible, feasible -> 1/2
  expect_equal(r$representation_threshold_accuracy, 0.5)
})

test_that("eta/rho metrics are NA when nothing matches", {
  truth <- list(
    beta = matrix(1, nrow = 1L, dimnames = list("M1", "A")),
    candidate_traits = "A",
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = mk_subset_tab("A", list("B"), list(c(B = 0.2)), 0.1)
  )
  est <- list(
    causal_marker = "M1",
    omnibus_summary = NULL,
    candidate_sets = list(M1 = "A"),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = mk_subset_tab(character(), list(), list(), numeric())
  )
  r <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_true(is.na(r$conditional_effect_bias))
  expect_true(is.na(r$representation_loss_bias))
  expect_true(is.na(r$representation_threshold_accuracy))
})

test_that("Rep is marked primary and Irr secondary in settings", {
  dir <- tempfile(); dir.create(dir)
  good <- list(
    truth = mk_truth("T1"),
    estimates = mk_est("T1"),
    settings = list(experiment = "sim6a", n = 500,
                    architecture = "rep_pair", locus_pve = 0.01,
                    correlation = "block", target_loss = 0.05,
                    tolerance = 0.10, alpha_omnibus = 0.05),
    status = list(ok = TRUE, code = "ok"),
    runtime = 1
  )
  saveRDS(good, file.path(dir, "r.rds"))
  res <- evaluate_condped_simulation(file.path(dir, "r.rds"))
  expect_true("representative_exact_recovery" %in%
                res$settings$metric_levels$primary)
  expect_true("irreducible_tpr" %in%
                res$settings$metric_levels$secondary)
  expect_false(any(res$settings$metric_levels$primary %in%
                     res$settings$metric_levels$secondary))
})
