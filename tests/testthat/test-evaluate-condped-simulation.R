# Tests for evaluate_condped_simulation() (v1.0 Stage 6C).
#
# The eight hand-constructed cases of Runbook 8.4 are checked against
# hand-computed values under the new multi-signal truth/estimate schema
# and the v1.0 formal metric names. Set families are always compared as
# sets of normalised keys. Anchors 30-35 cover signal matching and
# 44-50 cover key-based merges and key-uniqueness failures.

eval_one <- CondPED:::.evaluate_one_replicate
match_signals <- CondPED:::.match_signals

# helpers to build minimal truth/estimate objects (new schema) ---------------

mk_sets_tab <- function(sets, tol = 0.10, locus_id = "L1",
                        signal_id = "S1") {
  if (length(sets) == 0L) {
    return(data.frame(
      locus_id = character(), signal_id = character(),
      tolerance = numeric(), trait_set = I(list()),
      trait_key = character(), set_size = integer(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    locus_id = locus_id, signal_id = signal_id,
    tolerance = tol,
    trait_set = I(sets),
    trait_key = vapply(sets, paste, character(1), collapse = "|"),
    set_size = lengths(sets),
    stringsAsFactors = FALSE
  )
}

mk_subset_tab <- function(keys, comps, etas, losses, locus_id = "L1",
                          signal_id = "S1") {
  if (length(keys) == 0L) {
    return(data.frame(
      locus_id = character(), signal_id = character(),
      representing_key = character(), complement_set = I(list()),
      complement_key = character(), conditional_effect = I(list()),
      representation_loss = numeric(), stringsAsFactors = FALSE
    ))
  }
  data.frame(
    locus_id = locus_id, signal_id = signal_id,
    representing_key = keys,
    complement_set = I(comps),
    complement_key = vapply(comps, function(x) {
      paste(x, collapse = "|")
    }, character(1)),
    conditional_effect = I(etas),
    representation_loss = losses,
    stringsAsFactors = FALSE
  )
}

mk_sigs <- function(locus_id, signal_id, order, snp, chr, pos) {
  data.frame(
    locus_id = locus_id, signal_id = signal_id, signal_order = order,
    representative_snp = snp, chromosome = chr, position = pos,
    stringsAsFactors = FALSE
  )
}

# single-signal truth: locus L1 on chr1:900-1100, signal S1 at M1 (1000)
mk_truth <- function(cand, reps = list(), mods = list(),
                     subset_table = NULL,
                     beta_vals = c(T1 = 1),
                     direction = "single_trait") {
  list(
    loci = data.frame(
      locus_id = "L1", chromosome = "1", start = 900, end = 1100,
      n_region_markers = 10L, n_signals = 1L, stringsAsFactors = FALSE
    ),
    signals = mk_sigs("L1", "S1", 1L, "M1", "1", 1000),
    causal_markers = "M1",
    signal_count = 1L,
    beta = data.frame(
      locus_id = "L1", signal_id = "S1", representative_snp = "M1",
      causal_marker = "M1", trait = names(beta_vals),
      beta = unname(beta_vals), stringsAsFactors = FALSE
    ),
    candidate_traits = list(S1 = cand),
    effect_breadth = c(S1 = length(cand)),
    effect_direction = c(S1 = direction),
    minimum_representative_sets = mk_sets_tab(reps),
    irreducible_modules = mk_sets_tab(mods),
    subset_table = subset_table
  )
}

# null architecture: no loci, no signals, no causal markers
mk_null_truth <- function() {
  list(
    loci = data.frame(
      locus_id = character(), chromosome = character(),
      start = numeric(), end = numeric(), n_region_markers = integer(),
      n_signals = integer(), stringsAsFactors = FALSE
    ),
    signals = mk_sigs(character(), character(), integer(),
                      character(), character(), numeric()),
    causal_markers = character(),
    signal_count = 0L,
    beta = data.frame(
      locus_id = character(), signal_id = character(),
      representative_snp = character(), causal_marker = character(),
      trait = character(), beta = numeric(), stringsAsFactors = FALSE
    ),
    candidate_traits = setNames(list(), character()),
    effect_breadth = setNames(integer(), character()),
    effect_direction = setNames(character(), character()),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = NULL
  )
}

# single-signal estimate: locus id "1:900-1100" (parsed), signal E1
mk_est <- function(cand, reps = list(), mods = list(), om = NULL,
                   subset_table = NULL, beta = NULL,
                   locus_id = "1:900-1100", signal_id = "E1",
                   snp = "M1", pos = 1000) {
  cs <- list()
  if (!is.null(cand)) cs <- setNames(list(cand), signal_id)
  list(
    signals = data.frame(
      locus_id = locus_id, signal_id = signal_id, signal_order = 1L,
      representative_snp = snp, position = pos, stringsAsFactors = FALSE
    ),
    beta = beta,
    candidate_sets = cs,
    minimum_representative_sets = mk_sets_tab(
      reps, locus_id = locus_id, signal_id = signal_id
    ),
    irreducible_modules = mk_sets_tab(
      mods, locus_id = locus_id, signal_id = signal_id
    ),
    subset_table = subset_table,
    omnibus_summary = om
  )
}

row_of <- function(truth, est, settings = list(tolerance = 0.10),
                   G = NULL) {
  eval_one(truth, est, settings = settings, runtime = 1, G = G)
}

# ---- cases 1-3: empty-set rules ---------------------------------------------

test_that("case 1: truth and estimate both empty", {
  r <- row_of(mk_truth(character()), mk_est(character()))
  expect_identical(r$candidate_exact_recovery, 1)
  expect_identical(r$candidate_jaccard, 1)
  expect_true(is.na(r$candidate_tpr))
  expect_identical(r$candidate_fdp, 0)
  expect_identical(r$candidate_precision, 1)
  expect_identical(r$rep_exact_family_recovery, 1)
  expect_identical(r$rep_tie_coverage, 1)
})

test_that("case 2: truth non-empty, estimate empty", {
  r <- row_of(mk_truth(c("T1", "T2")), mk_est(character()))
  expect_identical(r$candidate_tpr, 0)
  expect_identical(r$candidate_fdp, 0)
  expect_identical(r$candidate_jaccard, 0)
  expect_identical(r$candidate_exact_recovery, 0)
  expect_identical(r$rep_tie_coverage, 1)  # both rep families empty here
})

test_that("case 3: truth empty, estimate non-empty", {
  r <- row_of(mk_truth(character()), mk_est(c("T1")))
  expect_true(is.na(r$candidate_tpr))
  expect_identical(r$candidate_fdp, 1)
  expect_identical(r$candidate_precision, 0)
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
  expect_identical(r$rep_exact_family_recovery, 1)
  expect_identical(r$irr_family_recovery, 1)
  expect_identical(r$rep_min_cardinality_recovery, 1)
  expect_identical(r$rep_tie_coverage, 1)
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
  expect_identical(r$rep_exact_family_recovery, 0)
  expect_identical(r$rep_tie_coverage, 0.5)
  expect_identical(r$rep_min_cardinality_recovery, 1)  # both size 2
  expect_identical(r$overall_recovery, 0)
  # conditional on the (correct) candidate set the family was wrong
  expect_identical(r$conditional_recovery, 0)
})

# ---- case 6: joint module completely missed --------------------------------------

test_that("case 6: joint module completely missed", {
  truth <- mk_truth(c("T1", "T2", "T3"), mods = list(c("T1", "T2")))
  est <- mk_est(c("T1", "T2", "T3"), mods = list("T3"))
  r <- row_of(truth, est)
  expect_identical(r$irr_missed_rate, 1)
  expect_identical(r$irr_extra_split_rate, 0)  # T3 not a subset of T1|T2
  expect_identical(r$irr_family_recovery, 0)
})

# ---- case 7: joint module fragmented ----------------------------------------------

test_that("case 7: joint module split into sub-modules", {
  truth <- mk_truth(c("T1", "T2"), mods = list(c("T1", "T2")))
  est <- mk_est(c("T1", "T2"), mods = list("T1", "T2"))
  r <- row_of(truth, est)
  expect_identical(r$irr_missed_rate, 1)
  expect_identical(r$irr_extra_split_rate, 1)  # two strict sub-modules
  expect_identical(r$irr_family_recovery, 0)
})

# ---- case 8: multiple modules partially recovered ------------------------------------

test_that("case 8: multiple modules partially recovered", {
  truth <- mk_truth(c("T1", "T2", "T3"),
                    mods = list("T1", c("T2", "T3")))
  est <- mk_est(c("T1", "T2", "T3"),
                mods = list("T1"))
  r <- row_of(truth, est)
  expect_identical(r$irr_family_recovery, 0)
  expect_identical(r$irr_missed_rate, 1)         # T2|T3 missed
  expect_identical(r$irr_extra_split_rate, 0)    # no sub-modules of T2|T3
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
  # null scene: no causal markers, no power; every tested marker is null
  r0 <- row_of(mk_null_truth(), mk_est(character(), om = om))
  expect_true(is.na(r0$power_omnibus))
  expect_equal(r0$type1_omnibus, 2 / 4)  # M1 and N1 of M1,N1,N2,N3
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

# ---- Stage 6A quantitative metric anchors (new schema) ----------------------------

test_that("beta metrics match the hand anchor", {
  mk_beta <- function(se = NULL) {
    df <- data.frame(
      locus_id = "1:900-1100", signal_order = 1L,
      representative_snp = "M1", trait = c("T1", "T2"),
      beta = c(1.1, 0.4), stringsAsFactors = FALSE
    )
    if (!is.null(se)) df$se <- se
    df
  }
  truth <- mk_truth(c("T1", "T2"), beta_vals = c(T1 = 1.0, T2 = 0.5),
                    direction = "concordant",
                    subset_table = mk_subset_tab(character(), list(),
                                                 list(), numeric()))
  est <- mk_est(c("T1", "T2"),
                subset_table = mk_subset_tab(character(), list(), list(),
                                             numeric(),
                                             locus_id = "1:900-1100",
                                             signal_id = "E1"),
                beta = mk_beta(c(0.06, 0.06)))
  r <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_equal(r$beta_bias, 0)                    # (0.1 + -0.1) / 2
  expect_equal(r$beta_rmse, 0.1)                  # sqrt((0.01 + 0.01) / 2)
  expect_equal(r$beta_coverage, 1)                # 0.1 <= 1.96 * 0.06 twice
  expect_equal(r$beta_sign_accuracy, 1)

  est$beta <- mk_beta(c(0.05, 0.05))
  r2 <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_equal(r2$beta_coverage, 0)               # 0.1 > 1.96 * 0.05 twice
  est$beta <- mk_beta()
  r3 <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_true(is.na(r3$beta_coverage))
})

test_that("eta metrics match the hand anchor", {
  st_true <- mk_subset_tab(
    "A", list(c("B", "C")), list(c(B = 0.20, C = -0.10)), 0.08
  )
  st_est <- mk_subset_tab(
    "A", list(c("B", "C")), list(c(B = 0.25, C = -0.08)), 0.11,
    locus_id = "1:900-1100", signal_id = "E1"
  )
  truth <- mk_truth(c("A", "B", "C"), beta_vals = c(A = 1),
                    subset_table = st_true)
  est <- mk_est(c("A", "B", "C"), subset_table = st_est)
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
    list(c(B = 0.25, C = -0.08), NULL), c(0.11, 0.00),
    locus_id = "1:900-1100", signal_id = "E1"
  )
  truth <- mk_truth(c("A", "B", "C"), beta_vals = c(A = 1),
                    subset_table = st_true)
  est <- mk_est(c("A", "B", "C"), subset_table = st_est)
  r <- eval_one(truth, est, settings = list(tolerance = 0.10))
  expect_equal(r$representation_loss_bias, 0.015)   # (0.03 + 0)/2
  expect_equal(r$representation_loss_rmse, sqrt((0.03^2 + 0) / 2))
  expect_equal(r$representation_loss_mae, 0.015)
  expect_equal(r$representation_map_mae, 0.015)
  # truth side: feasible, feasible; est side: infeasible, feasible -> 1/2
  expect_equal(r$representation_threshold_accuracy, 0.5)
})

test_that("eta/rho metrics are NA when nothing matches", {
  st_true <- mk_subset_tab("A", list("B"), list(c(B = 0.2)), 0.1)
  st_est <- mk_subset_tab(character(), list(), list(), numeric(),
                          locus_id = "1:900-1100", signal_id = "E1")
  truth <- mk_truth("A", beta_vals = c(A = 1), subset_table = st_true)
  est <- mk_est("A", subset_table = st_est)
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
  expect_true("rep_exact_family_recovery" %in%
                res$settings$metric_levels$primary)
  expect_true("irr_family_recovery" %in%
                res$settings$metric_levels$secondary)
  expect_false(any(res$settings$metric_levels$primary %in%
                     res$settings$metric_levels$secondary))
})

# ---- Stage 6C anchors 30-35: signal matching -------------------------------------

# two-signal truth: L1 chr1:900-5100, S1 at M1 (1000), S2 at M2 (5000)
mk_truth2 <- function() {
  list(
    loci = data.frame(
      locus_id = "L1", chromosome = "1", start = 900, end = 5100,
      n_region_markers = 50L, n_signals = 2L, stringsAsFactors = FALSE
    ),
    signals = mk_sigs("L1", c("S1", "S2"), 1:2, c("M1", "M2"), "1",
                      c(1000, 5000)),
    causal_markers = c("M1", "M2"),
    signal_count = 2L,
    beta = data.frame(
      locus_id = "L1", signal_id = c("S1", "S2"),
      representative_snp = c("M1", "M2"), causal_marker = c("M1", "M2"),
      trait = "T1", beta = 1, stringsAsFactors = FALSE
    ),
    candidate_traits = list(S1 = "T1", S2 = "T1"),
    effect_breadth = c(S1 = 1L, S2 = 1L),
    effect_direction = c(S1 = "single_trait", S2 = "single_trait"),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = NULL
  )
}

mk_est2 <- function(pos = c(1000, 5000), signal_ids = c("E1", "E2"),
                    snps = c("M1", "M2"),
                    locus_id = "1:950-5050") {
  list(
    signals = data.frame(
      locus_id = locus_id, signal_id = signal_ids,
      signal_order = seq_along(signal_ids),
      representative_snp = snps, position = pos,
      stringsAsFactors = FALSE
    ),
    beta = NULL,
    candidate_sets = setNames(rep(list("T1"), length(signal_ids)),
                              signal_ids),
    minimum_representative_sets = mk_sets_tab(list()),
    irreducible_modules = mk_sets_tab(list()),
    subset_table = NULL,
    omnibus_summary = NULL
  )
}

test_that("30: exact causal representative matches with r2 = 1", {
  x <- rep(c(0, 1, 2, 1, 0, 2), length.out = 60)
  G <- cbind(M1 = x, P1 = (x + 1) %% 3)
  truth <- mk_truth(c("T1", "T2"))
  est <- mk_est(c("T1", "T2"))
  m <- match_signals(truth, est, G = G)
  expect_identical(nrow(m$matches), 1L)
  expect_identical(m$matches$truth_signal_id, "S1")
  expect_identical(m$matches$est_signal_id, "E1")
  expect_equal(m$matches$representative_to_causal_r2, 1)
  expect_equal(m$matches$representative_to_causal_distance, 0)
  r <- row_of(truth, est, G = G)
  expect_equal(r$representative_to_causal_r2, 1)
  expect_equal(r$representative_to_causal_distance, 0)
  # without G the r2 metric is NA but distance is still recorded
  r0 <- row_of(truth, est)
  expect_true(is.na(r0$representative_to_causal_r2))
  expect_equal(r0$representative_to_causal_distance, 0)
})

test_that("31: high-LD proxy representative matches the causal via r2", {
  x <- rep(c(0, 1, 2, 1, 0, 2), length.out = 60)
  y <- x
  y[c(2, 17, 33)] <- (x[c(2, 17, 33)] + 1) %% 3
  G <- cbind(M1 = x, P1 = y)
  truth <- mk_truth("T1")
  est <- mk_est("T1", snp = "P1", pos = 1050)
  m <- match_signals(truth, est, G = G)
  expect_identical(nrow(m$matches), 1L)
  expect_identical(m$matches$est_signal_id, "E1")
  expect_equal(m$matches$representative_to_causal_r2,
               stats::cor(G[, "M1"], G[, "P1"])^2)
  expect_equal(m$matches$representative_to_causal_distance, 50)
})

test_that("32: two-signal one-to-one matching, no signal used twice", {
  truth <- mk_truth2()
  est <- mk_est2(pos = c(1100, 4900))
  m <- match_signals(truth, est)
  expect_identical(nrow(m$matches), 2L)
  expect_setequal(m$matches$truth_signal_id, c("S1", "S2"))
  expect_setequal(m$matches$est_signal_id, c("E1", "E2"))
  expect_length(unique(m$matches$truth_signal_id), 2L)
  expect_length(unique(m$matches$est_signal_id), 2L)
  # nearest pairs win: E1 -> S1, E2 -> S2
  mm <- m$matches[order(m$matches$truth_signal_id), ]
  expect_identical(mm$est_signal_id, c("E1", "E2"))
  expect_equal(mm$representative_to_causal_distance, c(100, 100))
  expect_length(m$missed, 0L)
  expect_length(m$extra, 0L)
  r <- row_of(truth, est)
  expect_identical(r$signal_coverage, 1)
  expect_identical(r$secondary_signal_power, 1)
  expect_identical(r$signal_count_exact_recovery, 1)
  expect_identical(r$locus_detection_rate, 1)
})

test_that("33: unmatched truth signal is counted as missed", {
  truth <- mk_truth2()
  est <- mk_est2(pos = 1010, signal_ids = "E1", snps = "P1",
                 locus_id = "1:950-2000")
  m <- match_signals(truth, est)
  expect_identical(m$missed, "S2")
  expect_length(m$extra, 0L)
  r <- row_of(truth, est)
  expect_identical(r$missed_signal_rate, 0.5)
  expect_identical(r$signal_coverage, 0.5)
  expect_identical(r$secondary_signal_power, 0)
  expect_identical(r$extra_signal_rate, 0)
})

test_that("34: unmatched estimated signal is counted as extra", {
  truth <- mk_truth("T1")
  est <- mk_est2(pos = c(1000, 3000), locus_id = "1:900-3100")
  est$candidate_sets <- list(E1 = "T1", E2 = "T1")
  m <- match_signals(truth, est)
  expect_identical(m$extra, "E2")
  expect_length(m$missed, 0L)
  r <- row_of(truth, est)
  expect_identical(r$extra_signal_rate, 0.5)
  expect_identical(r$missed_signal_rate, 0)
  expect_identical(r$signal_coverage, 1)
  expect_identical(r$residual_local_false_positive_rate, 1)
  # estimated signal count (2) != truth signal count (1)
  expect_identical(r$signal_count_exact_recovery, 0)
})

test_that("35: matching is invariant to estimated signal row order", {
  truth <- mk_truth2()
  est <- mk_est2(pos = c(1100, 4900))
  m1 <- match_signals(truth, est)
  perm <- c(2L, 1L)
  est$signals <- est$signals[perm, ]
  m2 <- match_signals(truth, est)
  o1 <- order(m1$matches$truth_signal_id)
  o2 <- order(m2$matches$truth_signal_id)
  expect_equal(m1$matches[o1, ], m2$matches[o2, ])
  expect_identical(m1$missed, m2$missed)
  expect_identical(m1$extra, m2$extra)
})

# ---- Stage 6C anchors 44-48: key-based merges -------------------------------------

test_that("44: beta rows merge on locus + signal_order + trait keys", {
  truth <- mk_truth2()
  truth$beta <- data.frame(
    locus_id = "L1", signal_id = rep(c("S1", "S2"), each = 2L),
    representative_snp = rep(c("M1", "M2"), each = 2L),
    causal_marker = rep(c("M1", "M2"), each = 2L),
    trait = rep(c("T1", "T2"), 2L),
    beta = c(1.0, 0.5, 0.2, -0.1),
    stringsAsFactors = FALSE
  )
  truth$candidate_traits <- list(S1 = c("T1", "T2"), S2 = c("T1", "T2"))
  truth$effect_breadth <- c(S1 = 2L, S2 = 2L)
  # shuffled rows; per-pair diffs: S1 (T1 +0.2, T2 -0.1), S2 (-0.1, -0.1)
  est <- mk_est2()
  est$beta <- data.frame(
    locus_id = "1:950-5050",
    signal_order = c(2L, 1L, 2L, 1L),
    representative_snp = c("M2", "M1", "M2", "M1"),
    trait = c("T2", "T2", "T1", "T1"),
    beta = c(-0.2, 0.4, 0.1, 1.2),
    stringsAsFactors = FALSE
  )
  est$candidate_sets <- list(E1 = c("T1", "T2"), E2 = c("T1", "T2"))
  r <- row_of(truth, est)
  pair_bias <- c(mean(c(0.2, -0.1)), mean(c(-0.1, -0.1)))
  pair_rmse <- c(sqrt(mean(c(0.2^2, 0.1^2))), sqrt(mean(c(0.1^2, 0.1^2))))
  expect_equal(r$beta_bias, mean(pair_bias))
  expect_equal(r$beta_rmse, mean(pair_rmse))
  expect_equal(r$beta_sign_accuracy, 1)
})

test_that("45: eta rows merge per signal via signal_id + representing_key", {
  # both signals use representing_key "A" with different etas
  st_true <- rbind(
    mk_subset_tab("A", list("B"), list(c(B = 0.2)), 0.05,
                  signal_id = "S1"),
    mk_subset_tab("A", list("B"), list(c(B = 0.5)), 0.05,
                  signal_id = "S2")
  )
  st_est <- rbind(
    mk_subset_tab("A", list("B"), list(c(B = 0.25)), 0.05,
                  locus_id = "1:950-5050", signal_id = "E1"),
    mk_subset_tab("A", list("B"), list(c(B = 0.45)), 0.05,
                  locus_id = "1:950-5050", signal_id = "E2")
  )
  truth <- mk_truth2()
  truth$subset_table <- st_true
  est <- mk_est2()
  est$subset_table <- st_est
  r <- row_of(truth, est)
  # per-pair biases +0.05 and -0.05 cancel; rmse 0.05 in both pairs
  expect_equal(r$conditional_effect_bias, 0)
  expect_equal(r$conditional_effect_rmse, 0.05)
  expect_equal(r$conditional_effect_sign_accuracy, 1)
})

test_that("46: rho rows merge per signal with threshold side per pair", {
  st_true <- rbind(
    mk_subset_tab("A", list("B"), list(c(B = 0.2)), 0.08, signal_id = "S1"),
    mk_subset_tab("A", list("B"), list(c(B = 0.5)), 0.20, signal_id = "S2")
  )
  st_est <- rbind(
    mk_subset_tab("A", list("B"), list(c(B = 0.2)), 0.11,
                  locus_id = "1:950-5050", signal_id = "E1"),
    mk_subset_tab("A", list("B"), list(c(B = 0.5)), 0.18,
                  locus_id = "1:950-5050", signal_id = "E2")
  )
  truth <- mk_truth2()
  truth$subset_table <- st_true
  est <- mk_est2()
  est$subset_table <- st_est
  r <- row_of(truth, est)
  expect_equal(r$representation_loss_bias, mean(c(0.03, -0.02)))
  expect_equal(r$representation_loss_rmse,
               mean(c(0.03, 0.02)))
  expect_equal(r$representation_loss_mae, mean(c(0.03, 0.02)))
  # pair 1: truth feasible (0.08), est infeasible (0.11) -> 0
  # pair 2: both infeasible (0.20, 0.18) -> 1
  expect_equal(r$representation_threshold_accuracy, 0.5)
})

test_that("47: Rep families match per signal even with shared trait keys", {
  # both signals have a "T1|T2" minimum set; only S1's family is recovered
  truth <- mk_truth2()
  truth$minimum_representative_sets <- rbind(
    mk_sets_tab(list(c("T1", "T2")), signal_id = "S1"),
    mk_sets_tab(list(c("T1", "T2"), c("T1", "T3")), signal_id = "S2")
  )
  est <- mk_est2()
  est$minimum_representative_sets <- rbind(
    mk_sets_tab(list(c("T1", "T2")), locus_id = "1:950-5050",
                signal_id = "E1"),
    mk_sets_tab(list(c("T1", "T2")), locus_id = "1:950-5050",
                signal_id = "E2")
  )
  r <- row_of(truth, est)
  expect_equal(r$rep_exact_family_recovery, 0.5)  # S1 exact, S2 not
  expect_equal(r$rep_tie_coverage, mean(c(1, 0.5)))
  expect_equal(r$rep_min_cardinality_recovery, 1)  # both min sizes 2
})

test_that("48: Irr families match per signal", {
  truth <- mk_truth2()
  truth$irreducible_modules <- rbind(
    mk_sets_tab(list("T1"), signal_id = "S1"),
    mk_sets_tab(list(c("T1", "T2")), signal_id = "S2")
  )
  est <- mk_est2()
  est$irreducible_modules <- rbind(
    mk_sets_tab(list("T1"), locus_id = "1:950-5050", signal_id = "E1"),
    mk_sets_tab(list("T1"), locus_id = "1:950-5050", signal_id = "E2")
  )
  r <- row_of(truth, est)
  expect_equal(r$irr_family_recovery, 0.5)  # S1 exact, S2 not
  # S2's joint module is missed; S1 has no joint module (NA, excluded)
  expect_equal(r$irr_missed_rate, 1)
  expect_equal(r$irr_extra_split_rate, 0)  # only one sub-module of T1|T2
})

# ---- Stage 6C anchors 49-50: key uniqueness and metric names --------------------

test_that("49: duplicate keys are rejected with condped_invalid_input", {
  truth <- mk_truth("T1")
  truth$signals <- rbind(truth$signals, truth$signals)
  expect_error(
    row_of(truth, mk_est("T1")),
    class = "condped_invalid_input"
  )
  truth2 <- mk_truth("T1")
  truth2$beta <- rbind(truth2$beta, truth2$beta)  # same signal + trait twice
  expect_error(
    row_of(truth2, mk_est("T1")),
    class = "condped_invalid_input"
  )
  truth3 <- mk_truth("T1")
  truth3$minimum_representative_sets <- rbind(
    mk_sets_tab(list("T1"), signal_id = "S1"),
    mk_sets_tab(list("T1"), signal_id = "S1")  # duplicate family key
  )
  expect_error(
    row_of(truth3, mk_est("T1")),
    class = "condped_invalid_input"
  )
})

test_that("50: no effect_magnitude anywhere in outputs", {
  dir <- tempfile(); dir.create(dir)
  good <- list(
    truth = mk_truth(c("T1", "T2")),
    estimates = mk_est(c("T1", "T2")),
    settings = list(experiment = "sim6c", n = 500,
                    architecture = "rep_pair", locus_pve = 0.01,
                    correlation = "block", target_loss = 0.05,
                    tolerance = 0.10, alpha_omnibus = 0.05),
    status = list(ok = TRUE, code = "ok"),
    runtime = 1
  )
  saveRDS(good, file.path(dir, "r.rds"))
  res <- evaluate_condped_simulation(file.path(dir, "r.rds"))
  expect_false(any(grepl("effect_magnitude",
                         names(res$replicate_metrics))))
  expect_false(any(grepl("effect_magnitude", names(res$summary))))
  expect_false(any(grepl("effect_magnitude", names(res$mcse))))
  expect_false(any(grepl("effect_magnitude", res$settings$metrics)))
})

# ---- ASSET / ablation metrics -----------------------------------------------------

test_that("ASSET metrics are NA without estimates$asset and computed with it", {
  truth <- mk_truth(c("T1", "T2"), direction = "concordant")
  est <- mk_est(c("T1", "T2"))
  r0 <- row_of(truth, est)
  asset_metrics <- c(
    "asset_power", "asset_trait_set_jaccard", "asset_trait_set_precision",
    "asset_trait_set_recall", "asset_direction_recovery",
    "lead_trait_breadth_inflation", "resolved_trait_breadth_inflation",
    "pseudo_multitrait_interpretation_rate"
  )
  expect_true(all(asset_metrics %in% names(r0)))
  expect_true(all(is.na(r0[asset_metrics])))

  est$asset <- list(
    lead = list(sets = list(S1 = c("T1", "T2", "T3")),
                p = c(S1 = 0.001)),
    resolved = list(sets = list(S1 = c("T1", "T2")),
                    p = c(S1 = 0.001)),
    condped = list(sets = list(S1 = c("T1", "T2")),
                   p = c(S1 = 0.001),
                   direction = c(S1 = "concordant"))
  )
  r <- row_of(truth, est)
  expect_identical(r$asset_power, 1)
  expect_identical(r$asset_trait_set_jaccard, 1)
  expect_identical(r$asset_trait_set_precision, 1)
  expect_identical(r$asset_trait_set_recall, 1)
  expect_identical(r$asset_direction_recovery, 1)
  expect_identical(r$lead_trait_breadth_inflation, 1)      # 3 - 2
  expect_identical(r$resolved_trait_breadth_inflation, 0)  # 2 - 2
  expect_identical(r$pseudo_multitrait_interpretation_rate, 0)
})

# ---- effect breadth and direction ---------------------------------------------------

test_that("effect_breadth_bias and direction_recovery follow the frozen rule", {
  truth <- mk_truth(c("T1", "T2"), beta_vals = c(T1 = 1, T2 = 0.5),
                    direction = "concordant")
  # estimate adds a spurious candidate trait and both betas are positive
  est <- mk_est(c("T1", "T2", "T3"),
                beta = data.frame(
                  locus_id = "1:900-1100", signal_order = 1L,
                  representative_snp = "M1", trait = c("T1", "T2", "T3"),
                  beta = c(1.1, 0.4, 0.2), stringsAsFactors = FALSE
                ))
  r <- row_of(truth, est)
  expect_identical(r$effect_breadth_bias, 1)  # 3 - 2
  expect_identical(r$direction_recovery, 1)   # concordant == concordant

  # two candidate traits with opposite signs -> antagonistic
  truth2 <- mk_truth(c("T1", "T2"), beta_vals = c(T1 = 1, T2 = -0.5),
                     direction = "antagonistic")
  est2 <- mk_est(c("T1", "T2"),
                 beta = data.frame(
                   locus_id = "1:900-1100", signal_order = 1L,
                   representative_snp = "M1", trait = c("T1", "T2"),
                   beta = c(0.9, -0.6), stringsAsFactors = FALSE
                 ))
  r2 <- row_of(truth2, est2)
  expect_identical(r2$direction_recovery, 1)
  # estimated concordant vs truth antagonistic -> 0
  est2$beta$beta <- c(0.9, 0.6)
  r3 <- row_of(truth2, est2)
  expect_identical(r3$direction_recovery, 0)
})
