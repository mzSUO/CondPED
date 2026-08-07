# Tests for the public decompose_conditional_effects() and
# derive_conditional_contrasts() (v1.0 Stage 2).
#
# The toy example is hand-verified: A and B are correlated (0.9) with
# beta_B exactly explained by A; C is independent.

Sigma_toy <- matrix(c(1.0, 0.9, 0.0,
                      0.9, 1.0, 0.0,
                      0.0, 0.0, 1.0), 3, 3,
                    dimnames = list(c("A", "B", "C"), c("A", "B", "C")))

toy_basis <- function() derive_conditional_contrasts(Sigma_toy)

toy_effects <- function() {
  data.frame(
    marker_id = rep(c("M1", "M2", "M3"), each = 3),
    trait = rep(c("A", "B", "C"), 3),
    beta = c(1.0, 0.9, 0.5,    # M1: C carries an independent effect
             1.0, 0.9, 0.02,   # M2: C ~no residual effect
             0.0, 0.0, 1.0),   # M3: single-trait locus
    se = 0.1, stringsAsFactors = FALSE
  )
}

toy_attribution <- function() {
  list(M1 = c("A", "B", "C"), M2 = c("A", "B", "C"), M3 = "C")
}

# ---- end-to-end structure -----------------------------------------------------

test_that("subset_table has 2^k rows per locus and contract columns", {
  res <- decompose_conditional_effects(
    toy_effects(), toy_attribution(), toy_basis(), tolerance = 0.10
  )
  expect_identical(res$status$code, "ok")
  tab <- res$subset_table
  expect_identical(names(tab),
                   c("marker_id", "set_id", "analysis_scope",
                     "representing_set", "representing_key",
                     "complement_set", "complement_key", "set_size",
                     "conditional_effect",
                     "full_qform", "subset_qform", "residual_qform",
                     "representation_loss", "decomposition_error",
                     "feasible_primary", "rank_Sigma_SS", "rank_Omega",
                     "condition_Sigma_SS", "condition_Omega",
                     "used_pseudoinverse", "status"))
  # M1/M2: 2^3 = 8 rows; M3 is single-candidate: 2 rows, not_applicable
  expect_equal(sum(tab$marker_id == "M1"), 8)
  expect_equal(sum(tab$marker_id == "M2"), 8)
  m3 <- tab[tab$marker_id == "M3", ]
  expect_equal(nrow(m3), 2)
  expect_true(all(m3$status == "not_applicable"))
  expect_identical(res$diagnostics$n_single_candidate, 1L)
  # no v0.6 profile/partition objects anywhere
  expect_false(any(c("profile_table", "partition_table", "R",
                     "persistent_modules") %in% names(res)))
})

test_that("losses match the unique Stage-1 path value by value", {
  res <- decompose_conditional_effects(
    toy_effects(), toy_attribution(), toy_basis(), tolerance = 0.10
  )
  tab <- res$subset_table[res$subset_table$marker_id == "M1", ]
  beta <- c(A = 1.0, B = 0.9, C = 0.5)
  for (i in seq_len(nrow(tab))) {
    ref <- CondPED:::.compute_subset_loss(
      beta, Sigma_toy, tab$representing_set[[i]], trait_names = c("A", "B", "C")
    )
    expect_equal(tab$representation_loss[i], ref$representation_loss,
                 tolerance = 1e-12)
    expect_equal(tab$full_qform[i], ref$full_qform, tolerance = 1e-12)
  }
  # hand check: S={A} leaves residual only on C: 0.25/1.25 = 0.2
  rowA <- tab[tab$representing_key == "A", ]
  expect_equal(rowA$representation_loss, 0.2, tolerance = 1e-10)
})

test_that("minimum representative sets and modules are correct end to end", {
  res <- decompose_conditional_effects(
    toy_effects(), toy_attribution(), toy_basis(),
    tolerance = 0.15, sensitivity_tolerance = 0.15
  )
  reps <- res$minimum_representative_sets
  m1_reps <- reps[reps$marker_id == "M1" & reps$tolerance == 0.15, ]
  # M1: feasible size-1 sets? loss(A)=0.2 > 0.15, loss(B) larger,
  # loss(C)=0.8; size-2: A|C loss 0 -> unique minimum "A|C"
  expect_identical(m1_reps$trait_key, "A|C")
  expect_identical(m1_reps$n_tied_solutions, 1L)
  m2_reps <- reps[reps$marker_id == "M2" & reps$tolerance == 0.15, ]
  expect_identical(m2_reps$trait_key, "A")   # C residual is tiny
  # M3 (single candidate) is excluded from extraction
  expect_false("M3" %in% reps$marker_id)

  mods <- res$irreducible_modules
  m1_mods <- mods[mods$marker_id == "M1" & mods$tolerance == 0.15, ]
  # loss(B|C) = 0.152 > 0.15 -> {A} module; loss(A|B) = 0.2 > 0.15 -> {C}
  expect_identical(sort(m1_mods$trait_key), c("A", "C"))
})

test_that("all and singleton modes agree on their common subsets", {
  common <- function(mode) {
    res <- decompose_conditional_effects(
      toy_effects(), toy_attribution(), toy_basis(),
      subset_mode = mode, tolerance = 0.10
    )
    tab <- res$subset_table[res$subset_table$marker_id == "M1", ]
    setNames(tab$representation_loss, tab$representing_key)
  }
  all_l <- common("all")
  sin_l <- common("singleton")
  shared <- intersect(names(all_l), names(sin_l))
  expect_gt(length(shared), 0L)
  expect_equal(sin_l[shared], all_l[shared], tolerance = 1e-12)
  # singleton mode performs no formal extraction
  res_sin <- decompose_conditional_effects(
    toy_effects(), toy_attribution(), toy_basis(),
    subset_mode = "singleton", tolerance = 0.10
  )
  expect_identical(nrow(res_sin$minimum_representative_sets), 0L)
  expect_identical(nrow(res_sin$irreducible_modules), 0L)
  expect_true(all(res_sin$tolerance_path$status %in%
                    c("not_applicable")))
})

test_that("custom mode evaluates exactly the requested sets plus boundaries", {
  res <- decompose_conditional_effects(
    toy_effects(), toy_attribution(), toy_basis(),
    subset_mode = "custom", custom_sets = list(c("A", "B")),
    tolerance = 0.10
  )
  tab <- res$subset_table[res$subset_table$marker_id == "M1", ]
  expect_identical(sort(tab$representing_key),
                   sort(c("<empty>", "A|B", "A|B|C")))
  # no formal extraction, and custom must not claim global optima
  expect_identical(nrow(res$minimum_representative_sets), 0L)
  expect_true(all(res$tolerance_path$status == "not_applicable"))
})

test_that("restrict_to_candidates = FALSE uses all traits (oracle scope)", {
  att <- list(M1 = "A", M2 = "A", M3 = "C")   # pretend tiny candidate sets
  res <- decompose_conditional_effects(
    toy_effects(), att, toy_basis(),
    restrict_to_candidates = FALSE, tolerance = 0.10
  )
  tab <- res$subset_table
  expect_true(all(tab$analysis_scope == "all_traits"))
  # all three loci now evaluate the full 3-trait space (M3 too: k = 3)
  expect_equal(sum(tab$marker_id == "M3"), 8)
})

test_that("too_many_traits is a legal status, never a silent demotion", {
  m <- 11L
  tn <- paste0("T", seq_len(m))
  S <- diag(m); dimnames(S) <- list(tn, tn)
  basis <- derive_conditional_contrasts(S)
  eff <- data.frame(marker_id = rep("BIG", m), trait = tn,
                    beta = 1, se = 0.1, stringsAsFactors = FALSE)
  res <- decompose_conditional_effects(
    eff, list(BIG = tn), basis, subset_mode = "all", tolerance = 0.1
  )
  expect_identical(nrow(res$subset_table), 0L)
  expect_true(all(res$tolerance_path$status == "too_many_traits"))
  expect_identical(res$diagnostics$n_too_many_traits, 1L)
  # singleton mode still works for the same locus
  res2 <- decompose_conditional_effects(
    eff, list(BIG = tn), basis, subset_mode = "singleton", tolerance = 0.1
  )
  expect_gt(nrow(res2$subset_table), 0L)
})

test_that("empty and missing candidate sets are handled explicitly", {
  # no marker has candidates -> legal empty result
  res <- decompose_conditional_effects(
    toy_effects(), list(M1 = character(), M2 = character(),
                        M3 = character()),
    toy_basis(), tolerance = 0.10
  )
  expect_identical(res$status$code, "empty_selection")
  expect_true(res$status$ok)
  expect_identical(nrow(res$subset_table), 0L)
  # markers absent from attribution count as empty candidates
  res2 <- decompose_conditional_effects(
    toy_effects(), list(M1 = c("A", "B", "C")), toy_basis(),
    tolerance = 0.10
  )
  expect_identical(res2$diagnostics$n_empty_candidates, 2L)
  expect_true(all(res2$subset_table$marker_id == "M1"))
})

test_that("return_subset_table = FALSE drops the table but keeps results", {
  res <- decompose_conditional_effects(
    toy_effects(), toy_attribution(), toy_basis(), tolerance = 0.15,
    return_subset_table = FALSE
  )
  expect_null(res$subset_table)
  expect_gt(nrow(res$minimum_representative_sets), 0L)
})

test_that("input validation", {
  expect_error(decompose_conditional_effects(
    toy_effects(), NULL, toy_basis()), class = "condped_invalid_input")
  expect_error(decompose_conditional_effects(
    toy_effects(), toy_attribution(), list(Sigma_P = Sigma_toy)),
    class = "condped_invalid_input")
  expect_error(decompose_conditional_effects(
    toy_effects(), toy_attribution(), toy_basis(), tolerance = 1.5),
    class = "condped_invalid_input")
  bad_eff <- toy_effects(); bad_eff$beta[1] <- NA
  expect_error(decompose_conditional_effects(
    bad_eff, toy_attribution(), toy_basis()),
    class = "condped_invalid_input")
  # attribution=NULL is legal in oracle mode
  expect_no_error(decompose_conditional_effects(
    toy_effects(), NULL, toy_basis(), restrict_to_candidates = FALSE
  ))
})

# ---- derive_conditional_contrasts ---------------------------------------------

test_that("basis object: fields, class and diagnostics", {
  basis <- toy_basis()
  expect_s3_class(basis, "condped_representation_basis")
  expect_identical(names(basis),
                   c("Sigma_P_ref", "trait_names", "basis", "status",
                     "diagnostics"))
  expect_identical(basis$trait_names, c("A", "B", "C"))
  expect_equal(basis$Sigma_P_ref, Sigma_toy)
  expect_identical(basis$diagnostics$rank, 3L)
  expect_false(basis$diagnostics$used_pseudoinverse)
  expect_identical(basis$basis$ridge, 0)
  expect_false(basis$basis$standardize)
})

test_that("basis: ridge only on request, standardize frozen off", {
  basis <- derive_conditional_contrasts(Sigma_toy, ridge = 0.5)
  expect_equal(diag(basis$Sigma_P_ref), diag(Sigma_toy) + 0.5)
  expect_match(paste(basis$status$warnings, collapse = " "), "ridge")
  corr_basis <- derive_conditional_contrasts(Sigma_toy, standardize = TRUE)
  expect_equal(unname(diag(corr_basis$Sigma_P_ref)), rep(1, 3))
  expect_equal(corr_basis$Sigma_P_ref["A", "B"], 0.9)
})

test_that("basis: rank-deficient Sigma_P is reported, not hidden", {
  S <- matrix(c(1, 1, 1, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  basis <- derive_conditional_contrasts(S)
  expect_identical(basis$status$code, "rank_deficient")
  expect_false(basis$status$ok)
  expect_identical(basis$diagnostics$rank, 1L)
  expect_true(basis$diagnostics$used_pseudoinverse)
})

test_that("basis input validation", {
  expect_error(derive_conditional_contrasts(matrix(c(1, 2, 3, 4), 2)),
               class = "condped_invalid_input")
  expect_error(derive_conditional_contrasts(Sigma_toy, ridge = -1),
               class = "condped_invalid_input")
  expect_error(derive_conditional_contrasts(Sigma_toy,
                                            trait_names = c("A", "A", "B")),
               class = "condped_invalid_input")
})

test_that("conditional_effect carries eta_R|S from the unique loss path", {
  res <- decompose_conditional_effects(
    toy_effects(), toy_attribution(), toy_basis(), tolerance = 0.10
  )
  tab <- res$subset_table[res$subset_table$marker_id == "M1", ]
  beta <- c(A = 1.0, B = 0.9, C = 0.5)
  for (i in seq_len(nrow(tab))) {
    ref <- CondPED:::.compute_subset_loss(
      beta, Sigma_toy, tab$representing_set[[i]], trait_names = c("A", "B", "C")
    )
    expect_identical(tab$conditional_effect[[i]], ref$eta,
                     label = tab$representing_key[i])
  }
  # hand values for S = {A}: eta_B = 0.9 - 0.9 * 1 = 0, eta_C = 0.5
  rowA <- tab[tab$representing_key == "A", ]
  expect_equal(unname(rowA$conditional_effect[[1]]), c(0, 0.5),
               tolerance = 1e-12)
  expect_identical(names(rowA$conditional_effect[[1]]), c("B", "C"))
  # boundary sets have NULL conditional_effect
  expect_null(tab$conditional_effect[
    tab$representing_key == "<empty>"][[1L]])
  expect_null(tab$conditional_effect[
    tab$representing_key == "A|B|C"][[1L]])
})
