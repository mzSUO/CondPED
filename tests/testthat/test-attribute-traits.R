# Tests for attribute_traits() (S4): hierarchical trait attribution.
#
# Covers, per the S4 task specification:
#   A. hand-computed Benjamini-Bogomolov on four families;
#   B. first-layer selection under BH / bonferroni / none and the meaning
#      of M;
#   C. holm_fwer pooled Holm against p.adjust(..., "holm");
#   D. none mode leaves p-values untouched;
#   E. empty selection and malformed inputs;
#   F. A is exactly the attributed rows of trait_table.

# ---- shared toy fixtures -----------------------------------------------------

# Four omnibus families; BH at 0.05 selects exactly L1, L2.
toy_omnibus <- function() {
  data.frame(
    marker_id = paste0("L", 1:4),
    p_value = c(0.001, 0.004, 0.20, 0.60),
    stringsAsFactors = FALSE
  )
}

toy_effects <- function() {
  data.frame(
    marker_id = rep(paste0("L", 1:4), each = 3),
    trait = rep(paste0("Trait", 1:3), times = 4),
    beta = c(0.36, 0.15, -0.67, 0.84, -0.51, -0.16,
             0.50, 0.50, 0.50, 0.50, 0.50, 0.50),
    se = rep(0.1, 12),
    p_value = c(0.001, 0.010, 0.040,   # L1 family
                0.020, 0.030, 0.500,   # L2 family
                0.001, 0.001, 0.001,   # L3 (never selected)
                0.001, 0.001, 0.001),  # L4 (never selected)
    stringsAsFactors = FALSE
  )
}

# ---- A. hand-computed BB -----------------------------------------------------

test_that("BB matches a four-family hand calculation", {
  res <- attribute_traits(
    toy_omnibus(), toy_effects(),
    omnibus_method = "BH", attribution_mode = "bb_fdr", q_target = 0.05
  )
  # Layer 1: BH over M = 4 families
  expect_identical(res$diagnostics$M, 4L)
  expect_identical(res$selected_loci, c("L1", "L2"))
  # q_within = q_target * |S| / M = 0.05 * 2 / 4
  expect_equal(res$diagnostics$q_within, 0.025)

  # Hand BH within family L1: p = (0.001, 0.010, 0.040)
  #   adjusted = (0.003, 0.015, 0.040) -> reject first two at 0.025
  # Hand BH within family L2: p = (0.020, 0.030, 0.500)
  #   adjusted = (0.045, 0.045, 0.500) -> reject none at 0.025
  tab <- res$trait_table
  expect_equal(tab$p_adjusted[tab$marker_id == "L1"],
               c(0.003, 0.015, 0.040))
  expect_equal(tab$p_adjusted[tab$marker_id == "L2"],
               c(0.045, 0.045, 0.500))
  expect_identical(tab$attributed[tab$marker_id == "L1"],
                   c(TRUE, TRUE, FALSE))
  expect_identical(tab$attributed[tab$marker_id == "L2"],
                   c(FALSE, FALSE, FALSE))
  expect_identical(res$A, list(L1 = c("Trait1", "Trait2")))
  expect_identical(res$diagnostics$n_attributed, 2L)
})

# ---- B. first layer ----------------------------------------------------------

test_that("first layer matches p.adjust for BH and bonferroni", {
  om <- data.frame(marker_id = paste0("L", 1:5),
                   p_value = c(0.001, 0.010, 0.020, 0.040, 0.060))
  eff <- data.frame(
    marker_id = rep(om$marker_id, each = 2),
    trait = rep(c("T1", "T2"), 5),
    beta = 0, se = 1, p_value = 0.9, stringsAsFactors = FALSE
  )
  for (method in c("BH", "bonferroni")) {
    res <- attribute_traits(om, eff, omnibus_method = method,
                            attribution_mode = "none", q_target = 0.5)
    expected <- om$marker_id[
      stats::p.adjust(om$p_value, method = method) <= 0.05
    ]
    expect_identical(res$selected_loci, expected, label = method)
  }
  # none: raw threshold
  res <- attribute_traits(om, eff, omnibus_method = "none",
                          attribution_mode = "none", q_target = 0.5)
  expect_identical(res$selected_loci, om$marker_id[om$p_value <= 0.05])
})

test_that("M counts only valid (non-NA) omnibus families", {
  om <- data.frame(marker_id = paste0("L", 1:6),
                   p_value = c(0.001, 0.004, 0.20, 0.60, NA, NA))
  eff <- toy_effects()
  eff <- rbind(eff, data.frame(
    marker_id = rep(c("L5", "L6"), each = 3),
    trait = rep(paste0("Trait", 1:3), 2),
    beta = 0, se = 1, p_value = NA, stringsAsFactors = FALSE
  ))
  res <- attribute_traits(om, eff, omnibus_method = "BH",
                          attribution_mode = "bb_fdr", q_target = 0.05)
  expect_identical(res$diagnostics$M, 4L)          # NA rows not counted
  expect_identical(res$diagnostics$n_omnibus_rows, 6L)
  expect_identical(res$diagnostics$n_na_omnibus_p, 2L)
  expect_identical(res$diagnostics$n_selected, 2L)
  expect_equal(res$diagnostics$q_within, 0.05 * 2 / 4)
  # unselected loci never reach layer 2
  expect_false(any(res$trait_table$marker_id %in% c("L3", "L4", "L5", "L6")))
  expect_false(any(names(res$A) %in% c("L3", "L4", "L5", "L6")))
})

# ---- C. holm_fwer ------------------------------------------------------------

test_that("holm_fwer pools all tested hypotheses and matches p.adjust holm", {
  res <- attribute_traits(
    toy_omnibus(), toy_effects(),
    omnibus_method = "bonferroni", attribution_mode = "holm_fwer",
    alpha_total = 0.05, alpha_split = c(0.025, 0.025)
  )
  # Layer 1 budget: alpha_split[1] = 0.025, bonferroni over M = 4
  #   adjusted = (0.004, 0.016, 0.80, 1) -> L1, L2 selected
  expect_identical(res$selected_loci, c("L1", "L2"))
  expect_equal(res$diagnostics$alpha_layer1, 0.025)
  expect_equal(res$diagnostics$alpha_layer2, 0.025)
  # Layer 2: pooled Holm over the 6 tested p-values at alpha_split[2]
  pool <- c(0.001, 0.010, 0.040, 0.020, 0.030, 0.500)
  expect_equal(res$trait_table$p_adjusted,
               stats::p.adjust(pool, method = "holm"))
  expect_identical(res$trait_table$attributed,
                   stats::p.adjust(pool, method = "holm") <= 0.025)
})

test_that("holm_fwer requires alpha_split to sum to alpha_total", {
  expect_error(
    attribute_traits(toy_omnibus(), toy_effects(),
                     attribution_mode = "holm_fwer",
                     alpha_total = 0.05, alpha_split = c(0.01, 0.01)),
    class = "condped_invalid_input"
  )
  # in non-holm modes the sum is not required
  expect_no_error(
    attribute_traits(toy_omnibus(), toy_effects(),
                     attribution_mode = "bb_fdr",
                     alpha_split = c(0.01, 0.01))
  )
})

# ---- D. none -----------------------------------------------------------------

test_that("none mode leaves p-values unchanged and uses an explicit rule", {
  res <- attribute_traits(
    toy_omnibus(), toy_effects(),
    omnibus_method = "BH", attribution_mode = "none", q_target = 0.025
  )
  expect_identical(res$trait_table$p_adjusted, res$trait_table$p_raw)
  expect_identical(res$trait_table$attributed,
                   res$trait_table$p_raw <= 0.025)
  expect_match(paste(res$assumptions, collapse = " "),
               "p_raw <= q_target", fixed = TRUE)
})

# ---- E. empty selection and malformed input ----------------------------------

test_that("empty first-layer selection returns a complete empty result", {
  om <- data.frame(marker_id = paste0("L", 1:4), p_value = c(0.5, 0.6, 0.7, 0.8))
  res <- attribute_traits(om, toy_effects(), omnibus_method = "BH",
                          attribution_mode = "bb_fdr")
  expect_true(res$status$ok)
  expect_match(paste(res$status$warnings, collapse = " "), "no loci")
  expect_identical(res$selected_loci, character())
  expect_identical(nrow(res$trait_table), 0L)
  expect_identical(names(res$trait_table),
                   c("marker_id", "trait", "beta", "se",
                     "p_raw", "p_adjusted", "attributed"))
  expect_identical(res$A, setNames(list(), character()))
  expect_identical(res$diagnostics$n_selected, 0L)
  expect_identical(res$diagnostics$n_attributed, 0L)
  expect_equal(res$diagnostics$q_within, 0)
})

test_that("selected locus missing from effects is recorded, never enters A", {
  eff <- toy_effects()
  eff <- eff[eff$marker_id != "L2", ]   # L2 selected but absent from effects
  res <- attribute_traits(toy_omnibus(), eff, omnibus_method = "BH",
                          attribution_mode = "bb_fdr")
  expect_true(res$status$ok)
  expect_match(paste(res$status$warnings, collapse = " "), "no rows in effects")
  expect_identical(res$diagnostics$n_missing_effect_loci, 1L)
  expect_false("L2" %in% res$trait_table$marker_id)
  expect_false("L2" %in% names(res$A))
})

test_that("NA per-trait p-values are kept, unattributed and recorded", {
  eff <- toy_effects()
  eff$p_value[eff$marker_id == "L1" & eff$trait == "Trait2"] <- NA
  res <- attribute_traits(toy_omnibus(), eff, omnibus_method = "BH",
                          attribution_mode = "bb_fdr")
  row_na <- res$trait_table$marker_id == "L1" &
    res$trait_table$trait == "Trait2"
  expect_true(is.na(res$trait_table$p_adjusted[row_na]))
  expect_false(res$trait_table$attributed[row_na])
  expect_identical(res$diagnostics$n_na_effect_p, 1L)
  expect_match(paste(res$status$warnings, collapse = " "), "NA")
  # remaining family members still BH-adjusted among themselves
  fam <- res$trait_table[res$trait_table$marker_id == "L1" & !row_na, ]
  expect_equal(fam$p_adjusted, stats::p.adjust(c(0.001, 0.040), "BH"))
})

test_that("duplicated marker_id x trait keys are rejected", {
  eff <- rbind(toy_effects(), toy_effects()[1, ])
  expect_error(
    attribute_traits(toy_omnibus(), eff),
    class = "condped_invalid_input"
  )
})

test_that("trait order in effects does not change the result", {
  eff <- toy_effects()
  perm <- sample.int(nrow(eff))
  res1 <- attribute_traits(toy_omnibus(), eff, omnibus_method = "BH",
                           attribution_mode = "bb_fdr")
  res2 <- attribute_traits(toy_omnibus(), eff[perm, ], omnibus_method = "BH",
                           attribution_mode = "bb_fdr")
  tab1 <- res1$trait_table[order(res1$trait_table$marker_id,
                                 res1$trait_table$trait), ]
  tab2 <- res2$trait_table[order(res2$trait_table$marker_id,
                                 res2$trait_table$trait), ]
  expect_identical(tab1$p_adjusted, tab2$p_adjusted)
  expect_identical(tab1$attributed, tab2$attributed)
  # A follows the effects row order within a locus; compare as sets
  expect_identical(names(res1$A), names(res2$A))
  expect_identical(lapply(res1$A, sort), lapply(res2$A, sort))
})

test_that("result objects are accepted and failed statuses are rejected", {
  om <- toy_omnibus()
  eff <- toy_effects()
  as_df <- attribute_traits(om, eff, omnibus_method = "BH",
                            attribution_mode = "bb_fdr")
  as_obj <- attribute_traits(
    list(omnibus = om, status = list(ok = TRUE)),
    list(effects_long = eff, status = list(ok = TRUE)),
    omnibus_method = "BH", attribution_mode = "bb_fdr"
  )
  expect_identical(as_df$trait_table, as_obj$trait_table)
  expect_identical(as_df$A, as_obj$A)
  expect_error(
    attribute_traits(list(omnibus = om, status = list(ok = FALSE)), eff),
    class = "condped_invalid_input"
  )
  expect_error(
    attribute_traits(om, list(effects_long = eff, status = list(ok = FALSE))),
    class = "condped_invalid_input"
  )
})

test_that("invalid levels and arguments are rejected", {
  expect_error(attribute_traits(toy_omnibus(), toy_effects(),
                                alpha_omnibus = 1.5),
               class = "condped_invalid_input")
  expect_error(attribute_traits(toy_omnibus(), toy_effects(), q_target = 0),
               class = "condped_invalid_input")
  expect_error(attribute_traits(toy_omnibus(), toy_effects(),
                                alpha_split = c(0.01, 0.01, 0.01)),
               class = "condped_invalid_input")
  expect_error(attribute_traits(toy_omnibus(), toy_effects(),
                                omnibus_method = "holm"),
               class = "error")
})

test_that("dependence = empirical keeps the interface but records a warning", {
  res <- attribute_traits(toy_omnibus(), toy_effects(),
                          attribution_mode = "bb_fdr",
                          dependence = "empirical")
  expect_true(res$status$ok)
  expect_match(paste(res$status$warnings, collapse = " "),
               "empirical")
  expect_match(paste(res$assumptions, collapse = " "),
               "unverified")
  ref <- attribute_traits(toy_omnibus(), toy_effects(),
                          attribution_mode = "bb_fdr",
                          dependence = "assumed")
  expect_identical(res$trait_table, ref$trait_table)
  expect_identical(res$A, ref$A)
})

# ---- F. A consistency ---------------------------------------------------------

test_that("A is exactly trait_table's attributed rows split by marker", {
  for (mode in c("bb_fdr", "holm_fwer", "none")) {
    res <- attribute_traits(toy_omnibus(), toy_effects(),
                            omnibus_method = "BH", attribution_mode = mode)
    att <- res$trait_table[res$trait_table$attributed %in% TRUE, ]
    expected_A <- lapply(split(att$trait, att$marker_id), unname)
    expect_identical(res$A, expected_A, label = mode)
  }
})

test_that("return_all = FALSE filters rows but never changes A", {
  full <- attribute_traits(toy_omnibus(), toy_effects(),
                           omnibus_method = "BH", attribution_mode = "bb_fdr",
                           return_all = TRUE)
  slim <- attribute_traits(toy_omnibus(), toy_effects(),
                           omnibus_method = "BH", attribution_mode = "bb_fdr",
                           return_all = FALSE)
  expect_identical(full$A, slim$A)
  expect_identical(full$selected_loci, slim$selected_loci)
  expect_true(all(slim$trait_table$attributed))
  expect_identical(nrow(slim$trait_table), sum(full$trait_table$attributed))
})
