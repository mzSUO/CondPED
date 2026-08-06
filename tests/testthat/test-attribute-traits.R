# Tests for attribute_traits() (v1.0 Stage 3): candidate trait sets via
# within-locus Holm. Replaces the v0.3 BB/pooled-Holm tests.

# ---- shared toy fixtures -----------------------------------------------------

# Four omnibus families; BH at 0.05 selects exactly L1, L2.
toy_omnibus <- function() {
  data.frame(
    marker_id = paste0("L", 1:4),
    p_value = c(0.001, 0.004, 0.20, 0.60),
    stringsAsFactors = FALSE
  )
}

# Layer-2 anchor (Runbook 5.4): four raw p-values per locus.
toy_effects <- function() {
  data.frame(
    marker_id = rep(paste0("L", 1:4), each = 4),
    trait = rep(paste0("Trait", 1:4), times = 4),
    beta = c(0.8, 0.5, -0.4, 0.2,   0.9, 0.3, 0.1, 0.1,
             0.5, 0.5, 0.5, 0.5,   0.5, 0.5, 0.5, 0.5),
    se = 0.2,
    p_value = c(0.001, 0.010, 0.030, 0.200,   # L1 anchor
                0.040, 0.040, 0.040, 0.040,   # L2: none passes Holm
                0.001, 0.001, 0.001, 0.001,   # L3 (never selected)
                0.001, 0.001, 0.001, 0.001),  # L4 (never selected)
    stringsAsFactors = FALSE
  )
}

# ---- Runbook 5.4 hand-calculation anchor -------------------------------------

test_that("within-locus Holm matches the hand calculation", {
  res <- attribute_traits(toy_omnibus(), toy_effects(),
                          omnibus_method = "BH",
                          candidate_mode = "holm_fwer", alpha_trait = 0.05)
  tab <- res$trait_table[res$trait_table$marker_id == "L1", ]
  # Holm on (0.001, 0.010, 0.030, 0.200):
  #   sorted x (4,3,2,1) = (0.004, 0.030, 0.060, 0.200), cummax same
  #   -> first two <= 0.05 are candidates
  expect_equal(tab$p_adjusted, c(0.004, 0.030, 0.060, 0.200))
  expect_identical(tab$candidate, c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(res$candidate_sets$L1, c("Trait1", "Trait2"))
  expect_identical(tab$adjust_scope, rep("within_locus", 4L))
})

test_that("adjustment scope is within_locus, not pooled across loci", {
  res <- attribute_traits(toy_omnibus(), toy_effects(),
                          omnibus_method = "BH",
                          candidate_mode = "holm_fwer", alpha_trait = 0.05)
  tab <- res$trait_table
  # per-locus adjusted values equal p.adjust on that locus alone
  p1 <- c(0.001, 0.010, 0.030, 0.200)
  p2 <- c(0.040, 0.040, 0.040, 0.040)
  expect_equal(tab$p_adjusted[tab$marker_id == "L1"],
               stats::p.adjust(p1, "holm"))
  expect_equal(tab$p_adjusted[tab$marker_id == "L2"],
               stats::p.adjust(p2, "holm"))
  # and differ from pooling all eight p-values
  pooled <- stats::p.adjust(c(p1, p2), "holm")
  expect_false(isTRUE(all.equal(tab$p_adjusted[tab$marker_id == "L2"],
                                pooled[5:8])))
  # L2: 0.04 * 4 = 0.16 > 0.05 -> omnibus_only
  expect_identical(res$locus_table$locus_status[res$locus_table$marker_id == "L2"],
                   "omnibus_only")
  expect_false("L2" %in% names(res$candidate_sets))
})

# ---- layer 1 ------------------------------------------------------------------

test_that("layer 1 matches p.adjust and M counts only valid families", {
  om <- data.frame(marker_id = paste0("L", 1:6),
                   p_value = c(0.001, 0.004, 0.20, 0.60, NA, NA))
  eff <- toy_effects()
  eff <- rbind(eff, data.frame(
    marker_id = rep(c("L5", "L6"), each = 4),
    trait = rep(paste0("Trait", 1:4), 2),
    beta = 0, se = 1, p_value = NA, stringsAsFactors = FALSE
  ))
  res <- attribute_traits(om, eff, omnibus_method = "BH")
  expect_identical(res$diagnostics$M, 4L)
  expect_identical(res$diagnostics$n_na_omnibus_p, 2L)
  expect_identical(res$selected_loci, c("L1", "L2"))
  # unselected loci never enter layer 2
  expect_false(any(res$trait_table$marker_id %in% c("L3", "L4", "L5", "L6")))
  expect_false(any(names(res$candidate_sets) %in% c("L3", "L4", "L5", "L6")))
  # bonferroni / none agree with the direct computation
  for (method in c("bonferroni", "none")) {
    r <- attribute_traits(om[, ], eff, omnibus_method = method,
                          alpha_omnibus = 0.05)
    p <- om$p_value[!is.na(om$p_value)]
    adj <- if (method == "none") p else stats::p.adjust(p, method)
    expect_identical(r$selected_loci,
                     om$marker_id[!is.na(om$p_value)][adj <= 0.05],
                     label = method)
  }
})

# ---- locus status ---------------------------------------------------------------

test_that("locus_status covers omnibus_only / trait_restricted / multi_trait", {
  eff <- data.frame(
    marker_id = rep(c("L1", "L2", "L3"), each = 3),
    trait = rep(c("T1", "T2", "T3"), 3),
    beta = 1, se = 0.2,
    p_value = c(0.001, 0.002, 0.010,   # L1: all pass -> multi_trait
                0.001, 0.500, 0.500,   # L2: one passes -> trait_restricted
                0.100, 0.200, 0.300),  # L3: none passes -> omnibus_only
    stringsAsFactors = FALSE
  )
  om <- data.frame(marker_id = c("L1", "L2", "L3"),
                   p_value = c(0.001, 0.001, 0.001))
  res <- attribute_traits(om, eff, omnibus_method = "none",
                          alpha_omnibus = 0.05, alpha_trait = 0.05)
  lt <- res$locus_table
  expect_identical(lt$locus_status[lt$marker_id == "L1"], "multi_trait")
  expect_identical(lt$locus_status[lt$marker_id == "L2"], "trait_restricted")
  expect_identical(lt$locus_status[lt$marker_id == "L3"], "omnibus_only")
  expect_identical(lt$n_candidates, c(3L, 1L, 0L))
  expect_identical(res$diagnostics$n_multi_trait, 1L)
  expect_identical(res$diagnostics$n_trait_restricted, 1L)
  expect_identical(res$diagnostics$n_omnibus_only, 1L)
  # omnibus_only locus has no candidate set
  expect_identical(res$candidate_sets$L3, NULL)
  expect_identical(res$candidate_sets$L2, "T1")
})

# ---- all_traits oracle mode -----------------------------------------------------

test_that("all_traits takes every trait without testing, selected loci only", {
  res <- attribute_traits(toy_omnibus(), toy_effects(),
                          omnibus_method = "BH",
                          candidate_mode = "all_traits")
  tab <- res$trait_table
  expect_true(all(tab$candidate))
  expect_identical(tab$p_adjusted, tab$p_raw)
  expect_true(all(tab$adjust_scope == "none"))
  expect_identical(res$candidate_sets$L1, paste0("Trait", 1:4))
  expect_identical(res$candidate_sets$L2, paste0("Trait", 1:4))
  # still only selected loci
  expect_false("L3" %in% names(res$candidate_sets))
  expect_true(all(res$locus_table$locus_status == "multi_trait"))
})

# ---- robustness -----------------------------------------------------------------

test_that("empty selection returns a complete empty result", {
  om <- data.frame(marker_id = paste0("L", 1:4), p_value = c(0.5, 0.6, 0.7, 0.8))
  res <- attribute_traits(om, toy_effects(), omnibus_method = "BH")
  expect_true(res$status$ok)
  expect_match(paste(res$status$warnings, collapse = " "), "no loci")
  expect_identical(res$selected_loci, character())
  expect_identical(nrow(res$trait_table), 0L)
  expect_identical(names(res$trait_table),
                   c("marker_id", "trait", "beta", "se", "Q", "p_raw",
                     "p_adjusted", "candidate", "adjust_scope"))
  expect_identical(res$candidate_sets, setNames(list(), character()))
  expect_identical(nrow(res$locus_table), 0L)
  expect_identical(res$diagnostics$n_selected, 0L)
})

test_that("NA p-values, missing loci and duplicate keys are handled", {
  eff <- toy_effects()
  eff$p_value[eff$marker_id == "L1" & eff$trait == "Trait2"] <- NA
  res <- attribute_traits(toy_omnibus(), eff, omnibus_method = "BH")
  row_na <- res$trait_table$marker_id == "L1" &
    res$trait_table$trait == "Trait2"
  expect_true(is.na(res$trait_table$p_adjusted[row_na]))
  expect_false(res$trait_table$candidate[row_na])
  expect_identical(res$diagnostics$n_na_effect_p, 1L)
  # remaining L1 members Holm-adjusted among themselves
  fam <- res$trait_table[res$trait_table$marker_id == "L1" & !row_na, ]
  expect_equal(fam$p_adjusted,
               stats::p.adjust(c(0.001, 0.030, 0.200), "holm"))

  # selected locus missing from effects
  eff2 <- toy_effects()[toy_effects()$marker_id != "L2", ]
  res2 <- attribute_traits(toy_omnibus(), eff2, omnibus_method = "BH")
  expect_match(paste(res2$status$warnings, collapse = " "), "no rows in effects")
  expect_identical(res2$diagnostics$n_missing_effect_loci, 1L)
  expect_false("L2" %in% res2$trait_table$marker_id)

  # duplicate keys rejected
  expect_error(attribute_traits(toy_omnibus(),
                                rbind(toy_effects(), toy_effects()[1, ])),
               class = "condped_invalid_input")
})

test_that("result objects accepted; failed statuses rejected; return_all", {
  om <- toy_omnibus(); eff <- toy_effects()
  as_obj <- attribute_traits(
    list(omnibus = om, status = list(ok = TRUE)),
    list(effects_long = eff, status = list(ok = TRUE)),
    omnibus_method = "BH"
  )
  as_df <- attribute_traits(om, eff, omnibus_method = "BH")
  expect_identical(as_obj$trait_table, as_df$trait_table)
  expect_error(attribute_traits(list(omnibus = om, status = list(ok = FALSE)),
                                eff),
               class = "condped_invalid_input")
  slim <- attribute_traits(om, eff, omnibus_method = "BH",
                           return_all = FALSE)
  expect_true(all(slim$trait_table$candidate))
  expect_identical(slim$candidate_sets, as_df$candidate_sets)
  expect_identical(slim$locus_table, as_df$locus_table)
})

test_that("Q equals the squared Wald statistic; trait order preserved", {
  res <- attribute_traits(toy_omnibus(), toy_effects(), omnibus_method = "BH")
  tab <- res$trait_table
  expect_equal(tab$Q, (tab$beta / tab$se)^2)
  # rows keep the effects order within each locus
  expect_identical(tab$trait[tab$marker_id == "L1"],
                   paste0("Trait", 1:4))
  expect_identical(names(res$candidate_sets), "L1")
})

test_that("invalid arguments are rejected", {
  expect_error(attribute_traits(toy_omnibus(), toy_effects(),
                                alpha_omnibus = 0),
               class = "condped_invalid_input")
  expect_error(attribute_traits(toy_omnibus(), toy_effects(),
                                alpha_trait = 1),
               class = "condped_invalid_input")
  expect_error(attribute_traits(toy_omnibus(), toy_effects(),
                                candidate_mode = "bb_fdr"),
               class = "error")
  expect_error(attribute_traits(toy_omnibus(), toy_effects(),
                                return_all = NA),
               class = "condped_invalid_input")
})
