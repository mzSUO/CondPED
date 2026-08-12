# Tests for attribute_traits() (Stage 6B-5): within-signal Holm
# candidate sets from final joint signal effects. No omnibus re-gate.

toy_signal_effects <- function() {
  # two signals in one locus + one SNP-level (degenerate) signal
  data.frame(
    locus_id = c(rep("chr1:100-200", 8), rep("M9", 4)),
    signal_id = c(rep("chr1:100-200::S1", 4),
                  rep("chr1:100-200::S2", 4), rep("M9::S1", 4)),
    representative_snp = c(rep("M1", 4), rep("M2", 4), rep("M9", 4)),
    trait = rep(paste0("Trait", 1:4), 3),
    beta = c(0.8, 0.5, -0.4, 0.2,   0.9, 0.3, 0.1, 0.1,  0.6, 0.5, 0.4, 0.3),
    se = 0.2,
    p_value = c(0.001, 0.010, 0.030, 0.200,   # S1: Holm -> first two
                0.040, 0.040, 0.040, 0.040,   # S2: none passes
                0.001, 0.001, 0.001, 0.001),  # M9: all pass
    stringsAsFactors = FALSE
  )
}

test_that("within-signal Holm matches hand calculation", {
  res <- attribute_traits(toy_signal_effects())
  tab <- res$trait_table
  s1 <- tab[tab$signal_id == "chr1:100-200::S1", ]
  # Holm on (0.001, 0.010, 0.030, 0.200): (0.004, 0.030, 0.060, 0.200)
  expect_equal(s1$p_adjusted, c(0.004, 0.030, 0.060, 0.200))
  expect_identical(s1$candidate, c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(res$candidate_sets[["chr1:100-200::S1"]],
                   c("Trait1", "Trait2"))
  expect_true(all(tab$adjust_scope == "within_signal"))
})

test_that("Holm scope is within_signal, never pooled across signals", {
  res <- attribute_traits(toy_signal_effects())
  tab <- res$trait_table
  for (sid in unique(tab$signal_id)) {
    rows <- tab[tab$signal_id == sid, ]
    expect_equal(rows$p_adjusted,
                 stats::p.adjust(rows$p_raw, "holm"), label = sid)
  }
  # pooled across all 12 rows would give different values
  pooled <- stats::p.adjust(tab$p_raw, "holm")
  expect_false(isTRUE(all.equal(tab$p_adjusted, pooled)))
  # S2 (secondary signal, could be marginally non-significant) still
  # goes through Holm and correctly gets no candidate
  expect_false("chr1:100-200::S2" %in% names(res$candidate_sets))
  expect_identical(res$signal_table$signal_status[
    res$signal_table$signal_id == "chr1:100-200::S2"], "omnibus_only")
})

test_that("signal_table: breadth, direction patterns, statuses", {
  res <- attribute_traits(toy_signal_effects())
  st <- res$signal_table
  expect_identical(st$effect_breadth, st$n_candidates)
  expect_identical(st$effect_breadth, c(2L, 0L, 4L))
  expect_identical(st$direction_pattern,
                   c("concordant", "not_applicable", "concordant"))
  expect_identical(st$signal_status,
                   c("multi_trait", "omnibus_only", "multi_trait"))

  # antagonistic: B = 2 with opposite signs
  eff <- toy_signal_effects()
  eff <- eff[eff$signal_id == "chr1:100-200::S1", ]
  eff$beta[2] <- -abs(eff$beta[2])
  r2 <- attribute_traits(eff)
  expect_identical(r2$signal_table$direction_pattern, "antagonistic")
  # single_trait
  eff2 <- eff
  eff2$p_value[2] <- 0.5
  r3 <- attribute_traits(eff2)
  expect_identical(r3$signal_table$direction_pattern, "single_trait")
  expect_identical(r3$signal_table$signal_status, "trait_restricted")
  # mixed: B >= 3 with both signs
  eff3 <- toy_signal_effects()
  eff3 <- eff3[eff3$signal_id == "M9::S1", ]
  eff3$beta[2] <- -abs(eff3$beta[2])
  r4 <- attribute_traits(eff3)
  expect_identical(r4$signal_table$direction_pattern, "mixed")
})

test_that("predefined mode uses the given sets; all_traits is not an oracle", {
  eff <- toy_signal_effects()
  sets <- list("chr1:100-200::S1" = c("Trait1", "Trait3"),
               "chr1:100-200::S2" = "Trait2",
               "M9::S1" = character())
  res <- attribute_traits(eff, candidate_mode = "predefined",
                          predefined_sets = sets)
  expect_identical(res$candidate_sets[["chr1:100-200::S1"]],
                   c("Trait1", "Trait3"))
  expect_identical(res$candidate_sets[["chr1:100-200::S2"]], "Trait2")
  expect_false("M9::S1" %in% names(res$candidate_sets))
  expect_true(all(res$trait_table$adjust_scope == "predefined"))
  # missing set -> error; unknown trait -> error
  expect_error(attribute_traits(eff, candidate_mode = "predefined",
                                predefined_sets = sets[1:2]),
               class = "condped_invalid_input")
  bad <- sets
  bad[["M9::S1"]] <- "TraitZ"
  expect_error(attribute_traits(eff, candidate_mode = "predefined",
                                predefined_sets = bad),
               class = "condped_invalid_input")
  # all_traits takes everything and differs from the oracle set
  res_all <- attribute_traits(eff, candidate_mode = "all_traits")
  expect_identical(res_all$candidate_sets[["chr1:100-200::S1"]],
                   paste0("Trait", 1:4))
  expect_false(identical(res_all$candidate_sets[["chr1:100-200::S1"]],
                         sets[["chr1:100-200::S1"]]))
})

test_that("SNP-level inputs map to degenerate one-signal loci", {
  df <- data.frame(
    marker_id = rep(c("M1", "M2"), each = 3),
    trait = rep(c("T1", "T2", "T3"), 2),
    beta = 1, se = 0.2,
    p_value = c(0.001, 0.4, 0.5, 0.6, 0.7, 0.8),
    stringsAsFactors = FALSE
  )
  res <- attribute_traits(df)
  expect_true(all(res$trait_table$locus_id == res$trait_table$representative_snp))
  expect_identical(unique(res$trait_table$signal_id),
                   c("M1::S1", "M2::S1"))
  expect_identical(res$candidate_sets[["M1::S1"]], "T1")
  expect_identical(res$signal_table$signal_status,
                   c("trait_restricted", "omnibus_only"))
})

test_that("resolve_locus_signals output feeds attribute directly", {
  fake_resolve <- list(
    beta = data.frame(
      locus_id = "chr1:1-9", signal_order = 1:2,
      representative_snp = c("S1", "S2"),
      trait = rep(c("T1", "T2"), each = 2),
      # rows: (S1,T1) (S2,T1) (S1,T2) (S2,T2) -> S1 strong, S2 weak
      beta = c(1, 0.1, 1, 0.1), se = 0.2,
      stringsAsFactors = FALSE
    ),
    status = list(ok = TRUE)
  )
  res <- attribute_traits(fake_resolve)
  expect_identical(unique(res$trait_table$signal_id),
                   c("chr1:1-9::S1", "chr1:1-9::S2"))
  expect_identical(res$candidate_sets[["chr1:1-9::S1"]], c("T1", "T2"))
})

test_that("empty input, NA p-values, return_all and validation", {
  empty <- data.frame(
    locus_id = character(), signal_id = character(),
    representative_snp = character(), trait = character(),
    beta = numeric(), se = numeric(), p_value = numeric(),
    stringsAsFactors = FALSE
  )
  res <- attribute_traits(empty)
  expect_identical(res$status$code, "empty_selection")
  expect_true(res$status$ok)
  expect_identical(nrow(res$signal_table), 0L)

  eff <- toy_signal_effects()
  eff$p_value[eff$signal_id == "chr1:100-200::S1" &
                eff$trait == "Trait2"] <- NA
  res2 <- attribute_traits(eff)
  row_na <- res2$trait_table$signal_id == "chr1:100-200::S1" &
    res2$trait_table$trait == "Trait2"
  expect_false(res2$trait_table$candidate[row_na])
  fam <- res2$trait_table[res2$trait_table$signal_id ==
                            "chr1:100-200::S1" & !row_na, ]
  expect_equal(fam$p_adjusted,
               stats::p.adjust(c(0.001, 0.030, 0.200), "holm"))

  slim <- attribute_traits(toy_signal_effects(), return_all = FALSE)
  expect_true(all(slim$trait_table$candidate))
  expect_identical(slim$candidate_sets,
                   attribute_traits(toy_signal_effects())$candidate_sets)

  expect_error(attribute_traits(toy_signal_effects(), alpha_trait = 2),
               class = "condped_invalid_input")
  expect_error(attribute_traits(list(status = list(ok = FALSE))),
               class = "condped_invalid_input")
  dup <- rbind(toy_signal_effects(), toy_signal_effects()[1, ])
  expect_error(attribute_traits(dup), class = "condped_invalid_input")
})
