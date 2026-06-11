# tests/testthat/test-classify.R
# Test suite for six-class classification system

# -----------------------------------------------------------------------------
# classify_locus
# -----------------------------------------------------------------------------

test_that("null: zero marginal significance", {
  result <- classify_locus(n_marg_sig = 0, cond_sig_by_trait = c(A = TRUE, B = FALSE))
  expect_equal(result, "null")
})

test_that("class1: single trait significance", {
  result <- classify_locus(n_marg_sig = 1, cond_sig_by_trait = c(A = TRUE, B = FALSE))
  expect_equal(result, "class1")
})

test_that("class3: multi-trait, multiple conditional sig, no MR", {
  result <- classify_locus(
    n_marg_sig = 2,
    cond_sig_by_trait = c(A = TRUE, B = TRUE),
    mr_results = NULL
  )
  expect_equal(result, "class3")
})

test_that("class2: multi-trait, few conditional sig, no MR", {
  result <- classify_locus(
    n_marg_sig = 2,
    cond_sig_by_trait = c(A = TRUE, B = FALSE),
    mr_results = NULL
  )
  expect_equal(result, "class2")
})

test_that("class6: bidirectional MR significant", {
  mr <- list(AB = list(sig = TRUE), BA = list(sig = TRUE))
  result <- classify_locus(
    n_marg_sig = 2,
    cond_sig_by_trait = c(A = TRUE, B = TRUE),
    mr_results = mr,
    traits = c("A", "B")
  )
  expect_equal(result, "class6")
})

test_that("class4: complete mediation (A->B, B condition not sig)", {
  # MR A->B significant, B's conditional effect not significant
  mr <- list(AB = list(sig = TRUE), BA = list(sig = FALSE))
  result <- classify_locus(
    n_marg_sig = 2,
    cond_sig_by_trait = c(A = TRUE, B = FALSE),  # B not significant after conditioning
    mr_results = mr,
    traits = c("A", "B")
  )
  expect_equal(result, "class4")
})

test_that("class5: partial mediation (A->B, B condition still sig)", {
  # MR A->B significant, B's conditional effect still significant
  mr <- list(AB = list(sig = TRUE), BA = list(sig = FALSE))
  result <- classify_locus(
    n_marg_sig = 2,
    cond_sig_by_trait = c(A = TRUE, B = TRUE),  # B still significant
    mr_results = mr,
    traits = c("A", "B")
  )
  expect_equal(result, "class5")
})

test_that("class2: no MR direction, few conditional sig", {
  mr <- list(AB = list(sig = FALSE), BA = list(sig = FALSE))
  result <- classify_locus(
    n_marg_sig = 2,
    cond_sig_by_trait = c(A = TRUE, B = FALSE),
    mr_results = mr,
    traits = c("A", "B")
  )
  expect_equal(result, "class2")
})

test_that("class3: no MR direction, multiple conditional sig", {
  mr <- list(AB = list(sig = FALSE), BA = list(sig = FALSE))
  result <- classify_locus(
    n_marg_sig = 2,
    cond_sig_by_trait = c(A = TRUE, B = TRUE),
    mr_results = mr,
    traits = c("A", "B")
  )
  expect_equal(result, "class3")
})

test_that("class4 with BA direction significant", {
  # MR B->A significant, A's conditional not significant
  mr <- list(AB = list(sig = FALSE), BA = list(sig = TRUE))
  result <- classify_locus(
    n_marg_sig = 2,
    cond_sig_by_trait = c(A = FALSE, B = TRUE),
    mr_results = mr,
    traits = c("A", "B")
  )
  expect_equal(result, "class4")  # A is outcome, A not significant after conditioning
})

test_that("handles NA in n_marg_sig", {
  result <- classify_locus(n_marg_sig = NA, cond_sig_by_trait = c(A = TRUE))
  expect_equal(result, "null")
})


# -----------------------------------------------------------------------------
# summarise_classification
# -----------------------------------------------------------------------------

test_that("summarise_classification counts correctly", {
  result <- data.frame(
    snp_id = c("S1", "S2", "S3", "S4", "S5"),
    class = c("class1", "class1", "class3", "class4", "class4"),
    stringsAsFactors = FALSE
  )

  summary <- summarise_classification(result, verbose = FALSE)
  expect_equal(nrow(summary), 3)
  expect_equal(summary$n[summary$class == "class1"], 2)
  expect_equal(summary$n[summary$class == "class4"], 2)
  expect_equal(summary$pct[summary$class == "class3"], 20.0)
})

test_that("summarise_classification handles all classes", {
  result <- data.frame(
    snp_id = paste0("S", 1:6),
    class = c("class1", "class2", "class3", "class4", "class5", "class6"),
    stringsAsFactors = FALSE
  )

  summary <- summarise_classification(result, verbose = FALSE)
  expect_equal(nrow(summary), 6)
  expect_true(all(summary$n == 1))
})

test_that("summarise_classification errors without class column", {
  result <- data.frame(snp_id = "S1", foo = "bar")
  expect_error(summarise_classification(result))
})

test_that("summarise_classification handles empty data", {
  result <- data.frame(snp_id = character(), class = character(), stringsAsFactors = FALSE)
  summary <- summarise_classification(result, verbose = FALSE)
  expect_equal(nrow(summary), 0)
})