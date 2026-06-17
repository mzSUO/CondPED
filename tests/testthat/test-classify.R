library(testthat)

traits <- c("TraitA", "TraitB")

# ── 1. 基础路径 ─────────────────────────────────────────────────────────────

test_that("null: marginal non-significant", {
  expect_equal(classify_locus(0L, c(TraitA = TRUE, TraitB = TRUE)), "null")
  expect_equal(classify_locus(NA_integer_, c(TraitA = TRUE, TraitB = TRUE)), "null")
})

test_that("class1: single trait marginal significant", {
  expect_equal(classify_locus(1L, c(TraitA = TRUE, TraitB = FALSE)), "class1")
  expect_equal(classify_locus(1L, c(TraitA = FALSE, TraitB = TRUE)), "class1")
})

# ── 2. 消融实验（mr_results = NULL）────────────────────────────────────────

test_that("ablation: n_cond >= 2 → class2 (horizontal)", {
  expect_equal(
    classify_locus(2L, c(TraitA = TRUE, TraitB = TRUE), mr_results = NULL),
    "class2"
  )
})

test_that("ablation: n_cond < 2 → class1 (conservative downgrade)", {
  # 完全中介指纹：n_cond=1，无 MR 时无法区分，保守降级为 class1
  expect_equal(
    classify_locus(2L, c(TraitA = TRUE, TraitB = FALSE), mr_results = NULL),
    "class1"
  )
})

# ── 3. 双向 MR ─────────────────────────────────────────────────────────────

test_that("class5: bidirectional MR significant", {
  mr <- list(AB = list(sig = TRUE), BA = list(sig = TRUE))
  expect_equal(classify_locus(2L, c(TraitA = TRUE, TraitB = TRUE), mr, traits), "class5")
})

# ── 4. MR 均不显著 ─────────────────────────────────────────────────────────

test_that("class2: no MR + n_cond >= 2 → horizontal pleiotropy", {
  mr <- list(AB = list(sig = FALSE), BA = list(sig = FALSE))
  expect_equal(classify_locus(2L, c(TraitA = TRUE, TraitB = TRUE), mr, traits), "class2")
})

test_that("class1: no MR + n_cond < 2 → conservative", {
  mr <- list(AB = list(sig = FALSE), BA = list(sig = FALSE))
  expect_equal(classify_locus(2L, c(TraitA = TRUE, TraitB = FALSE), mr, traits), "class1")
})

# ── 5. 单向 MR：A→B 显著 ───────────────────────────────────────────────────

test_that("class3: A→B MR + outcome(B) NOT sig → complete mediation", {
  # SNP → A → B，条件投影后 B 被解释掉
  mr <- list(AB = list(sig = TRUE), BA = list(sig = FALSE))
  expect_equal(classify_locus(2L, c(TraitA = TRUE, TraitB = FALSE), mr, traits), "class3")
})

test_that("class4: A→B MR + outcome(B) sig → partial mediation", {
  # SNP → A, SNP → B, A → B，条件投影后 B 仍有直接效应
  mr <- list(AB = list(sig = TRUE), BA = list(sig = FALSE))
  expect_equal(classify_locus(2L, c(TraitA = TRUE, TraitB = TRUE), mr, traits), "class4")
})

# ── 6. 单向 MR：B→A 显著 ───────────────────────────────────────────────────

test_that("class3: B→A MR + outcome(A) NOT sig → complete mediation", {
  mr <- list(AB = list(sig = FALSE), BA = list(sig = TRUE))
  expect_equal(classify_locus(2L, c(TraitA = FALSE, TraitB = TRUE), mr, traits), "class3")
})

test_that("class4: B→A MR + outcome(A) sig → partial mediation", {
  mr <- list(AB = list(sig = FALSE), BA = list(sig = TRUE))
  expect_equal(classify_locus(2L, c(TraitA = TRUE, TraitB = TRUE), mr, traits), "class4")
})

# ── 7. 边界防御 ─────────────────────────────────────────────────────────────

test_that("defensive: n_cond=0 with unidirectional MR → class1", {
  # 理论上不应发生：投影后都不显著，但 MR 显著
  mr <- list(AB = list(sig = TRUE), BA = list(sig = FALSE))
  expect_equal(classify_locus(2L, c(TraitA = FALSE, TraitB = FALSE), mr, traits), "class1")
})

test_that("defensive: input validation", {
  expect_error(classify_locus(2L, c(TRUE, FALSE), traits = NULL),
               "traits 与 cond_sig_by_trait 名称不能同时为空")
  expect_error(classify_locus(2L, c(TraitA = TRUE, TraitB = FALSE), traits = c("A")),
               "traits 必须是长度为 2 的字符向量")
})