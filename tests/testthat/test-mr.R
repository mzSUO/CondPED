# tests/testthat/test-mr.R
# Test suite for Layer 3 Mendelian Randomization

# -----------------------------------------------------------------------------
# compute_f_stat
# -----------------------------------------------------------------------------

test_that("F-statistic is (A/SE)^2", {
  expect_equal(compute_f_stat(2.0, 0.5), 16.0)
  expect_equal(compute_f_stat(1.0, 1.0), 1.0)
  expect_equal(compute_f_stat(0.0, 0.5), 0.0)
})

test_that("F-statistic handles character input", {
  expect_equal(compute_f_stat("2.0", "0.5"), 16.0)
})


# -----------------------------------------------------------------------------
# ld_prune_ivs
# -----------------------------------------------------------------------------

test_that("LD pruning keeps independent SNPs", {
  ld_mat <- diag(3)
  rownames(ld_mat) <- colnames(ld_mat) <- c("S1", "S2", "S3")
  result <- ld_prune_ivs(c("S1", "S2", "S3"), ld_mat, r2_threshold = 0.1)
  expect_equal(sort(result), c("S1", "S2", "S3"))
})

test_that("LD pruning removes correlated SNPs", {
  ld_mat <- matrix(c(1, 0.9, 0, 0.9, 1, 0, 0, 0, 1), 3, 3)
  rownames(ld_mat) <- colnames(ld_mat) <- c("S1", "S2", "S3")
  result <- ld_prune_ivs(c("S1", "S2", "S3"), ld_mat, r2_threshold = 0.1)
  # S1 and S2 are correlated (r=0.9), should keep only one
  expect_true(length(result) >= 2)
  expect_false(all(c("S1", "S2") %in% result))
})

test_that("LD pruning handles single SNP", {
  expect_equal(ld_prune_ivs("S1", NULL, 0.1), "S1")
})


# -----------------------------------------------------------------------------
# select_instruments
# -----------------------------------------------------------------------------

test_that("instrument selection applies all four criteria", {
  marginal <- data.frame(
    TRAIT = c("A", "A", "A", "B", "B"),
    QTL = c("S1", "S2", "S3", "S1", "S4"),
    A = c("1.0", "0.5", "0.3", "0.8", "0.2"),
    SE = c("0.1", "0.1", "0.1", "0.1", "0.1"),
    P_Value = c("1e-10", "1e-8", "1e-5", "1e-9", "1e-3"),
    stringsAsFactors = FALSE
  )

  cond <- data.frame(
    trait = c("B", "B", "B"),
    locus = c("S1", "S2", "S3"),
    pval_cond = c(0.5, 0.1, 0.01),  # S1 and S2 pass exclusivity (p > 0.05)
    stringsAsFactors = FALSE
  )

  ivs <- select_instruments(marginal, cond, "A", "B", alpha1 = 0.01, alpha2 = 0.05, F_threshold = 10)

  # S1: marginal sig (1e-10 < 0.01), cond p=0.5 > 0.05, F=(1/0.1)^2=100 > 10 -> PASS
  # S2: marginal sig (1e-8 < 0.01), cond p=0.1 > 0.05, F=(0.5/0.1)^2=25 > 10 -> PASS
  # S3: marginal sig (1e-5 < 0.01), cond p=0.01 < 0.05 -> FAIL exclusivity
  expect_true("S1" %in% ivs)
  expect_true("S2" %in% ivs)
  expect_false("S3" %in% ivs)
})

test_that("instrument selection returns empty with no significant SNPs", {
  marginal <- data.frame(
    TRAIT = "A", QTL = "S1", A = "0.1", SE = "0.1", P_Value = "0.5",
    stringsAsFactors = FALSE
  )
  cond <- data.frame(trait = "B", locus = "S1", pval_cond = 0.5, stringsAsFactors = FALSE)

  expect_warning(
    ivs <- select_instruments(marginal, cond, "A", "B", alpha1 = 0.05),
    "无边际显著位点"
  )
  expect_equal(length(ivs), 0)
})

test_that("instrument selection warns when < 3 IVs", {
  marginal <- data.frame(
    TRAIT = c("A", "A"),
    QTL = c("S1", "S2"),
    A = c("1.0", "0.5"),
    SE = c("0.1", "0.1"),
    P_Value = c("1e-10", "1e-8"),
    stringsAsFactors = FALSE
  )
  cond <- data.frame(
    trait = c("B", "B"),
    locus = c("S1", "S2"),
    pval_cond = c(0.5, 0.5),
    stringsAsFactors = FALSE
  )

  expect_warning(
    ivs <- select_instruments(marginal, cond, "A", "B", alpha1 = 0.01, F_threshold = 1000),
    "有效 IV 数 = 2 < 3"
  )
})


# -----------------------------------------------------------------------------
# estimate_causal_gls
# -----------------------------------------------------------------------------

test_that("GLS estimate matches true causal effect", {
  # Simulate: theta_A = 1.0, theta_B = 0.5, true gamma = 0.5
  theta_A <- c(S1 = 1.0, S2 = 1.0, S3 = 1.0)
  theta_B <- c(S1 = 0.5, S2 = 0.5, S3 = 0.5)

  result <- estimate_causal_gls(theta_A, theta_B, ld_matrix = NULL)

  expect_equal(result$gamma, 0.5, tolerance = 1e-10)
  expect_equal(result$n_iv, 3)
})

test_that("GLS with no LD matrix uses identity", {
  theta_A <- c(S1 = 1.0, S2 = 0.5)
  theta_B <- c(S1 = 0.5, S2 = 0.25)

  result <- estimate_causal_gls(theta_A, theta_B)
  # gamma = sum(theta_A * theta_B) / sum(theta_A^2)
  expected <- sum(theta_A * theta_B) / sum(theta_A^2)
  expect_equal(result$gamma, expected, tolerance = 1e-10)
})

test_that("GLS requires at least 3 IVs", {
  theta_A <- c(S1 = 1.0, S2 = 0.5)
  theta_B <- c(S1 = 0.5, S2 = 0.25)
  expect_error(estimate_causal_gls(theta_A, theta_B), "k >= 3L")
})


# -----------------------------------------------------------------------------
# bidirectional_mr
# -----------------------------------------------------------------------------

test_that("bidirectional_mr detects A->B causality", {
  # Create mock data where A->B is true
  set.seed(42)
  n <- 1000
  X <- matrix(sample(c(-1, 1), n * 5, replace = TRUE), ncol = 5)
  colnames(X) <- paste0("S", 1:5)

  # True model: A = S1*1.0 + noise, B = A*0.3 + noise
  A <- X[, 1] + rnorm(n, 0, 1)
  B <- 0.3 * A + rnorm(n, 0, 1)
  Y <- cbind(A, B)
  colnames(Y) <- c("A", "B")

  marginal <- data.frame(
    TRAIT = rep(c("A", "B"), each = 5),
    QTL = rep(paste0("S", 1:5), 2),
    A = c(1.0, 0.1, 0.1, 0.1, 0.1,  # A effects
          0.3, 0.05, 0.05, 0.05, 0.05),  # B effects (mediated)
    SE = rep(0.1, 10),
    P_Value = c(1e-20, 0.1, 0.1, 0.1, 0.1,
                1e-10, 0.1, 0.1, 0.1, 0.1),
    stringsAsFactors = FALSE
  )

  cond <- data.frame(
    trait = rep(c("A", "B"), each = 5),
    locus = rep(paste0("S", 1:5), 2),
    pval_cond = c(0.5, 0.5, 0.5, 0.5, 0.5,  # A|B not significant
                  0.5, 0.5, 0.5, 0.5, 0.5),  # B|A not significant
    stringsAsFactors = FALSE
  )

  mr <- bidirectional_mr(marginal, cond, traits = c("A", "B"), 
                          Y = Y, X_all = X, alpha1 = 0.05, n_boot = 100)

  # A->B should be significant (gamma ~ 0.3)
  expect_true(!is.na(mr$AB$gamma))
  # B->A should not be significant (no reverse causality)
  expect_true(is.na(mr$BA$gamma) || mr$BA$pval > 0.05)
})

test_that("bidirectional_mr handles no valid IVs", {
  marginal <- data.frame(
    TRAIT = "A", QTL = "S1", A = "0.1", SE = "0.1", P_Value = "0.5",
    stringsAsFactors = FALSE
  )
  cond <- data.frame(trait = "B", locus = "S1", pval_cond = 0.5, stringsAsFactors = FALSE)

  mr <- bidirectional_mr(marginal, cond, traits = c("A", "B"), alpha1 = 0.05)
  expect_true(is.na(mr$AB$gamma))
  expect_true(is.na(mr$BA$gamma))
})


# -----------------------------------------------------------------------------
# test_mediation_consistency
# -----------------------------------------------------------------------------

test_that("mediation consistency detects effect attenuation", {
  # Marginal effect > conditional effect (mediation present)
  theta_marg <- c(S1 = 0.8, S2 = 0.7, S3 = 0.9)
  theta_cond <- c(S1 = 0.1, S2 = 0.2, S3 = 0.1)

  result <- test_mediation_consistency(theta_marg, theta_cond)

  expect_true(result$consistent)  # Effects systematically attenuated
  expect_true(result$delta_median > 0)
  expect_true(result$pval < 0.05)
})

test_that("mediation consistency handles no attenuation", {
  # No systematic difference
  theta_marg <- c(S1 = 0.5, S2 = 0.5, S3 = 0.5)
  theta_cond <- c(S1 = 0.5, S2 = 0.5, S3 = 0.5)

  result <- test_mediation_consistency(theta_marg, theta_cond)

  expect_false(result$consistent)
  expect_equal(result$delta_median, 0, tolerance = 1e-10)
})

test_that("mediation consistency handles shared IV names", {
  theta_marg <- c(S1 = 0.8, S2 = 0.7, S3 = 0.9, S4 = 0.6)
  theta_cond <- c(S1 = 0.1, S2 = 0.2, S3 = 0.1, S5 = 0.5)  # S4 missing, S5 extra

  result <- test_mediation_consistency(theta_marg, theta_cond)
  expect_equal(length(result$delta), 3)  # Only shared S1, S2, S3
})

test_that("mediation consistency warns with < 3 shared IVs", {
  theta_marg <- c(S1 = 0.8, S2 = 0.7)
  theta_cond <- c(S1 = 0.1, S2 = 0.2)

  expect_warning(
    result <- test_mediation_consistency(theta_marg, theta_cond),
    "共享工具变量 < 3"
  )
  expect_true(is.na(result$pval))
})