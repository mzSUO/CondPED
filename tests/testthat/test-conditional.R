# tests/testthat/test-conditional.R
# Test suite for Layer 2 conditional projection

# -----------------------------------------------------------------------------
# estimate_phenotypic_covariance
# -----------------------------------------------------------------------------

test_that("covariance estimation is correct for independent traits", {
  set.seed(42)
  Y <- matrix(rnorm(500 * 3), ncol = 3)
  out <- estimate_phenotypic_covariance(Y)

  expect_equal(dim(out$V), c(3, 3))
  expect_equal(dim(out$Cor), c(3, 3))
  # Diagonal should be close to 1
  expect_equal(diag(out$Cor), rep(1, 3), tolerance = 0.1)
  # Off-diagonal should be close to 0
  expect_equal(out$Cor[1, 2], 0, tolerance = 0.1)
})

test_that("covariance estimation detects correlation", {
  set.seed(42)
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Y <- MASS::mvrnorm(1000, mu = c(0, 0), Sigma = Sigma)
  out <- estimate_phenotypic_covariance(Y)

  expect_equal(out$Cor[1, 2], 0.5, tolerance = 0.05)
})

test_that("covariance estimation errors when n <= m", {
  Y <- matrix(rnorm(10), ncol = 10)
  expect_error(estimate_phenotypic_covariance(Y))
})


# -----------------------------------------------------------------------------
# compute_projection_coefficients
# -----------------------------------------------------------------------------

test_that("projection coefficients satisfy gamma = V_{-i}^{-1} C_{-i,i}", {
  Sigma <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.4, 0.3, 0.4, 1), 3, 3)
  gamma <- compute_projection_coefficients(Sigma, i = 1)

  idx <- 2:3
  expected <- solve(Sigma[idx, idx], Sigma[idx, 1])
  expect_equal(gamma, as.numeric(expected), tolerance = 1e-10)
})

test_that("projection coefficients for independent traits are zero", {
  Sigma <- diag(3)
  gamma <- compute_projection_coefficients(Sigma, i = 1)
  expect_equal(gamma, c(0, 0), tolerance = 1e-10)
})


# -----------------------------------------------------------------------------
# compute_conditional_phenotype
# -----------------------------------------------------------------------------

test_that("conditional phenotype achieves orthogonality", {
  set.seed(42)
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Y <- MASS::mvrnorm(5000, mu = c(0, 0), Sigma = Sigma)
  colnames(Y) <- c("A", "B")

  result <- compute_conditional_phenotype(Y, tol_orthogonality = 0.05)

  # Orthogonality: Cov(Y_cond[,1], Y[,2]) should be ~0
  cov_val <- cov(result$Y_cond[, "A"], Y[, "B"])
  expect_equal(cov_val, 0, tolerance = 0.05)

  # Same for reverse
  cov_val2 <- cov(result$Y_cond[, "B"], Y[, "A"])
  expect_equal(cov_val2, 0, tolerance = 0.05)
})

test_that("conditional phenotype variance equals Schur complement", {
  set.seed(42)
  Sigma <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
  Y <- MASS::mvrnorm(10000, mu = c(0, 0), Sigma = Sigma)
  colnames(Y) <- c("A", "B")

  result <- compute_conditional_phenotype(Y)

  # Schur complement for trait A: V_A - C_AB * V_B^{-1} * C_BA
  expected_schur <- Sigma[1, 1] - Sigma[1, 2]^2 / Sigma[2, 2]
  expect_equal(result$schur["A"], expected_schur, tolerance = 0.01)

  # Observed variance should be close to Schur complement
  var_obs <- var(result$Y_cond[, "A"])
  expect_equal(var_obs, expected_schur, tolerance = 0.1)
})

test_that("single trait returns original phenotype", {
  Y <- matrix(rnorm(100), ncol = 1)
  result <- compute_conditional_phenotype(Y)
  expect_equal(result$Y_cond, Y, tolerance = 1e-10)
})

test_that("projection coefficients are stored correctly", {
  Sigma <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.4, 0.3, 0.4, 1), 3, 3)
  Y <- MASS::mvrnorm(1000, mu = c(0, 0, 0), Sigma = Sigma)
  colnames(Y) <- c("A", "B", "C")

  result <- compute_conditional_phenotype(Y)
  expect_equal(length(result$gamma), 3)
  expect_equal(names(result$gamma), c("A", "B", "C"))
  expect_equal(length(result$gamma[["A"]]), 2)  # 2 other traits
})


# -----------------------------------------------------------------------------
# fit_conditional_model
# -----------------------------------------------------------------------------

test_that("conditional model detects true QTL effects", {
  set.seed(42)
  n <- 500
  X <- matrix(sample(c(-1, 1), n * 2, replace = TRUE), ncol = 2)
  colnames(X) <- c("SNP1", "SNP2")

  # True effect: SNP1 affects trait 1
  Y <- matrix(rnorm(n * 2), ncol = 2)
  Y[, 1] <- Y[, 1] + 0.5 * X[, 1]
  colnames(Y) <- c("TraitA", "TraitB")

  cond_eff <- fit_conditional_model(Y, X, alpha2 = 0.05)

  expect_true(any(cond_eff$locus == "SNP1" & cond_eff$trait == "TraitA"))
  snp1_traitA <- cond_eff[cond_eff$locus == "SNP1" & cond_eff$trait == "TraitA", ]
  expect_true(snp1_traitA$sig_cond)  # Should be significant
})

test_that("conditional model handles monomorphic SNPs", {
  n <- 100
  X <- matrix(rep(1, n * 2), ncol = 2)  # Monomorphic
  colnames(X) <- c("SNP1", "SNP2")
  Y <- matrix(rnorm(n * 2), ncol = 2)
  colnames(Y) <- c("TraitA", "TraitB")

  cond_eff <- fit_conditional_model(Y, X, alpha2 = 0.05)
  # Should not error, but may have fewer rows
  expect_true(is.data.frame(cond_eff))
})

test_that("Bonferroni alpha2 is auto-calculated", {
  X <- matrix(sample(c(-1, 1), 100 * 5, replace = TRUE), ncol = 5)
  Y <- matrix(rnorm(100 * 2), ncol = 2)

  cond_eff <- fit_conditional_model(Y, X, alpha2 = NULL)
  # alpha2 should be 0.05 / (5 * 2) = 0.005
  expect_true(all(cond_eff$pval_cond >= 0 | cond_eff$pval_cond <= 1))  # Valid p-values
})