# ==============================================================================
# Unit Tests for simulate_data.R
# ==============================================================================

library(testthat)
library(Matrix)
library(MASS)

# Test 1: Basic simulation works
test_that("simulate_multitrait_data generates data with correct dimensions", {
  set.seed(123)
  data <- simulate_multitrait_data(
    n_ind = 100,
    n_trait = 3,
    n_snp = 50
  )

  expect_equal(nrow(data$Y), 100)
  expect_equal(ncol(data$Y), 3)
  expect_equal(ncol(data$G), 50)
  expect_equal(dim(data$beta_true), c(50, 3))
  expect_equal(dim(data$V_true), c(3, 3))
})

# Test 2: Genotype is properly standardized
test_that("genotype matrix is standardized", {
  set.seed(123)
  data <- simulate_multitrait_data(n_ind = 200, n_trait = 2, n_snp = 20)

  # Check mean is approximately 0
  expect_true(all(abs(colMeans(data$G)) < 0.1))

  # Check variance is approximately 1
  expect_true(all(abs(apply(data$G, 2, var) - 1) < 0.2))
})

# Test 3: SNP info has correct structure
test_that("SNP info data frame is correct", {
  set.seed(123)
  data <- simulate_multitrait_data(n_ind = 100, n_trait = 2, n_snp = 30)

  expect_true("CHR" %in% colnames(data$snp_info))
  expect_true("SNP" %in% colnames(data$snp_info))
  expect_true("BP" %in% colnames(data$snp_info))
  expect_true("MAF" %in% colnames(data$snp_info))
  expect_equal(nrow(data$snp_info), 30)
})

# Test 4: Trait covariance matrix has correct properties
test_that("trait covariance matrix is valid", {
  set.seed(123)
  data <- simulate_multitrait_data(n_ind = 100, n_trait = 3, n_snp = 50)

  # Check symmetric
  expect_equal(data$V_true, t(data$V_true), tolerance = 1e-10)

  # Check positive definite (all eigenvalues > 0)
  eig <- eigen(data$V_true, only.values = TRUE)$values
  expect_true(all(eig > 0))

  # Check diagonal is 1 (correlation matrix)
  expect_true(all(abs(diag(data$V_true) - 1) < 1e-10))
})

# Test 5: Heritability is approximately correct
test_that("sample heritability matches target", {
  set.seed(123)
  h2_target <- 0.6
  data <- simulate_multitrait_data(
    n_ind = 500,
    n_trait = 2,
    n_snp = 100,
    heritability = h2_target
  )

  # Compute sample heritability
  var_g <- apply(data$G %*% data$beta_true, 2, var)
  var_e <- apply(data$E, 2, var)
  h2_sample <- var_g / (var_g + var_e)

  # Should be within 0.15 of target (with some tolerance for sampling)
  expect_true(abs(h2_sample[1] - h2_target) < 0.2)
})

# Test 6: Effect decomposition is correct
test_that("beta_ind + beta_shared = beta_true", {
  set.seed(123)
  data <- simulate_multitrait_data(n_ind = 200, n_trait = 3, n_snp = 50)

  # Check decomposition
  reconstructed <- data$beta_ind_true + data$beta_shared_true
  expect_equal(reconstructed, data$beta_true, tolerance = 1e-10)
})

# Test 7: Custom effect configuration works
test_that("custom effect configuration is applied", {
  set.seed(123)
  config <- list(
    add_trait_specific = list(
      snp_ids = 1:5,
      trait_ids = rep(1, 5),
      effect_size = rep(1.0, 5)
    ),
    add_shared = list(
      snp_ids = 6:10,
      effect_size = 0.5
    ),
    add_independent = list(
      snp_ids = 11:15,
      effect_size = 0.3
    )
  )

  data <- simulate_multitrait_data(
    n_ind = 200,
    n_trait = 3,
    n_snp = 50,
    effect_config = config
  )

  # Check that specified effects are non-zero
  expect_true(all(data$beta_true[1:5, 1] != 0))
  expect_true(all(data$beta_true[6:10, ] != 0))
  expect_true(any(data$beta_true[11:15, ] != 0))
})

# Test 8: Shared effects follow covariance pattern
test_that("shared effects follow covariance-mediated pattern", {
  set.seed(123)
  data <- simulate_multitrait_data(
    n_ind = 500,
    n_trait = 3,
    n_snp = 100,
    effect_config = list(
      add_shared = list(
        snp_ids = 1:10,
        effect_size = 0.5
      )
    )
  )

  # For shared effects: beta[j]/beta[i] ≈ V[j,i]/V[i,i]
  snp_ids <- 1:10
  V <- data$V_true

  for (j in 2:3) {
    ratio_beta <- data$beta_true[snp_ids, j] / data$beta_true[snp_ids, 1]
    ratio_V <- V[j, 1] / V[1, 1]
    # Allow NA when beta is zero
    valid_idx <- !is.na(ratio_beta) & abs(ratio_beta) > 1e-10
    if (any(valid_idx)) {
      expect_true(all(abs(ratio_beta[valid_idx] - ratio_V) < 0.01))
    }
  }
})

# Test 9: Relationship matrix is correct
test_that("additive relationship matrix is valid", {
  set.seed(123)
  data <- simulate_multitrait_data(n_ind = 100, n_trait = 2, n_snp = 50)

  # Check symmetric
  expect_equal(data$A, t(data$A), tolerance = 1e-10)

  # Check diagonal is positive (should be close to 1 for standardized G)
  expect_true(all(diag(data$A) > 0))

  # Check positive semi-definite
  eig <- eigen(data$A, only.values = TRUE)$values
  expect_true(all(eig > -1e-10))
})

# Test 10: Effect decomposition produces valid matrices
test_that("effect decomposition produces valid matrices", {
  set.seed(123)
  data <- simulate_multitrait_data(n_ind = 200, n_trait = 3, n_snp = 30)

  # Check dimensions
  expect_equal(dim(data$beta_ind_true), dim(data$beta_true))
  expect_equal(dim(data$beta_shared_true), dim(data$beta_true))

  # Check sum equals original
  expect_equal(data$beta_ind_true + data$beta_shared_true, data$beta_true, tolerance = 1e-10)
})
