test_that("Data simulation produces correct dimensions", {
  
  data <- simulate_multitrait_data(
    n_ind = 100,
    n_trait = 3,
    n_snp = 50
  )
  
  # Check dimensions
  expect_equal(nrow(data$Y), 100)
  expect_equal(ncol(data$Y), 3)
  expect_equal(nrow(data$G), 100)
  expect_equal(ncol(data$G), 50)
  expect_equal(dim(data$beta_true), c(50, 3))
  
  # Check covariance matrix is symmetric
  expect_true(isSymmetric(data$V_true))
})

test_that("Effect decomposition is correct", {
  
  V <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  beta <- matrix(c(0.3, 0, 0.15, 0), 2, 2)
  
  decomp <- CondPED:::decompose_true_effects(beta, V)
  
  # Check reconstruction
  reconstructed <- decomp$beta_ind + decomp$beta_shared
  expect_equal(beta, reconstructed, tolerance = 1e-10)
})