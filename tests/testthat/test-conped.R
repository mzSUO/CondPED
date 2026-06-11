# tests/testthat/test-condped.R
# Test suite for main controller function

# -----------------------------------------------------------------------------
# .remove_covariates (internal helper)
# -----------------------------------------------------------------------------

test_that("covariate removal produces residuals", {
  n <- 100
  Y <- matrix(rnorm(n * 2), ncol = 2)
  cov <- matrix(rnorm(n * 2), ncol = 2)

  resid <- CondPED:::.remove_covariates(Y, cov)

  expect_equal(dim(resid), dim(Y))
  # Residuals should be uncorrelated with covariates
  for (j in 1:2) {
    for (k in 1:2) {
      expect_equal(cor(resid[, j], cov[, k]), 0, tolerance = 0.1)
    }
  }
})


# -----------------------------------------------------------------------------
# condped (integration test with simulation data)
# -----------------------------------------------------------------------------

test_that("condped runs end-to-end on simulated data", {
  skip_if_not(file.exists("/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork"),
              "QTXNetwork not available")

  sim <- generate_class4(n = 100, p = 50, pve_A = 0.02, tau = 0.3, seed = 42)

  result <- condped(
    geno = sim$X,
    pheno = as.data.frame(sim$Y),
    qtxnetwork_path = "/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork",
    traits = c("TraitA", "TraitB"),
    output_dir = tempfile(),
    verbose = FALSE
  )

  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("layer1", "layer2", "cond_effects", "layer3", "classification", "params") %in% names(result)))

  # Check classification
  expect_true(is.data.frame(result$classification))
  expect_true("class" %in% colnames(result$classification))

  # Check output files exist
  expect_true(file.exists(file.path(result$params$output_dir, "condped_result.rds")))
  expect_true(file.exists(file.path(result$params$output_dir, "classification.csv")))
})

test_that("condped handles class1-only data (no multi-trait SNPs)", {
  skip_if_not(file.exists("/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork"),
              "QTXNetwork not available")

  sim <- generate_class1(n = 100, p = 50, pve = 0.01, seed = 42)

  result <- condped(
    geno = sim$X,
    pheno = as.data.frame(sim$Y),
    qtxnetwork_path = "/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork",
    output_dir = tempfile(),
    verbose = FALSE
  )

  # Should return early with only layer1 and classification
  expect_true(is.list(result))
  expect_true("layer1" %in% names(result))
  expect_true("classification" %in% names(result))
})

test_that("condped handles pheno with ID column", {
  skip_if_not(file.exists("/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork"),
              "QTXNetwork not available")

  sim <- generate_class1(n = 50, p = 30, seed = 42)
  pheno_with_id <- data.frame(
    id = paste0("Ind", 1:50),
    sim$Y,
    stringsAsFactors = FALSE
  )

  # Should not error
  expect_error(condped(
    geno = sim$X,
    pheno = pheno_with_id,
    qtxnetwork_path = "/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork",
    output_dir = tempfile(),
    verbose = FALSE
  ), NA)
})

test_that("condped rejects invalid geno encoding", {
  sim <- generate_class1(n = 50, p = 30, seed = 42)
  sim$X[1, 1] <- 2  # Invalid encoding

  expect_error(condped(
    geno = sim$X,
    pheno = as.data.frame(sim$Y),
    qtxnetwork_path = "/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork",
    output_dir = tempfile()
  ))
})

test_that("condped requires qtxnetwork_path for layer1", {
  sim <- generate_class1(n = 50, p = 30, seed = 42)

  expect_error(condped(
    geno = sim$X,
    pheno = as.data.frame(sim$Y),
    output_dir = tempfile()
  ), "qtxnetwork_path is required")
})

test_that("condped output is reproducible with same seed", {
  skip_if_not(file.exists("/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork"),
              "QTXNetwork not available")

  sim <- generate_class4(n = 100, p = 50, pve_A = 0.02, tau = 0.3, seed = 42)

  res1 <- condped(geno = sim$X, pheno = as.data.frame(sim$Y),
                  qtxnetwork_path = "/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork",
                  output_dir = tempfile(), verbose = FALSE)

  res2 <- condped(geno = sim$X, pheno = as.data.frame(sim$Y),
                  qtxnetwork_path = "/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork",
                  output_dir = tempfile(), verbose = FALSE)

  # Same input should produce same classification
  expect_equal(res1$classification$class, res2$classification$class)
})


# -----------------------------------------------------------------------------
# condped with covariates
# -----------------------------------------------------------------------------

test_that("condped handles covariates correctly", {
  skip_if_not(file.exists("/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork"),
              "QTXNetwork not available")

  sim <- generate_class4(n = 100, p = 50, pve_A = 0.02, tau = 0.3, seed = 42)
  covariates <- data.frame(
    sex = sample(c(0, 1), 100, replace = TRUE),
    age = rnorm(100)
  )

  result <- condped(
    geno = sim$X,
    pheno = as.data.frame(sim$Y),
    qtxnetwork_path = "/data2/smz/QTXNetwork_4.0_MT/build/QTXNetwork",
    covariates = covariates,
    output_dir = tempfile(),
    verbose = FALSE
  )

  expect_true(is.list(result))
  expect_true("classification" %in% names(result))
})