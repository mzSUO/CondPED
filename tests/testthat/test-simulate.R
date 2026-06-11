# tests/testthat/test-simulate_data.R
# Test suite for CondPED simulation system

# -----------------------------------------------------------------------------
# sim_genotype
# -----------------------------------------------------------------------------

test_that("sim_genotype produces RIL with only -1 and 1", {
  X <- sim_genotype(n = 100, p = 50, maf = 0.3, population = "RIL", seed = 42)
  expect_equal(dim(X), c(100, 50))
  expect_true(all(X %in% c(-1, 1)))
  expect_equal(colnames(X), paste0("SNP", 1:50))
  expect_equal(rownames(X), paste0("Ind", 1:100))
})

test_that("sim_genotype produces F2 with -1, 0, 1", {
  X <- sim_genotype(n = 100, p = 50, maf = 0.3, population = "F2", seed = 42)
  expect_true(all(X %in% c(-1, 0, 1)))
})

test_that("sim_genotype respects maf parameter", {
  X <- sim_genotype(n = 10000, p = 1000, maf = 0.5, population = "RIL", seed = 42)
  # For RIL, P(X=1) = maf, so mean should be close to 2*maf - 1 = 0
  expect_equal(mean(X), 0, tolerance = 0.02)
})

test_that("sim_genotype rejects invalid inputs", {
  expect_error(sim_genotype(n = 0, p = 10))
  expect_error(sim_genotype(n = 10, p = 0))
  expect_error(sim_genotype(n = 10, p = 10, maf = 1.5))
  expect_error(sim_genotype(n = 10, p = 10, maf = -0.1))
})


# -----------------------------------------------------------------------------
# sd_from_pve
# -----------------------------------------------------------------------------

test_that("sd_from_pve returns correct formula for RIL", {
  # RIL: Var(X) = 4*maf*(1-maf), so sd = sqrt(pve / (4*maf*(1-maf)))
  expect_equal(sd_from_pve(0.01, 0.3, "RIL"), sqrt(0.01 / (4 * 0.3 * 0.7)))
})

test_that("sd_from_pve returns correct formula for F2", {
  # F2: Var(X) = 2*maf*(1-maf)
  expect_equal(sd_from_pve(0.01, 0.3, "F2"), sqrt(0.01 / (2 * 0.3 * 0.7)))
})


# -----------------------------------------------------------------------------
# sim_phenotype_M1
# -----------------------------------------------------------------------------

test_that("M1 generates correct additive effects only", {
  set.seed(42)
  X <- sim_genotype(n = 500, p = 100, maf = 0.3, population = "RIL")
  B <- matrix(0, 100, 2)
  B[1:10, 1] <- 1.0
  out <- sim_phenotype_M1(X, B, Sigma_E = diag(2), maf = 0.3, population = "RIL")

  expect_equal(dim(out$Y), c(500, 2))
  expect_equal(dim(out$G), c(500, 2))
  # G should equal X %*% B
  expect_equal(out$G, X %*% B, tolerance = 1e-10)
  # Variance should be close to 1
  expect_equal(var(out$Y[, 1]), 1, tolerance = 0.2)
})

test_that("M1 includes epistasis when epi_pairs provided", {
  set.seed(42)
  X <- sim_genotype(n = 500, p = 100, maf = 0.3, population = "RIL")
  B <- matrix(0, 100, 2)
  out <- sim_phenotype_M1(X, B, Sigma_E = diag(2),
                          epi_pairs = list(c(1, 2)), epi_effects = c(0.5),
                          maf = 0.3, population = "RIL")

  # G_epi should be non-zero
  expect_true(any(out$G_epi != 0))
  # G should equal G_add + G_epi
  expect_equal(out$G, out$G_add + out$G_epi, tolerance = 1e-10)
})

test_that("M1 errors when epi_effects missing but epi_pairs present", {
  X <- sim_genotype(n = 100, p = 50, population = "RIL")
  B <- matrix(0, 50, 2)
  expect_error(sim_phenotype_M1(X, B, diag(2), epi_pairs = list(c(1, 2))))
})

test_that("M1 errors when epi_pairs and epi_effects length mismatch", {
  X <- sim_genotype(n = 100, p = 50, population = "RIL")
  B <- matrix(0, 50, 2)
  expect_error(sim_phenotype_M1(X, B, diag(2), 
                                 epi_pairs = list(c(1, 2), c(3, 4)),
                                 epi_effects = c(0.5)))
})


# -----------------------------------------------------------------------------
# sim_phenotype_M2
# -----------------------------------------------------------------------------

test_that("M2 generates correct causal chain", {
  set.seed(42)
  X <- sim_genotype(n = 2000, p = 100, maf = 0.3, population = "RIL")
  B <- matrix(0, 100, 2)
  B[1:10, 1] <- 1.0
  Tau <- matrix(0, 2, 2); Tau[2, 1] <- 0.3

  out <- sim_phenotype_M2(X, B, Tau, Sigma_E = diag(2), maf = 0.3, population = "RIL")

  # Verify causal effect: Y[,2] ~ Y[,1] should have slope ~ 0.3
  fit <- lm(out$Y[, 2] ~ out$Y[, 1])
  expect_equal(coef(fit)[2], 0.3, tolerance = 0.05)

  # Y_base should equal G + E (before causal transformation)
  expect_equal(out$Y_base, out$G + out$E, tolerance = 1e-10)
})

test_that("M2 includes epistasis and causal path", {
  set.seed(42)
  X <- sim_genotype(n = 500, p = 100, maf = 0.3, population = "RIL")
  B <- matrix(0, 100, 2)
  Tau <- matrix(0, 2, 2); Tau[2, 1] <- 0.3

  out <- sim_phenotype_M2(X, B, Tau, Sigma_E = diag(2),
                           epi_pairs = list(c(1, 2)), epi_effects = c(0.5),
                           maf = 0.3, population = "RIL")

  expect_true(any(out$G_epi != 0))
  expect_equal(out$G, out$G_add + out$G_epi, tolerance = 1e-10)
})


# -----------------------------------------------------------------------------
# generate_class*
# -----------------------------------------------------------------------------

test_that("generate_class1 produces single-trait QTL", {
  sim <- generate_class1(n = 200, p = 100, pve = 0.01, seed = 42)
  expect_equal(sim$class_label, "class1")
  # Only first 30 SNPs should affect trait 1
  expect_true(all(sim$B[31:100, 1] == 0))
  expect_true(all(sim$B[, 2] == 0))
  # Phenotype should be generated
  expect_equal(dim(sim$Y), c(200, 2))
})

test_that("generate_class3 produces independent multi-trait QTL", {
  sim <- generate_class3(n = 200, p = 100, pve = 0.01, seed = 42)
  expect_equal(sim$class_label, "class3")
  # Same SNPs affect both traits independently
  expect_true(all(sim$B[1:30, 1] != 0))
  expect_true(all(sim$B[1:30, 2] != 0))
})

test_that("generate_class4 produces causal chain A->B", {
  sim <- generate_class4(n = 2000, p = 100, pve_A = 0.02, tau = 0.3, seed = 42)
  expect_equal(sim$class_label, "class4")
  # Trait A should have QTL effects
  expect_true(any(sim$B[, 1] != 0))
  # Trait B should have no direct QTL effects (only mediated)
  expect_true(all(sim$B[, 2] == 0))
  # Verify causal structure
  fit <- lm(sim$Y[, 2] ~ sim$Y[, 1])
  expect_equal(coef(fit)[2], 0.3, tolerance = 0.05)
})

test_that("generate_class5 produces partial mediation", {
  sim <- generate_class5(n = 2000, p = 100, pve_A = 0.02, pve_B = 0.003, 
                          tau = 0.3, seed = 42)
  expect_equal(sim$class_label, "class5")
  # Both traits have direct QTL effects
  expect_true(any(sim$B[, 1] != 0))
  expect_true(any(sim$B[, 2] != 0))
})

test_that("generate_class6 produces bidirectional causality", {
  sim <- generate_class6(n = 2000, p = 100, tau_AB = 0.2, tau_BA = 0.2, seed = 42)
  expect_equal(sim$class_label, "class6")
  expect_true(all(sim$Tau == matrix(c(0, 0.2, 0.2, 0), 2, 2)))
})

test_that("generate_covariance_stress produces residual correlation", {
  sim <- generate_covariance_stress(n = 200, p = 100, rho = 0.5, seed = 42)
  expect_equal(sim$class_label, "stress_test")
  # Check residual correlation structure
  cor_mat <- cor(sim$E)
  expect_equal(cor_mat[1, 2], 0.5, tolerance = 0.1)
})


# -----------------------------------------------------------------------------
# generate_condped (unified entry)
# -----------------------------------------------------------------------------

test_that("generate_condped dispatches correctly", {
  sim1 <- generate_condped("class1", n = 100, p = 50, seed = 42)
  expect_equal(sim1$class_label, "class1")

  sim4 <- generate_condped("class4", n = 100, p = 50, pve_A = 0.02, tau = 0.3, seed = 42)
  expect_equal(sim4$class_label, "class4")
})

test_that("generate_condped rejects invalid scenario", {
  expect_error(generate_condped("invalid", n = 100, p = 50))
})