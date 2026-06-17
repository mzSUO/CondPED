# ==============================================================================
# Test Suite for simulate_data.R
# ==============================================================================

library(testthat)

# 加载被测代码
source("/Users/mingzhe/Desktop/CondPED/R/simulate_data.R")

# ------------------------------------------------------------------------------
# 1. sim_genotype 测试
# ------------------------------------------------------------------------------

test_that("sim_genotype produces correct dimensions and values", {
  X <- sim_genotype(n = 100, p = 50, maf = 0.3, population = "RIL", seed = 42)
  
  expect_equal(dim(X), c(100, 50))
  expect_true(all(X %in% c(-1, 1)))  # RIL 只有 -1/1
  expect_equal(colnames(X), paste0("SNP", 1:50))
  
  # F2 测试
  X2 <- sim_genotype(n = 100, p = 50, maf = 0.3, population = "F2", seed = 42)
  expect_true(all(X2 %in% c(-1, 0, 1)))
})

test_that("sim_genotype is reproducible with seed", {
  X1 <- sim_genotype(n = 100, p = 50, seed = 123)
  X2 <- sim_genotype(n = 100, p = 50, seed = 123)
  expect_equal(X1, X2)
})

# ------------------------------------------------------------------------------
# 2. sd_from_pve 测试
# ------------------------------------------------------------------------------

test_that("sd_from_pve matches analytical formula", {
  pve <- 0.01; maf <- 0.3
  var_x_ril <- 4 * maf * (1 - maf)
  var_x_f2  <- 2 * maf * (1 - maf)
  
  expect_equal(sd_from_pve(pve, maf, "RIL"), sqrt(pve / var_x_ril))
  expect_equal(sd_from_pve(pve, maf, "F2"),  sqrt(pve / var_x_f2))
})

# ------------------------------------------------------------------------------
# 3. compute_residual_variance 测试（核心）
# ------------------------------------------------------------------------------

test_that("M1 variance standardization yields Var(y) ≈ 1", {
  n <- 1000; p <- 100; maf <- 0.3
  X <- sim_genotype(n, p, maf, "RIL")
  B <- matrix(0, p, 2)
  B[1:2, 1] <- rnorm(2, 0, sd_from_pve(0.02, maf, "RIL"))
  B[3:4, 2] <- rnorm(2, 0, sd_from_pve(0.02, maf, "RIL"))
  
  var_x <- 4 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, matrix(0, 2, 2), var_x, rho = 0)
  out <- sim_phenotype_M1(X, B, Sigma_E)
  
  # 允许 5% 抽样误差
  expect_equal(var(out$Y[, 1]), 1, tolerance = 0.1)
  expect_equal(var(out$Y[, 2]), 1, tolerance = 0.1)
})

test_that("M2 (unidirectional) variance standardization yields Var(y) ≈ 1", {
  n <- 1000; p <- 100; maf <- 0.3; tau <- 0.3
  X <- sim_genotype(n, p, maf, "RIL")
  B <- matrix(0, p, 2)
  B[1:2, 1] <- rnorm(2, 0, sd_from_pve(0.02, maf, "RIL"))
  
  Tau <- matrix(0, 2, 2); Tau[2, 1] <- tau
  var_x <- 4 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, Tau, var_x, rho = 0)
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E)
  
  expect_equal(var(out$Y[, 1]), 1, tolerance = 0.1)
  expect_equal(var(out$Y[, 2]), 1, tolerance = 0.1)
})

test_that("M2 (bidirectional) variance standardization yields Var(y) ≈ 1", {
  n <- 1000; p <- 100; maf <- 0.3
  X <- sim_genotype(n, p, maf, "RIL")
  B <- matrix(0, p, 2)
  B[1:2, 1] <- rnorm(2, 0, sd_from_pve(0.01, maf, "RIL"))
  B[1:2, 2] <- rnorm(2, 0, sd_from_pve(0.01, maf, "RIL"))
  
  Tau <- matrix(0, 2, 2)
  Tau[2, 1] <- 0.2
  Tau[1, 2] <- 0.2
  
  var_x <- 4 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, Tau, var_x, rho = 0)
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E)
  
  expect_equal(var(out$Y[, 1]), 1, tolerance = 0.1)
  expect_equal(var(out$Y[, 2]), 1, tolerance = 0.1)
})

test_that("M2 with rho=0.5 residual correlation works", {
  n <- 1000; p <- 100; maf <- 0.3; tau <- 0.3
  X <- sim_genotype(n, p, maf, "RIL")
  B <- matrix(0, p, 2)
  B[1:2, 1] <- rnorm(2, 0, sd_from_pve(0.02, maf, "RIL"))
  
  Tau <- matrix(0, 2, 2); Tau[2, 1] <- tau
  var_x <- 4 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, Tau, var_x, rho = 0.5)
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E)
  
  # 方差接近 1
  expect_equal(var(out$Y[, 1]), 1, tolerance = 0.1)
  expect_equal(var(out$Y[, 2]), 1, tolerance = 0.1)
  # 残差相关接近 0.5
  cor_eps <- cor(out$eps[, 1], out$eps[, 2])
  expect_equal(cor_eps, 0.5, tolerance = 0.1)
})

# ------------------------------------------------------------------------------
# 4. 三数据集生成函数测试
# ------------------------------------------------------------------------------

test_that("Dataset I structure is correct", {
  dat <- generate_dataset_I(n = 500, p = 100, seed = 1)
  
  expect_equal(dim(dat$Y), c(500, 2))
  expect_equal(nrow(dat$truth), 100)
  expect_equal(sum(dat$truth$class == "class1"), 2)
  expect_equal(sum(dat$truth$class == "class2"), 2)
  expect_equal(sum(dat$truth$class == "null"), 96)
  expect_equal(sum(dat$truth$is_IV), 0)
  expect_equal(dat$dataset, "I")
})

test_that("Dataset II structure is correct", {
  dat <- generate_dataset_II(n = 500, p = 100, seed = 1)
  
  expect_equal(sum(dat$truth$class == "class1"), 2)
  expect_equal(sum(dat$truth$class == "class3"), 2)
  expect_equal(sum(dat$truth$class == "class4"), 2)
  expect_equal(sum(dat$truth$class == "IV_A"), 10)
  expect_equal(sum(dat$truth$is_IV), 10)
  expect_equal(dat$dataset, "II")
  expect_equal(dat$Tau[2, 1], 0.3)
})

test_that("Dataset III structure is correct", {
  dat <- generate_dataset_III(n = 500, p = 100, seed = 1)
  
  expect_equal(sum(dat$truth$class == "class5"), 2)
  expect_equal(sum(dat$truth$class == "IV_A"), 10)
  expect_equal(sum(dat$truth$class == "IV_B"), 10)
  expect_equal(sum(dat$truth$is_IV), 20)
  expect_equal(dat$dataset, "III")
  expect_equal(dat$Tau[2, 1], 0.2)
  expect_equal(dat$Tau[1, 2], 0.2)
})

test_that("Class 5 loci have both beta_A and beta_B non-zero", {
  dat <- generate_dataset_III(n = 500, p = 100, seed = 1)
  class5_idx <- which(dat$truth$class == "class5")
  
  expect_true(all(dat$truth$beta_A[class5_idx] != 0))
  expect_true(all(dat$truth$beta_B[class5_idx] != 0))
})

# ------------------------------------------------------------------------------
# 5. 边界条件测试
# ------------------------------------------------------------------------------

test_that("small n and p work", {
  dat <- generate_dataset_I(n = 50, p = 20, seed = 1)
  expect_equal(dim(dat$Y), c(50, 2))
  expect_equal(nrow(dat$truth), 20)
})

test_that("large p works", {
  dat <- generate_dataset_I(n = 200, p = 10000, seed = 1)
  expect_equal(dim(dat$Y), c(200, 2))
  expect_equal(nrow(dat$truth), 10000)
})

cat("\n=== All tests passed ===\n")