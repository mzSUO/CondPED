# 保存为: tests/testthat/test-simulate.R

# ========== 测试 1: M1 方差标准化 ==========
test_that("M1 produces variance near 1", {
  set.seed(42)
  X <- sim_genotype(n = 500, p = 100)
  B <- matrix(0, 100, 2)
  B[1:10, 1] <- rnorm(10, mean = 0, sd = 0.3)
  
  out <- sim_phenotype_M1(X, B, Sigma_E = diag(2), maf = 0.3)
  
  vars <- apply(out$Y, 2, var)
  expect_lt(abs(vars[1] - 1), 0.15)
  expect_lt(abs(vars[2] - 1), 0.15)
})

# ========== 测试 2: M2 保留真实 tau（Class 4 完全中介） ==========
test_that("M2 preserves true tau for complete mediation", {
  set.seed(42)
  dat <- generate_class4(n = 5000, pve_A = 0.02, tau = 0.3)
  
  fit <- lm(dat$Y[, 2] ~ dat$Y[, 1])
  tau_est <- coef(fit)[2]
  
  expect_lt(abs(tau_est - 0.3), 0.03)
  expect_lt(abs(dat$params$tau_observed - 0.3), 0.03)
})

# ========== 测试 3: Class 5 部分中介保留直接效应（修正版） ==========
test_that("Class 5 preserves direct effect on B", {
  set.seed(42)
  dat <- generate_class5(n = 5000, pve_A = 0.02, pve_B = 0.003, tau = 0.3)
  
  # 修正：控制 Y_A 后，检验 B 的纯直接效应（与间接效应分离）
  resid_B <- residuals(lm(dat$Y[, 2] ~ dat$Y[, 1]))
  pvals <- apply(dat$X[, 1:30], 2, function(x) summary(lm(resid_B ~ x))$coef[2, 4])
  
  # 弱直接效应（PVE=0.3%）在 n=5000 时应有 10-60% 检出率
  expect_gt(mean(pvals < 0.05), 0.10)
  expect_lt(mean(pvals < 0.05), 0.80)
})

# ========== 测试 4: Class 6 双向因果（修正版：检验方向显著性） ==========
test_that("M2 handles bidirectional causality", {
  set.seed(42)
  dat <- generate_class6(n = 5000, pve = 0.01, tau_AB = 0.2, tau_BA = 0.2)
  
  fit_ab <- lm(dat$Y[, 2] ~ dat$Y[, 1])
  fit_ba <- lm(dat$Y[, 1] ~ dat$Y[, 2])
  
  tau_obs_ab <- coef(fit_ab)[2]
  tau_obs_ba <- coef(fit_ba)[2]
  
  # 修正：双向因果下，两个回归都应显著且为正
  # 简化式系数不等于结构参数 tau（反馈回路放大），故不检验数值
  expect_lt(summary(fit_ab)$coef[2, 4], 0.05)  # B~A 显著
  expect_lt(summary(fit_ba)$coef[2, 4], 0.05)  # A~B 显著
  expect_gt(tau_obs_ab, 0.25)  # 反馈放大后应 > 0.25
  expect_gt(tau_obs_ba, 0.25)
})

# ========== 测试 5: Class 2 压力测试协方差结构（修正版） ==========
test_that("Covariance stress test produces expected phenotypic correlation", {
  set.seed(42)
  dat <- generate_covariance_stress(n = 1000, rho = 0.5, pve = 0.01)
  
  # 修正：Sigma_E 的 rho=0.5 是残差相关。
  # 由于性状 A 有遗传方差（h^2 ≈ 0.30），表型相关 ≈ rho * sqrt(1-h^2_A) ≈ 0.42
  rho_obs <- cor(dat$Y[, 1], dat$Y[, 2])
  expect_lt(abs(rho_obs - 0.42), 0.08)
})

# ========== 测试 6: 统一入口函数 ==========
test_that("generate_condped returns correct class labels", {
  set.seed(42)
  
  dat1 <- generate_condped("class1", n = 100, p = 50)
  expect_equal(dat1$class_label, "class1")
  expect_true(!is.null(dat1$X))
  
  dat4 <- generate_condped("class4", n = 100, p = 50, pve_A = 0.02)
  expect_equal(dat4$class_label, "class4")
  expect_true(!is.null(dat4$X))
  
  dat_stress <- generate_condped("stress_test", n = 100, p = 50, rho = 0.7)
  expect_equal(dat_stress$class_label, "stress_test")
  expect_true(!is.null(dat_stress$X))
})

# ========== 测试 7: M1 带残差相关时方差仍受控（新增） ==========
test_that("M1 with rho=0.5 still standardizes variance", {
  set.seed(42)
  X <- sim_genotype(500, 100)
  B <- matrix(0, 100, 2)
  B[1:10, 1] <- rnorm(10, 0, 0.3)
  Sigma_E <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  out <- sim_phenotype_M1(X, B, Sigma_E, maf = 0.3)
  
  expect_lt(abs(var(out$Y[, 1]) - 1), 0.15)
  expect_lt(abs(var(out$Y[, 2]) - 1), 0.15)
})

# ========== 测试 8: Class 4 条件效应接近 0（CondPED 核心假设，新增） ==========
test_that("Class 4 has near-zero conditional B effect", {
  set.seed(42)
  dat <- generate_class4(n = 2000, pve_A = 0.02, tau = 0.3)
  
  # 条件效应：控制 Y_A 后，X 对 Y_B 的回归系数应接近 0（真实 beta_B = 0）
  cond_effects <- apply(dat$X[, 1:30], 2, function(x) {
    coef(lm(dat$Y[, 2] ~ dat$Y[, 1] + x))[3]
  })
  
  expect_lt(mean(abs(cond_effects)), 0.05)
})