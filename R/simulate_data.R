# ==============================================================================
# CondPED Simulation System (三数据集 m=2 方案)
# ==============================================================================
# 修正：M2 方差标准化公式 Sigma_Y = A^T %*% Sigma_base %*% A
# ==============================================================================

library(MASS)

# ------------------------------------------------------------------------------
# 1. 基础函数
# ------------------------------------------------------------------------------

#' 模拟基因型矩阵
#' @param n 样本量
#' @param p SNP 数
#' @param maf 次等位基因频率
#' @param population "RIL"(-1/1) 或 "F2"(-1/0/1)
sim_genotype <- function(n, p, maf = 0.3, population = c("RIL", "F2"), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  population <- match.arg(population)
  
  if (population == "RIL") {
    X <- matrix(sample(c(-1, 1), n * p, replace = TRUE, prob = c(1 - maf, maf)),
                nrow = n, ncol = p)
  } else {
    X <- matrix(rbinom(n * p, size = 2, prob = maf), nrow = n, ncol = p)
    X[X == 0] <- -1; X[X == 1] <- 0; X[X == 2] <- 1
  }
  
  colnames(X) <- paste0("SNP", seq_len(p))
  rownames(X) <- paste0("Ind", seq_len(n))
  X
}

#' 由 PVE 反推效应量标准差
#' RIL: Var(X)=4*maf*(1-maf); F2: Var(X)=2*maf*(1-maf)
sd_from_pve <- function(pve, maf = 0.3, population = c("RIL", "F2")) {
  population <- match.arg(population)
  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  sqrt(pve / var_x)
}

# ------------------------------------------------------------------------------
# 2. 核心新增：方差标准化函数（修正矩阵顺序）
# ------------------------------------------------------------------------------

#' 计算使表型方差为 1 的残差协方差矩阵 Sigma_E
#' 
#' 数学原理：
#'   M1 (Tau=0): Sigma_Y = Sigma_G + Sigma_E
#'   M2 (Tau≠0): Y = Y_base %*% A, A = (I - t(Tau))^{-1}
#'               Sigma_Y = A^T %*% Sigma_base %*% A
#'                       = (I-Tau)^{-1} %*% Sigma_base %*% (I-t(Tau))^{-1}
#' 
#' @param B p×m 效应矩阵
#' @param Tau m×m 因果矩阵
#' @param var_x 基因型方差（RIL: 4*maf*(1-maf); F2: 2*maf*(1-maf)）
#' @param rho 残差相关系数
#' @param max_iter 最大迭代次数
#' @param tol 收敛容差
compute_residual_variance <- function(B, Tau, var_x, rho = 0, 
                                        max_iter = 200, tol = 1e-8) {
  m <- ncol(B)
  Sigma_G <- var_x * crossprod(B)  # m×m，基因型贡献的协方差
  
  # M1: 直接解析解
  if (all(abs(Tau) < 1e-12)) {
    diag_vals <- pmax(1e-6, 1 - diag(Sigma_G))
    Sigma_E <- diag(diag_vals)
    if (abs(rho) > 1e-10 && m > 1) {
      s <- sqrt(diag_vals)
      Sigma_E <- matrix(rho * outer(s, s), m, m)
      diag(Sigma_E) <- diag_vals
    }
    return(Sigma_E)
  }
  
  # M2: 迭代求解
  I_m <- diag(m)
  A <- solve(I_m - t(Tau))  # A = (I - t(Tau))^{-1}
  
  # 初始 Sigma_E
  Sigma_E <- diag(pmax(1e-6, 1 - diag(Sigma_G)))
  if (abs(rho) > 1e-10 && m > 1) {
    s <- sqrt(diag(Sigma_E))
    Sigma_E <- matrix(rho * outer(s, s), m, m)
    diag(Sigma_E) <- s^2             # ✅ 对角线还原为 s²
  }
  
  for (iter in seq_len(max_iter)) {
    Sigma_base <- Sigma_G + Sigma_E
    # ✅ 修正：Sigma_Y = t(A) %*% Sigma_base %*% A
    #    t(A) = (I - Tau)^{-1}
    Sigma_Y <- t(A) %*% Sigma_base %*% A
    
    diff <- 1 - diag(Sigma_Y)
    diag(Sigma_E) <- diag(Sigma_E) + diff
    
    if (abs(rho) > 1e-10 && m > 1) {
      s <- sqrt(pmax(1e-6, diag(Sigma_E)))
      for (i in 1:m) for (j in 1:m) if (i != j) Sigma_E[i, j] <- rho * s[i] * s[j]
    }
    
    if (max(abs(diff)) < tol) break
    if (iter == max_iter) warning("Variance standardization did not converge")
  }
  
  Sigma_E
}

# ------------------------------------------------------------------------------
# 3. 表型生成模型
# ------------------------------------------------------------------------------

#' M1: 基础模型（无因果链）
sim_phenotype_M1 <- function(X, B, Sigma_E,
                              epi_pairs = NULL, epi_effects = NULL,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X); m <- ncol(B)
  
  G <- X %*% B
  if (!is.null(epi_pairs) && !is.null(epi_effects)) {
    for (k in seq_along(epi_pairs)) {
      pair <- epi_pairs[[k]]
      G <- G + (X[, pair[1]] * X[, pair[2]]) %*% epi_effects[k, , drop = FALSE]
    }
  }
  
  eps <- MASS::mvrnorm(n, mu = rep(0, m), Sigma = Sigma_E)
  Y <- G + eps
  colnames(Y) <- paste0("Trait", 1:m)
  list(Y = Y, G = G, eps = eps)
}

#' M2: 因果链模型
#' Y = (X*B + eps) %*% solve(I - t(Tau))
#' 协方差传播: Sigma_Y = (I-Tau)^{-1} %*% Sigma_base %*% (I-t(Tau))^{-1}
sim_phenotype_M2 <- function(X, B, Tau, Sigma_E,
                              epi_pairs = NULL, epi_effects = NULL,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- ncol(B)
  base <- sim_phenotype_M1(X, B, Sigma_E, epi_pairs, epi_effects, seed = NULL)
  
  I_T_inv <- solve(diag(m) - t(Tau))
  Y <- base$Y %*% I_T_inv
  
  list(Y = Y, G = base$G, eps = base$eps, Tau = Tau, I_T_inv = I_T_inv)
}

# ------------------------------------------------------------------------------
# 4. 三数据集生成函数（对应 Figure 1/2）
# ------------------------------------------------------------------------------

#' Dataset I: 无因果 (τ=0, ρ=0)，验证 Class 1 / Class 2
generate_dataset_I <- function(n, p = 1000,
                                class1_pve = 0.01, class2_pve = 0.01,
                                maf = 0.3, population = "RIL", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf, population)
  B <- matrix(0, p, 2)
  
  # Class 1: 位点 1-2，仅影响 Trait A
  B[1:2, 1] <- rnorm(2, 0, sd_from_pve(class1_pve, maf, population))
  # Class 2: 位点 3-4，独立影响 A 和 B
  B[3:4, 1] <- rnorm(2, 0, sd_from_pve(class2_pve, maf, population))
  B[3:4, 2] <- rnorm(2, 0, sd_from_pve(class2_pve, maf, population))
  
  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, matrix(0, 2, 2), var_x, rho = 0)
  out <- sim_phenotype_M1(X, B, Sigma_E)
  
  truth <- data.frame(
    SNP = colnames(X),
    class = c(rep("class1", 2), rep("class2", 2), rep("null", p - 4)),
    beta_A = B[, 1], beta_B = B[, 2],
    is_IV = FALSE,
    stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, dataset = "I", seed = seed))
}

#' Dataset II: 单向因果 A→B (τ=0.3)，验证 Class 1 / Class 3 / Class 4
#' 含 10 个 IV 辅助位点（仅影响 A，不纳入分类）
generate_dataset_II <- function(n, p = 1000,
                                 class1_pve = 0.01,
                                 class3_pve_A = 0.02,
                                 class4_pve_A = 0.02, class4_pve_B = 0.003,
                                 tau = 0.3, iv_pve = 0.02,
                                 maf = 0.3, population = "RIL", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf, population)
  B <- matrix(0, p, 2)
  
  # Class 1: 位点 1-2，仅影响 B
  B[1:2, 2] <- rnorm(2, 0, sd_from_pve(class1_pve, maf, population))
  # Class 3: 位点 3-4，仅影响 A（B 完全由 τ 中介）
  B[3:4, 1] <- rnorm(2, 0, sd_from_pve(class3_pve_A, maf, population))
  # Class 4: 位点 5-6，影响 A(强) + B(弱)
  B[5:6, 1] <- rnorm(2, 0, sd_from_pve(class4_pve_A, maf, population))
  B[5:6, 2] <- rnorm(2, 0, sd_from_pve(class4_pve_B, maf, population))
  # IV 辅助: 位点 7-16，仅影响 A，PVE=2%
  B[7:16, 1] <- rnorm(10, 0, sd_from_pve(iv_pve, maf, population))
  
  Tau <- matrix(0, 2, 2); Tau[2, 1] <- tau  # A -> B
  
  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, Tau, var_x, rho = 0)
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E)
  
  truth <- data.frame(
    SNP = colnames(X),
    class = c(rep("class1", 2), rep("class3", 2), rep("class4", 2),
              rep("IV_A", 10), rep("null", p - 16)),
    beta_A = B[, 1], beta_B = B[, 2],
    tau = c(rep(0, 6), rep(tau, 10), rep(0, p - 16)),
    is_IV = c(rep(FALSE, 6), rep(TRUE, 10), rep(FALSE, p - 16)),
    stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, dataset = "II", Tau = Tau, seed = seed))
}

#' Dataset III: 双向因果 A↔B (τ_AB=τ_BA=0.2)，验证 Class 5
#' 含 IV-仅A(10) + IV-仅B(10)
generate_dataset_III <- function(n, p = 1000,
                                    class5_pve = 0.01,
                                    tau_AB = 0.2, tau_BA = 0.2,
                                    iv_pve = 0.02,
                                    maf = 0.3, population = "RIL", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf, population)
  B <- matrix(0, p, 2)
  
  # Class 5: 位点 1-2，直接影响 A 和 B（确保 n_cond=2）
  B[1:2, 1] <- rnorm(2, 0, sd_from_pve(class5_pve, maf, population))
  B[1:2, 2] <- rnorm(2, 0, sd_from_pve(class5_pve, maf, population))
  # IV-仅A: 位点 3-12
  B[3:12, 1] <- rnorm(10, 0, sd_from_pve(iv_pve, maf, population))
  # IV-仅B: 位点 13-22
  B[13:22, 2] <- rnorm(10, 0, sd_from_pve(iv_pve, maf, population))
  
  Tau <- matrix(0, 2, 2)
  Tau[2, 1] <- tau_AB  # A -> B
  Tau[1, 2] <- tau_BA  # B -> A
  
  var_x <- if (population == "RIL") 4 * maf * (1 - maf) else 2 * maf * (1 - maf)
  Sigma_E <- compute_residual_variance(B, Tau, var_x, rho = 0)
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E)
  
  truth <- data.frame(
    SNP = colnames(X),
    class = c(rep("class5", 2), rep("IV_A", 10), rep("IV_B", 10),
              rep("null", p - 22)),
    beta_A = B[, 1], beta_B = B[, 2],
    tau_AB = c(rep(tau_AB, 22), rep(0, p - 22)),
    tau_BA = c(rep(tau_BA, 22), rep(0, p - 22)),
    is_IV = c(rep(FALSE, 2), rep(TRUE, 10), rep(TRUE, 10), rep(FALSE, p - 22)),
    stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, dataset = "III", Tau = Tau, seed = seed))
}

# ------------------------------------------------------------------------------
# 5. 统一入口
# ------------------------------------------------------------------------------

generate_condped <- function(dataset = c("I", "II", "III"), ...) {
  dataset <- match.arg(dataset)
  switch(dataset,
    I   = generate_dataset_I(...),
    II  = generate_dataset_II(...),
    III = generate_dataset_III(...)
  )
}

# ------------------------------------------------------------------------------
# 6. 验证脚本（跑完确认方差标准化正确）
# ------------------------------------------------------------------------------
# 
# cat("=== 验证 Dataset II (tau=0.3) ===\n")
# set.seed(123)
# dat2 <- generate_dataset_II(n = 1000, p = 1000)
# cat("Trait variances:", round(apply(dat2$Y, 2, var), 4), "\n")
# cat("Expected: ~1.0, ~1.0\n\n")
# 
# cat("=== 验证 Dataset III (tau_AB=0.2, tau_BA=0.2) ===\n")
# set.seed(456)
# dat3 <- generate_dataset_III(n = 1000, p = 1000)
# cat("Trait variances:", round(apply(dat3$Y, 2, var), 4), "\n")
# cat("Expected: ~1.0, ~1.0\n\n")
# 
# cat("=== 验证 Dataset I (tau=0) ===\n")
# set.seed(789)
# dat1 <- generate_dataset_I(n = 1000, p = 1000)
# cat("Trait variances:", round(apply(dat1$Y, 2, var), 4), "\n")
# cat("Expected: ~1.0, ~1.0\n")