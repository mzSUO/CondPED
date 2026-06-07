# ==============================================================================
# CondPED 模拟数据生成系统（修正版）
# ==============================================================================
# 修正内容：
#   1. M1/M2 使用期望方差（deterministic）代替样本方差控制残差，避免小样本波动；
#   2. M2 浮点比较改为安全阈值；
#   3. 所有 generate_* 函数将 maf 传入底层模型。
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 基础辅助函数
# ------------------------------------------------------------------------------

#' 模拟 SNP 基因型矩阵
#'
#' 基于二项分布生成 Hardy-Weinberg 平衡下的基因型剂量编码。
#'
#' @param n 样本量。
#' @param p SNP 位点数。
#' @param maf 次等位基因频率。默认 0.3。
#' @param seed 随机种子。
#' @return 整数矩阵，维度为 n x p，元素取值 0, 1, 2。
#' @export
sim_genotype <- function(n, p, maf = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(n > 0, p > 0, maf > 0, maf < 1)

  X <- matrix(stats::rbinom(n * p, size = 2, prob = maf), nrow = n, ncol = p)
  colnames(X) <- paste0("SNP", seq_len(p))
  rownames(X) <- paste0("Ind", seq_len(n))
  X
}

#' 根据 PVE 计算效应大小的标准差
#'
#' 公式：sd = sqrt(PVE / Var(X))，Var(X) = 2 * maf * (1 - maf)。
#'
#' @param pve 目标表型方差解释率。
#' @param maf 次等位基因频率。默认 0.3。
#' @return 效应大小的标准差。
#' @keywords internal
sd_from_pve <- function(pve, maf = 0.3) {
  var_x <- 2 * maf * (1 - maf)
  sqrt(pve / var_x)
}

#' 校验表型方差标准化
#'
#' @param Y 表型矩阵，n x m。
#' @param tol 容差阈值。默认 0.15。
#' @return 逻辑值。
#' @keywords internal
check_variance_standardization <- function(Y, tol = 0.15) {
  vars <- apply(Y, 2, stats::var)
  all(abs(vars - 1) < tol)
}

# ------------------------------------------------------------------------------
# 2. 核心表型生成模型（修正版：期望方差控制）
# ------------------------------------------------------------------------------

#' 基础模型 M1
#'
#' 生成式：Y = X * B + E
#' 使用期望方差（deterministic）控制残差，避免小样本随机波动导致截断失真。
#'
#' @param X 基因型矩阵，n x p。
#' @param B SNP 效应矩阵，p x m。
#' @param Sigma_E m x m 残差协方差矩阵（相关结构保留，方差被动态调整）。
#' @param maf 次等位基因频率（用于计算期望遗传方差）。默认 0.3。
#' @param seed 随机种子。
#' @export
sim_phenotype_M1 <- function(X, B, Sigma_E, maf = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  m <- ncol(B)

  stopifnot(ncol(X) == nrow(B))
  stopifnot(dim(Sigma_E) == c(m, m))

  # 1. 遗传效应
  G <- X %*% B

  # 2. 期望方差控制（修正：用 deterministic 期望代替样本方差）
  expected_var_G <- colSums(B^2) * 2 * maf * (1 - maf)
  sigma_eps <- pmax(sqrt(1 - expected_var_G), 0.1)

  D <- diag(sigma_eps)
  Cor_E <- stats::cov2cor(Sigma_E)
  Sigma_E_actual <- D %*% Cor_E %*% D

  E <- mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = Sigma_E_actual)

  # 3. 基础表型
  Y <- G + E

  # 4. 校验
  vars <- apply(Y, 2, stats::var)
  if (any(abs(vars - 1) > 0.15)) {
    warning(sprintf("M1: 性状方差偏离 1: %s", paste(round(vars, 3), collapse = ", ")))
  }

  list(
    Y = Y,
    G = G,
    E = E,
    B = B,
    Sigma_E = Sigma_E_actual,
    expected_var_G = expected_var_G,
    params = list(model = "M1", target_var = 1, actual_var = as.numeric(vars))
  )
}

#' 因果链模型 M2
#'
#' 生成式：
#'   Y_base = X * B + E
#'   Y = Y_base %*% t(solve(I_m - Tau))
#'
#' 方差控制在因果变换之前完成，保留真实 Tau 尺度。
#'
#' @param X 基因型矩阵。
#' @param B SNP 效应矩阵，p x m。
#' @param Tau m x m 因果路径矩阵。
#' @param Sigma_E m x m 残差协方差矩阵。
#' @param maf 次等位基因频率。默认 0.3。
#' @param seed 随机种子。
#' @export
sim_phenotype_M2 <- function(X, B, Tau, Sigma_E, maf = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  m <- ncol(B)

  stopifnot(ncol(X) == nrow(B))
  stopifnot(dim(Tau) == c(m, m))
  stopifnot(dim(Sigma_E) == c(m, m))

  # 1. 遗传效应
  G <- X %*% B

  # 2. 期望方差控制（修正）
  expected_var_G <- colSums(B^2) * 2 * maf * (1 - maf)
  sigma_eps <- pmax(sqrt(1 - expected_var_G), 0.1)

  D <- diag(sigma_eps)
  Cor_E <- stats::cov2cor(Sigma_E)
  Sigma_E_actual <- D %*% Cor_E %*% D

  E <- mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = Sigma_E_actual)

  # 3. 基础表型
  Y_base <- G + E

  # 4. 因果链变换
  I_T <- diag(m) - Tau
  I_T_inv <- solve(I_T)
  Y <- Y_base %*% t(I_T_inv)

  # 5. 校验
  vars <- apply(Y, 2, stats::var)
  if (any(abs(vars - 1) > 0.5)) {
    warning(sprintf(
      "M2: 因果变换后性状方差: %s。tau 可能导致尺度变化，属正常现象。",
      paste(round(vars, 3), collapse = ", ")
    ))
  }

  # 6. 内部诊断：单向因果时验证观测 tau（修正：浮点安全比较）
  tau_observed <- NULL
  if (m == 2 && abs(Tau[2, 1]) > 1e-10 && abs(Tau[1, 2]) < 1e-10) {
    fit_check <- stats::lm(Y[, 2] ~ Y[, 1])
    tau_obs <- unname(stats::coef(fit_check)[2])
    tau_true <- Tau[2, 1]
    tau_observed <- tau_obs
    if (abs(tau_obs - tau_true) > 0.15) {
      warning(sprintf(
        "M2: 观测 tau (%.3f) 偏离真实 tau (%.3f)。请检查方差控制逻辑。",
        tau_obs, tau_true
      ))
    }
  }

  list(
    Y = Y,
    Y_base = Y_base,
    G = G,
    E = E,
    B = B,
    Tau = Tau,
    Sigma_E = Sigma_E_actual,
    expected_var_G = expected_var_G,
    params = list(
      model = "M2",
      target_var = 1,
      actual_var = as.numeric(vars),
      tau_observed = tau_observed
    )
  )
}

# ------------------------------------------------------------------------------
# 3. 五类遗传机制生成器（修正：传入 maf）
# ------------------------------------------------------------------------------

#' @export
generate_class1 <- function(n, p = 200, pve = 0.01, maf = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf))
  out <- sim_phenotype_M1(X, B, Sigma_E = diag(2), maf = maf)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class1", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class1", seed = seed))
}

#' @export
generate_class3 <- function(n, p = 200, pve = 0.01, maf = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf))
  B[1:30, 2] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf))
  out <- sim_phenotype_M1(X, B, Sigma_E = diag(2), maf = maf)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class3", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class3", seed = seed))
}

#' @export
generate_class4 <- function(n, p = 200, pve_A = 0.02, tau = 0.3, maf = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve_A, maf))
  Tau <- matrix(0, 2, 2); Tau[2, 1] <- tau
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E = diag(2), maf = maf)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class4", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], tau = tau, stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class4", seed = seed))
}

#' @export
generate_class5 <- function(n, p = 200, pve_A = 0.02, pve_B = 0.003,
                             tau = 0.3, maf = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve_A, maf))
  B[1:30, 2] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve_B, maf))
  Tau <- matrix(0, 2, 2); Tau[2, 1] <- tau
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E = diag(2), maf = maf)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class5", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], tau = tau, stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class5", seed = seed))
}

#' @export
generate_class6 <- function(n, p = 200, pve = 0.01,
                             tau_AB = 0.2, tau_BA = 0.2,
                             maf = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf))
  Tau <- matrix(0, 2, 2)
  Tau[2, 1] <- tau_AB
  Tau[1, 2] <- tau_BA
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E = diag(2), maf = maf)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class6", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2],
    tau_AB = tau_AB, tau_BA = tau_BA, stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class6", seed = seed))
}

# ------------------------------------------------------------------------------
# 4. Class 2 压力测试生成器（修正：传入 maf）
# ------------------------------------------------------------------------------

#' @export
generate_covariance_stress <- function(n, p = 200, pve = 0.01,
                                        rho = 0.5, maf = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf))
  Sigma_E <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
  out <- sim_phenotype_M1(X, B, Sigma_E, maf = maf)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class1", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], rho = rho, stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "stress_test", seed = seed))
}

# ------------------------------------------------------------------------------
# 5. 统一入口函数
# ------------------------------------------------------------------------------

#' @export
generate_condped <- function(scenario = c("class1", "class3", "class4", "class5", "class6", "stress_test"),
                              n, p = 200, ...) {
  scenario <- match.arg(scenario)
  switch(scenario,
    class1 = generate_class1(n = n, p = p, ...),
    class3 = generate_class3(n = n, p = p, ...),
    class4 = generate_class4(n = n, p = p, ...),
    class5 = generate_class5(n = n, p = p, ...),
    class6 = generate_class6(n = n, p = p, ...),
    stress_test = generate_covariance_stress(n = n, p = p, ...)
  )
}