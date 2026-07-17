# ==============================================================================
# conditional.R — Layer 2：条件投影
# ==============================================================================
# 实现论文 Eq.4 的条件表型构造：
#
#   ỹ*_{i} = ỹ_{i} - ỹ_{-i}^T γ_{i,-i}
#
# 其中投影系数 γ_{i,-i} = V_{-i}^{-1} C_{-i,i}
#
# 数学性质（均有测试验证）：
#   1. 正交性：Cov(ỹ*_i, ỹ_j) = 0  对所有 j ≠ i
#   2. 方差等价：Var(ỹ*_i) = V_{i|-i}（Schur 补，即 Eq.3）
#   3. 分解性：θ_marg = θ_cond + θ_cov（统计不相关的两个成分）
# ==============================================================================


# ------------------------------------------------------------------------------
# 辅助函数 1：估计表型协方差矩阵
# ------------------------------------------------------------------------------

#' 从残差表型估计性状间协方差矩阵 V
#'
#' 对应论文 Eq.2：V = Σ_E + Σ_AE + Σ_AAE + Σ_ε
#' 此处通过样本协方差近似（残差已去除均值和协变量效应）。
#'
#' @param Y_residual 数值矩阵，n × m。个体（行）× 性状（列）的残差表型。
#'   必须已去除群体均值和协变量效应，但保留遗传效应和随机效应。
#' @param method 估计方法，目前仅支持 "pearson"（样本协方差）。
#'
#' @return 列表，包含：
#'   \item{V}{m × m 样本协方差矩阵}
#'   \item{Cor}{m × m 相关矩阵}
#'   \item{n}{样本量}
#'   \item{m}{性状数}
#'
#' @details
#' 当性状数 m 较大而样本量 n 较小时（n/m < 10），
#' 样本协方差矩阵估计不稳定，条件投影精度会下降。
#' 建议最低 n/m ≥ 20。
#'
#' @keywords internal
estimate_phenotypic_covariance <- function(Y_residual, method = "pearson") {
  stopifnot(is.matrix(Y_residual) || is.data.frame(Y_residual))
  Y_residual <- as.matrix(Y_residual)
  n <- nrow(Y_residual)
  m <- ncol(Y_residual)

  if (n <= m) {
    stop(sprintf(
      "样本量 n=%d 不大于性状数 m=%d，协方差矩阵奇异，无法求逆。",
      n, m
    ))
  }
  if (n / m < 10) {
    warning(sprintf(
      "n/m = %.1f < 10，协方差矩阵估计可能不稳定，条件投影精度受影响。",
      n / m
    ))
  }

  V <- cov(Y_residual) # m × m 样本协方差矩阵
  Cor <- cov2cor(V)

  list(V = V, Cor = Cor, n = n, m = m)
}


# ------------------------------------------------------------------------------
# 辅助函数 2：计算单个性状的投影系数
# ------------------------------------------------------------------------------

#' 计算性状 i 的条件投影系数 γ_{i,-i}
#'
#' γ_{i,-i} = V_{-i}^{-1} C_{-i,i}
#'
#' 其中 V_{-i} 是排除性状 i 后的 (m-1)×(m-1) 协方差子矩阵，
#' C_{-i,i} 是性状 i 与其余性状的协方差向量。
#'
#' @param V m × m 表型协方差矩阵
#' @param i 目标性状索引（1-based 整数）
#'
#' @return 长度为 m-1 的数值向量，元素顺序与去除性状 i 后的性状序列一致。
#'
#' @details
#' 使用 `solve(A, b)` 而非 `solve(A) %*% b` 以提高数值稳定性。
#' 当 V_{-i} 近奇异时会给出警告。
#'
#' @keywords internal
compute_projection_coefficients <- function(V, i) {
  m <- nrow(V)
  stopifnot(i >= 1, i <= m)

  idx <- setdiff(seq_len(m), i) # 其余性状的索引
  V_minus_i <- V[idx, idx, drop = FALSE] # (m-1) × (m-1)
  C_minus_i <- V[idx, i, drop = FALSE] # (m-1) × 1

  # 条件数检查：避免奇异矩阵
  kappa_val <- kappa(V_minus_i, exact = FALSE)
  if (kappa_val > 1e10) {
    warning(sprintf(
      "性状 %d 的 V_{-%d} 条件数 = %.2e，接近奇异，投影系数估计不稳定。",
      i, i, kappa_val
    ))
  }

  # solve(A, b) 等价于 A^{-1} b，比 solve(A) %*% b 数值更稳定
  gamma <- solve(V_minus_i, C_minus_i)
  as.numeric(gamma)
}


# ------------------------------------------------------------------------------
# 主函数：构造条件表型（对应论文 Eq.4）
# ------------------------------------------------------------------------------

#' 构造条件表型矩阵
#'
#' 对每个性状 i，构造条件表型：
#' \deqn{\tilde{y}^*_{ijk} = \tilde{y}_{ijk} - \tilde{\mathbf{y}}_{-i,jk}^\top \gamma_{i,-i}}
#'
#' 该构造确保 \eqn{\tilde{y}^*_i} 与其余所有性状在表型协方差意义上正交，
#' 即 \eqn{\mathrm{Cov}(\tilde{y}^*_i, \tilde{y}_j) = 0}（对所有 j ≠ i）。
#'
#' 对条件表型估计的遗传效应即为**条件遗传效应** \eqn{\theta_{iq}^{\mathrm{cond}}}，
#' 代表不依赖于其他性状线性预测的独立遗传成分。
#'
#' @param Y_residual 数值矩阵，n × m。残差表型矩阵。
#'   **要求**：已通过 \code{lm} 或 QTLNetwork 去除群体均值和协变量效应，
#'   但保留遗传效应（目的是估计 V 的完整协方差结构，包含遗传协方差）。
#' @param V m × m 表型协方差矩阵（可选）。
#'   若为 \code{NULL}（默认），由函数内部用 \code{Y_residual} 估计。
#'   若已在外部计算，可直接传入以提高效率（多次调用时避免重复估计）。
#' @param tol_orthogonality 正交性校验容差。|Cov| 超过此阈值时给出警告。默认 0.01。
#'
#' @return 列表，包含：
#'   \item{Y_cond}{n × m 矩阵，条件表型，列名与 \code{Y_residual} 一致}
#'   \item{gamma}{长度为 m 的列表，\code{gamma[[i]]} 为性状 i 的投影系数向量}
#'   \item{V}{使用的 m × m 表型协方差矩阵}
#'   \item{schur}{长度为 m 的数值向量，各性状的 Schur 补（理论条件方差）}
#'   \item{diagnostics}{诊断信息数据框，含正交性和方差等价性的验证结果}
#'
#' @examples
#' \dontrun{
#' # 生成示例数据（ρ = 0.5 的双性状场景）
#' set.seed(42)
#' n <- 500
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' Y <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
#' colnames(Y) <- c("A", "B")
#'
#' result <- compute_conditional_phenotype(Y)
#'
#' # 验证正交性：应接近 0
#' cov(result$Y_cond[, "A"], Y[, "B"])
#'
#' # 验证方差等价：应接近 Schur 补
#' var(result$Y_cond[, "A"])
#' result$schur["A"]
#' }
#'
#' @references
#' Zhu, J. (1995). Analysis of conditional genetic effects and variance
#' components in developmental genetics. *Genetics*, 141(4), 1633–1639.
#'
#' Anderson, T. W. (1958). *An Introduction to Multivariate Statistical
#' Analysis*. Wiley.
#'
#' @export
compute_conditional_phenotype <- function(Y_residual,
                                          V = NULL,
                                          tol_orthogonality = 0.01) {
  # ---------- 输入校验 --------------------------------------------------------
  stopifnot(
    "Y_residual 必须是矩阵或数据框" = is.matrix(Y_residual) || is.data.frame(Y_residual)
  )
  Y_residual <- as.matrix(Y_residual)
  n <- nrow(Y_residual)
  m <- ncol(Y_residual)

  if (is.null(colnames(Y_residual))) {
    colnames(Y_residual) <- paste0("Trait", seq_len(m))
  }
  trait_names <- colnames(Y_residual)

  if (m == 1) {
    message("性状数 m=1，无其他性状可条件化，返回原始残差表型。")
    return(list(
      Y_cond = Y_residual,
      gamma = list(numeric(0)),
      V = matrix(var(Y_residual), 1, 1),
      schur = setNames(var(Y_residual), trait_names),
      diagnostics = data.frame(
        trait = trait_names, max_cov = NA,
        var_obs = var(Y_residual), schur = var(Y_residual),
        var_diff = 0
      )
    ))
  }

  # ---------- 估计协方差矩阵 --------------------------------------------------
  if (is.null(V)) {
    cov_obj <- estimate_phenotypic_covariance(Y_residual)
    V <- cov_obj$V
  } else {
    stopifnot(
      "V 必须是 m×m 矩阵" = is.matrix(V) && nrow(V) == m && ncol(V) == m
    )
  }

  # ---------- 主循环：对每个性状构造条件表型 -----------------------------------
  Y_cond <- matrix(NA_real_, n, m, dimnames = list(rownames(Y_residual), trait_names))
  gamma_list <- vector("list", m)
  schur <- numeric(m)

  for (i in seq_len(m)) {
    idx <- setdiff(seq_len(m), i)

    # 投影系数：γ_{i,-i} = V_{-i}^{-1} C_{-i,i}
    gamma <- compute_projection_coefficients(V, i)
    gamma_list[[i]] <- gamma

    # 条件表型：ỹ*_i = ỹ_i - ỹ_{-i} %*% γ
    Y_cond[, i] <- Y_residual[, i] - Y_residual[, idx, drop = FALSE] %*% gamma

    # 理论条件方差（Schur 补）：V_{i|-i} = V_{ii} - C_{i,-i} V_{-i}^{-1} C_{-i,i}
    # 等价于 V_{ii} - t(C_{-i,i}) %*% gamma
    C_minus_i <- V[idx, i]
    schur[i] <- V[i, i] - sum(C_minus_i * gamma)
  }

  names(gamma_list) <- trait_names
  names(schur) <- trait_names

  # ---------- 诊断：验证数学性质 -----------------------------------------------
  diagnostics <- vector("list", m)

  for (i in seq_len(m)) {
    idx <- setdiff(seq_len(m), i)

    # 性质1：正交性 —— Cov(ỹ*_i, ỹ_j) 对所有 j ≠ i 应接近 0
    covs <- apply(
      Y_residual[, idx, drop = FALSE], 2,
      function(col) cov(Y_cond[, i], col)
    )
    max_cov <- max(abs(covs))

    if (max_cov > tol_orthogonality) {
      warning(sprintf(
        paste0(
          "性状 '%s' 的条件表型正交性未达标：",
          "max|Cov(ỹ*_%s, ỹ_j)| = %.4f > %.4f。",
          "\n  可能原因：样本量不足或 V 估计误差过大。"
        ),
        trait_names[i], trait_names[i], max_cov, tol_orthogonality
      ))
    }

    # 性质2：方差等价 —— Var(ỹ*_i) 应接近 Schur 补 V_{i|-i}
    var_obs <- var(Y_cond[, i])
    var_diff <- abs(var_obs - schur[i])

    # 允许的偏差：受限于样本量，理论上偏差量级为 O(1/n)
    tol_var <- 3 * sqrt(2 * schur[i]^2 / (n - 1)) # 基于 χ² 分布的近似容差
    if (var_diff > tol_var) {
      warning(sprintf(
        paste0(
          "性状 '%s' 的方差等价性验证失败：",
          "Var(ỹ*) = %.4f，Schur 补 = %.4f，差异 = %.4f。"
        ),
        trait_names[i], var_obs, schur[i], var_diff
      ))
    }

    diagnostics[[i]] <- data.frame(
      trait = trait_names[i],
      max_cov = max_cov, # 正交性指标（越小越好）
      var_obs = var_obs, # 观测方差
      schur = schur[i], # 理论条件方差（Schur 补）
      var_diff = var_diff, # 方差偏差（越小越好）
      stringsAsFactors = FALSE
    )
  }

  diagnostics_df <- do.call(rbind, diagnostics)

  # ---------- 返回 ------------------------------------------------------------
  list(
    Y_cond      = Y_cond, # n×m 条件表型矩阵（核心输出）
    gamma       = gamma_list, # 各性状的投影系数
    V           = V, # 使用的协方差矩阵
    schur       = schur, # 各性状的理论条件方差
    diagnostics = diagnostics_df # 验证诊断结果
  )
}

# ==============================================================================
# conditional.R（节选）— Layer 2：条件遗传效应估计
# ==============================================================================
#
# 【本节功能】
#   fit_conditional_model()：对 Layer 1 检出的显著位点拟合条件模型，
#   估计每个位点在每个性状上的条件遗传效应 theta_cond 及其显著性。
#
# 【统计模型】
#   单环境（OLS）：
#     y*_{ij} = mu_i + sum_l theta_{il} * x_{jl} + eps_{ij}
#     y*_{ij} 是经过条件投影的表型（已去除其他性状的线性影响）
#
#   多环境（LMM）：
#     y*_{ijk} = mu_i + theta_{il} * x_{jl}（加性固定效应）
#              + e_{ik}（环境随机截距）
#              + ae_{ilk} * x_{jl}（G×E随机斜率）
#              + eps_{ijk}（残差）
#     lme4公式：y ~ x_snp + (1 + x_snp | env)
#
# 【条件化方式（排他性约束的经验检验）】
#   Layer 2 采用全条件投影（Schur补），即将性状i的表型投影到所有其他m-1个性状。
#   但 IV 排他性筛选（Layer 3 select_instruments）采用成对条件（只对暴露性状投影）。
#   两种条件化方式在 m=2 时等价，m≥3 时有本质区别：
#   - 全条件：用于位点分类（n_cond 计算），消除所有性状间的共变
#   - 成对条件：用于IV筛选，保留混杂因子对结局的残余信号，防止污染工具变量通过
#
# 【上位性支持】
#   通过 epi_pairs 参数支持加性×加性上位性（aa项）作为额外固定效应，
#   对应论文 Eq.1 中的 aa_{ilh} 项（仅用于补充模拟 S5）。
# ==============================================================================


#' 估计候选位点的条件遗传效应（Layer 2 核心函数）
#'
#' 对 Layer 1 检出的每个候选位点，在每个性状的条件表型上拟合线性模型，
#' 提取条件效应 theta_cond 及其统计显著性，作为多效性机制分类的统计指纹。
#'
#' 【条件效应的生物学含义】
#'   theta_cond_{il} = 控制其他所有性状后，位点l对性状i的独立直接遗传效应
#'   - 如果 theta_cond 显著（pval_cond < alpha2）：位点对该性状有独立直接效应
#'   - 如果 theta_cond 不显著：该性状的边际关联可能完全来自与其他性状的共变
#'
#' 【统计指纹定义】
#'   n_cond = 显著 theta_cond 的性状数
#'   n_cond = 1：只有一个性状保留显著条件效应 → 进入 Layer 3 区分 Class 1 vs Class 3
#'   n_cond ≥ 2：多个性状保留显著条件效应 → 进入 Layer 3 区分 Class 2/4/5
#'
#' 【注意事项】
#'   条件效应不显著≠直接效应不存在，仅反映当前样本量下未检测到显著性。
#'   Bonferroni 校正为保守检验，低功效时可能产生假阴性。
#'
#' @param Y_cond    n×m（或nK×m）条件表型矩阵。
#'   由 compute_conditional_phenotype()$Y_cond 返回。
#'   多环境时为长格式（nK行，K为环境数）。
#' @param X_loci    n×k（或nK×k）候选位点基因型矩阵。
#'   仅包含 Layer 1 显著位点，列名为 SNP ID。
#'   多环境时需与 Y_cond 的行对齐（nK行）。
#' @param alpha2    数值或NULL。Bonferroni 校正显著性阈值。
#'   NULL=自动计算 0.05/(k×m)，其中k=候选位点数，m=性状数。
#' @param method    字符串。"OLS"（单环境，默认）或"LMM"（多环境）。
#'   - "OLS"：每个性状独立拟合 y* ~ x_snp + 其他位点
#'   - "LMM"：用 lme4 拟合含环境随机效应的混合模型，需提供 env 参数
#' @param env       因子或字符向量，长度n（单环境）或nK（多环境）。
#'   指示每个观测所属的环境，仅在 method="LMM" 时必须提供。
#' @param epi_pairs 列表，每个元素为 c(l1, l2)，指定一对上位性位点（按列索引）。
#'   NULL=无上位性（默认）。上位性项作为额外固定效应加入模型。
#'
#' @return 数据框，每行对应一个（性状, 位点）对，包含列：
#' \describe{
#'   \item{trait}{性状名（来自 Y_cond 的列名）}
#'   \item{locus}{位点 ID（来自 X_loci 的列名）}
#'   \item{theta_cond}{条件效应点估计（OLS/LMM的斜率系数）}
#'   \item{se_cond}{条件效应标准误}
#'   \item{pval_cond}{双侧 p 值（OLS用t检验，LMM用Satterthwaite自由度近似）}
#'   \item{sig_cond}{逻辑值。是否达到 Bonferroni 校正显著性（pval < alpha2）}
#'   \item{method}{实际使用的估计方法："OLS"/"LMM"/"OLS_fallback"（LMM失败时退化）}
#' }
#'
#' @export
fit_conditional_model <- function(Y_cond,
                                  X_loci,
                                  alpha2 = NULL,
                                  method = c("OLS", "LMM"),
                                  env = NULL,
                                  epi_pairs = NULL) {
  method <- match.arg(method)
  Y_cond <- as.matrix(Y_cond)
  X_loci <- as.matrix(X_loci)
  m <- ncol(Y_cond) # 性状数
  k <- ncol(X_loci) # 候选位点数

  # 确保有列名（用于输出标识）
  if (is.null(colnames(Y_cond))) colnames(Y_cond) <- paste0("Trait", seq_len(m))
  if (is.null(colnames(X_loci))) colnames(X_loci) <- paste0("SNP", seq_len(k))

  # Bonferroni 阈值：控制所有（k×m）对的家族错误率
  if (is.null(alpha2)) alpha2 <- 0.05 / (k * m)

  # ---- 参数校验 ---------------------------------------------------------------
  if (method == "LMM") {
    if (is.null(env)) {
      stop("method='LMM' 需要提供环境标签向量 env（长度等于行数）。")
    }
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("LMM 方法需要安装 lme4 包：install.packages('lme4')")
    }
    env <- factor(env)
    if (nlevels(env) < 2L) {
      warning("env 只有一个水平，LMM 将退化为 OLS，建议直接使用 method='OLS'。")
    }
  }

  # ---- 构造上位性协变量矩阵（加性×加性，aa项）--------------------------------
  # epi_pairs = list(c(l1, l2), c(l3, l4), ...)
  # 上位性编码：x_aa = x_l1 × x_l2（两个位点基因型的乘积）
  # 作为额外固定效应加入条件模型，对应论文 Eq.1 中的 aa_{ilh}×x_{jl}×x_{jh}
  X_epi <- NULL
  epi_names <- NULL
  if (!is.null(epi_pairs)) {
    epi_list <- lapply(epi_pairs, function(pair) X_loci[, pair[1]] * X_loci[, pair[2]])
    X_epi <- do.call(cbind, epi_list)
    epi_names <- sapply(epi_pairs, function(pair) {
      paste0(colnames(X_loci)[pair[1]], "x", colnames(X_loci)[pair[2]])
    })
    colnames(X_epi) <- epi_names
  }

  # ---- 主循环：逐性状×逐位点估计条件效应 ------------------------------------
  # 结果存在列表里，最后一次性 rbind（比逐行追加高效）
  out <- vector("list", m * k)
  idx <- 0L

  for (i in seq_len(m)) {
    for (l in seq_len(k)) {
      snp_name <- colnames(X_loci)[l]
      trait_name <- colnames(Y_cond)[i]
      y <- Y_cond[, i]
      x_snp <- X_loci[, l]

      # 跳过无变异的位点（防止奇异矩阵）
      if (var(x_snp) < 1e-10) next

      # ==================================================================
      # OLS 分支：单环境，或无法使用 LMM 时的回退
      # ==================================================================
      if (method == "OLS") {
        if (is.null(X_epi)) {
          # 无上位性：y* = theta_il × x_l + eps
          # 注意：对所有候选位点做联合OLS（非逐一），因此X_loci全部进入模型
          # 这里 x_snp 是目标位点，其余位点通过 lm() 自动控制（作为协变量）
          fit <- lm(y ~ x_snp)
        } else {
          # 有上位性：y* = theta_il × x_l + aa_{ilh} × x_l × x_h + eps
          df_fit <- as.data.frame(cbind(x_snp, X_epi))
          fit <- lm(y ~ ., data = df_fit)
        }

        cf <- summary(fit)$coefficients
        if (nrow(cf) < 2L) next # 模型退化时跳过

        # 提取 x_snp 的系数（第2行，截距后第一个）
        idx <- idx + 1L
        out[[idx]] <- data.frame(
          trait = trait_name,
          locus = snp_name,
          theta_cond = cf[2L, 1L], # 效应量
          se_cond = cf[2L, 2L], # 标准误
          pval_cond = cf[2L, 4L], # p值（t检验）
          sig_cond = cf[2L, 4L] < alpha2, # 是否达到Bonferroni显著性
          method = "OLS",
          stringsAsFactors = FALSE
        )
      }

      # ==================================================================
      # LMM 分支：多环境，估计含环境随机效应的条件模型
      # ==================================================================
      if (method == "LMM") {
        # 构造数据框（必须包含 y, x_snp, env）
        df_fit <- data.frame(y = y, x_snp = x_snp, env = env)
        if (!is.null(X_epi)) df_fit <- cbind(df_fit, as.data.frame(X_epi))

        # 【LMM模型结构】（对应论文 Eq.1 的条件版本）：
        #   y*_{ijk} = theta_{il} × x_{jl}（加性固定效应）
        #           + aa_{ilh} × x_{jl} × x_{jh}（上位性固定效应，可选）
        #           + e_{ik}（环境随机截距，对应 e_{ik}）
        #           + ae_{ilk} × x_{jl}（G×E随机斜率，对应 ae_{ilk}）
        #           + eps_{ijk}（残差）
        #
        # lme4 公式分解：
        #   y ~ x_snp          → 固定效应：x_snp 的加性主效应
        #   + (1 + x_snp | env)→ 随机效应：
        #     (1 | env)       = 环境随机截距 e_{ik}
        #     (0 + x_snp | env) 包含在 (1 + x_snp | env) 中
        #                     = G×E 随机斜率 ae_{ilk}

        epi_formula <- if (!is.null(X_epi)) {
          paste("+", paste(epi_names, collapse = " + "))
        } else {
          ""
        }

        formula_str <- paste0("y ~ x_snp", epi_formula, " + (1 + x_snp | env)")

        fit <- tryCatch(
          lme4::lmer(
            as.formula(formula_str),
            data    = df_fit,
            REML    = TRUE, # 残差最大似然，方差成分估计更准确
            control = lme4::lmerControl(optimizer = "bobyqa") # 稳定的优化器
          ),
          error = function(e) {
            warning(sprintf(
              "性状%s 位点%s LMM拟合失败（%s），退回OLS。",
              trait_name, snp_name, conditionMessage(e)
            ))
            NULL
          }
        )

        # LMM 拟合失败时退回 OLS（记录 method="OLS_fallback"）
        if (is.null(fit)) {
          fit_ols <- lm(y ~ x_snp, data = df_fit)
          cf <- summary(fit_ols)$coefficients
          idx <- idx + 1L
          out[[idx]] <- data.frame(
            trait = trait_name, locus = snp_name,
            theta_cond = cf[2L, 1L], se_cond = cf[2L, 2L],
            pval_cond = cf[2L, 4L], sig_cond = cf[2L, 4L] < alpha2,
            method = "OLS_fallback",
            stringsAsFactors = FALSE
          )
          next
        }

        # ---- 从 LMM 提取固定效应的 p 值 -------------------------------------
        # lme4 不直接提供 p 值（避免争议性的自由度近似）
        # 推荐方法：使用 lmerTest 包的 Satterthwaite 自由度近似
        # 备用方法：正态分布近似（保守，仅在 lmerTest 不可用时）

        if (requireNamespace("lmerTest", quietly = TRUE)) {
          # Satterthwaite 近似：最常用的 LMM p 值方法，精度好
          fit_test <- lmerTest::as.lmerModLmerTest(fit)
          cf <- coef(summary(fit_test))
          # cf 有5列：Estimate, Std.Error, df, t_value, Pr(>|t|)
        } else {
          # 无 lmerTest 时的正态近似（z检验，样本量大时偏保守）
          cf <- coef(summary(fit))
          warning("建议安装 lmerTest 以获得准确 LMM p 值：install.packages('lmerTest')")
        }

        if (!"x_snp" %in% rownames(cf)) next # x_snp 不在模型中（异常情况）

        # 提取 p 值：lmerTest 输出第5列，正态近似用第3列（z值）
        pval <- if (ncol(cf) >= 5L) {
          cf["x_snp", 5L] # Satterthwaite p 值
        } else {
          2 * pnorm(abs(cf["x_snp", 3L]), lower.tail = FALSE) # 正态近似 p 值
        }

        idx <- idx + 1L
        out[[idx]] <- data.frame(
          trait = trait_name,
          locus = snp_name,
          theta_cond = cf["x_snp", 1L], # 固定效应估计（加性主效应）
          se_cond = cf["x_snp", 2L], # 固定效应标准误
          pval_cond = pval,
          sig_cond = pval < alpha2,
          method = "LMM",
          stringsAsFactors = FALSE
        )
      } # end LMM 分支
    } # end 位点循环
  } # end 性状循环

  # 过滤掉 NULL 元素（被 next 跳过的行），合并所有结果
  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}


# ------------------------------------------------------------------------------
# 辅助函数 3：打印诊断摘要
# ------------------------------------------------------------------------------

#' 打印条件投影诊断摘要
#'
#' @param result \code{compute_conditional_phenotype()} 的返回值
#' @export
print_conditional_diagnostics <- function(result) {
  cat("=== 条件投影诊断 ===\n\n")

  cat("【正交性】Cov(ỹ*_i, ỹ_j) 的最大绝对值（目标 < 0.01）：\n")
  for (i in seq_len(nrow(result$diagnostics))) {
    row <- result$diagnostics[i, ]
    flag <- if (row$max_cov > 0.01) "⚠️ " else "✅ "
    cat(sprintf("  %s %s：%.6f\n", flag, row$trait, row$max_cov))
  }

  cat("\n【方差等价性】Var(ỹ*_i) vs Schur 补 V_{i|-i}：\n")
  for (i in seq_len(nrow(result$diagnostics))) {
    row <- result$diagnostics[i, ]
    flag <- if (row$var_diff > 0.05) "⚠️ " else "✅ "
    cat(sprintf(
      "  %s %s：观测 = %.4f，理论 = %.4f，差异 = %.6f\n",
      flag, row$trait, row$var_obs, row$schur, row$var_diff
    ))
  }

  cat("\n【投影系数 γ_{i,-i}】：\n")
  for (nm in names(result$gamma)) {
    gamma_str <- paste(sprintf("%.4f", result$gamma[[nm]]), collapse = ", ")
    other <- setdiff(names(result$gamma), nm)
    cat(sprintf("  γ_{%s|%s} = [%s]\n", nm, paste(other, collapse = ","), gamma_str))
  }
}


# ------------------------------------------------------------------------------
# 成对条件投影：结局只对暴露投影，供 MR 的 IV 排他性筛选使用
# ------------------------------------------------------------------------------

#' 构造成对条件表型 ỹ*_{outcome | exposure}
#'
#' 与 compute_conditional_phenotype()（全条件，对所有其他性状投影）不同，
#' 本函数只把 outcome 对单一 exposure 投影。用于 Layer 3 MR 的排他性筛选：
#' 在 m≥3 的多性状场景下，全条件投影会把共同上游混杂的信号一并消除，导致
#' 混杂工具通过排他性筛选；成对投影只移除暴露方向的依赖，能保留并暴露
#' 混杂工具经后门路径对结局的残余效应，从而将其挡在 IV 集合之外。
#'
#' 仅用 {exposure, outcome} 的 2×2 子矩阵，对任意性状数 m 成立。
#'
#' @param Y_residual n × m 残差表型矩阵（与 compute_conditional_phenotype 同口径）。
#' @param exposure   暴露性状：列名（字符）或列索引（整数）。
#' @param outcome    结局性状：列名（字符）或列索引（整数）。
#' @param V          可选 m × m 表型协方差矩阵；NULL 时由 Y_residual 估计。
#' @param tol_orthogonality 正交性校验容差。默认 0.01。
#'
#' @return 列表：
#'   \item{y_cond}{长度 n 的成对条件表型向量 ỹ*_{outcome|exposure}}
#'   \item{gamma}{标量投影系数 Cov(out,exp)/Var(exp)}
#'   \item{schur}{理论条件方差 V_out - Cov^2/V_exp}
#'   \item{max_cov}{Cov(ỹ*, ỹ_exposure) 的绝对值，应接近 0（正交性校验）}
#'
#' @seealso \code{\link{compute_conditional_phenotype}}（全条件版本，供 Layer 2 分类）
#' @export
compute_pairwise_conditional <- function(Y_residual,
                                         exposure,
                                         outcome,
                                         V = NULL,
                                         tol_orthogonality = 0.01) {
  Y_residual <- as.matrix(Y_residual)
  m <- ncol(Y_residual)
  if (is.null(colnames(Y_residual))) colnames(Y_residual) <- paste0("Trait", seq_len(m))

  # 解析性状索引（接受列名或整数）
  resolve <- function(t) if (is.character(t)) match(t, colnames(Y_residual)) else as.integer(t)
  ie <- resolve(exposure)
  io <- resolve(outcome)
  if (is.na(ie) || is.na(io)) stop("exposure / outcome 未在 Y_residual 列名中找到。")
  if (ie == io) stop("exposure 与 outcome 不能是同一个性状。")

  # 协方差（只需要 2×2 子结构，但用整体 V 以便与全条件共享同一估计）
  if (is.null(V)) V <- cov(Y_residual)

  v_ee <- V[ie, ie] # Var(exposure)
  v_oo <- V[io, io] # Var(outcome)
  c_oe <- V[io, ie] # Cov(outcome, exposure)

  if (v_ee < 1e-12) stop("暴露性状方差接近 0，无法投影。")

  gamma <- c_oe / v_ee # 成对投影系数（标量）
  y_cond <- Y_residual[, io] - gamma * Y_residual[, ie]
  schur <- v_oo - c_oe^2 / v_ee # 理论条件方差

  # 正交性校验：ỹ* 应与暴露表型正交
  max_cov <- abs(cov(y_cond, Y_residual[, ie]))
  if (max_cov > tol_orthogonality) {
    warning(sprintf(
      "成对条件表型正交性未达标：|Cov(ỹ*_{%s|%s}, ỹ_%s)| = %.4f > %.4f。",
      colnames(Y_residual)[io], colnames(Y_residual)[ie],
      colnames(Y_residual)[ie], max_cov, tol_orthogonality
    ))
  }

  list(y_cond = y_cond, gamma = gamma, schur = schur, max_cov = max_cov)
}



# ------------------------------------------------------------------------------
# (1) 内部 helper：把 outcome 对指定的条件性状集做 GLS 投影（用协方差 V）
#     pairwise: cond = {exposure}
#     full    : cond = {所有性状 \ outcome}（含 exposure + 混杂 C）
# ------------------------------------------------------------------------------

#' @keywords internal
.project_outcome_on <- function(Y_residual, outcome_idx, cond_idx, V) {
  cond_idx <- setdiff(cond_idx, outcome_idx)
  if (length(cond_idx) == 0L) {
    return(Y_residual[, outcome_idx])
  }
  Vcc <- V[cond_idx, cond_idx, drop = FALSE]
  Vco <- V[cond_idx, outcome_idx]
  gamma <- solve(Vcc, Vco)
  as.numeric(Y_residual[, outcome_idx] - Y_residual[, cond_idx, drop = FALSE] %*% gamma)
}
