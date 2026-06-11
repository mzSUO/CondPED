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

  V   <- cov(Y_residual)  # m × m 样本协方差矩阵
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

  idx       <- setdiff(seq_len(m), i)   # 其余性状的索引
  V_minus_i <- V[idx, idx, drop = FALSE] # (m-1) × (m-1)
  C_minus_i <- V[idx, i,   drop = FALSE] # (m-1) × 1

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
#'   # 生成示例数据（ρ = 0.5 的双性状场景）
#'   set.seed(42)
#'   n  <- 500
#'   Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#'   Y  <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
#'   colnames(Y) <- c("A", "B")
#'
#'   result <- compute_conditional_phenotype(Y)
#'
#'   # 验证正交性：应接近 0
#'   cov(result$Y_cond[, "A"], Y[, "B"])
#'
#'   # 验证方差等价：应接近 Schur 补
#'   var(result$Y_cond[, "A"])
#'   result$schur["A"]
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
                                          V                  = NULL,
                                          tol_orthogonality  = 0.01) {

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
      Y_cond      = Y_residual,
      gamma       = list(numeric(0)),
      V           = matrix(var(Y_residual), 1, 1),
      schur       = setNames(var(Y_residual), trait_names),
      diagnostics = data.frame(trait = trait_names, max_cov = NA,
                               var_obs = var(Y_residual), schur = var(Y_residual),
                               var_diff = 0)
    ))
  }

  # ---------- 估计协方差矩阵 --------------------------------------------------
  if (is.null(V)) {
    cov_obj <- estimate_phenotypic_covariance(Y_residual)
    V       <- cov_obj$V
  } else {
    stopifnot(
      "V 必须是 m×m 矩阵" = is.matrix(V) && nrow(V) == m && ncol(V) == m
    )
  }

  # ---------- 主循环：对每个性状构造条件表型 -----------------------------------
  Y_cond     <- matrix(NA_real_, n, m, dimnames = list(rownames(Y_residual), trait_names))
  gamma_list <- vector("list", m)
  schur      <- numeric(m)

  for (i in seq_len(m)) {
    idx <- setdiff(seq_len(m), i)

    # 投影系数：γ_{i,-i} = V_{-i}^{-1} C_{-i,i}
    gamma     <- compute_projection_coefficients(V, i)
    gamma_list[[i]] <- gamma

    # 条件表型：ỹ*_i = ỹ_i - ỹ_{-i} %*% γ
    Y_cond[, i] <- Y_residual[, i] - Y_residual[, idx, drop = FALSE] %*% gamma

    # 理论条件方差（Schur 补）：V_{i|-i} = V_{ii} - C_{i,-i} V_{-i}^{-1} C_{-i,i}
    # 等价于 V_{ii} - t(C_{-i,i}) %*% gamma
    C_minus_i  <- V[idx, i]
    schur[i]   <- V[i, i] - sum(C_minus_i * gamma)
  }

  names(gamma_list) <- trait_names
  names(schur)      <- trait_names

  # ---------- 诊断：验证数学性质 -----------------------------------------------
  diagnostics <- vector("list", m)

  for (i in seq_len(m)) {
    idx <- setdiff(seq_len(m), i)

    # 性质1：正交性 —— Cov(ỹ*_i, ỹ_j) 对所有 j ≠ i 应接近 0
    covs     <- apply(Y_residual[, idx, drop = FALSE], 2,
                      function(col) cov(Y_cond[, i], col))
    max_cov  <- max(abs(covs))

    if (max_cov > tol_orthogonality) {
      warning(sprintf(
        paste0("性状 '%s' 的条件表型正交性未达标：",
               "max|Cov(ỹ*_%s, ỹ_j)| = %.4f > %.4f。",
               "\n  可能原因：样本量不足或 V 估计误差过大。"),
        trait_names[i], trait_names[i], max_cov, tol_orthogonality
      ))
    }

    # 性质2：方差等价 —— Var(ỹ*_i) 应接近 Schur 补 V_{i|-i}
    var_obs  <- var(Y_cond[, i])
    var_diff <- abs(var_obs - schur[i])

    # 允许的偏差：受限于样本量，理论上偏差量级为 O(1/n)
    tol_var  <- 3 * sqrt(2 * schur[i]^2 / (n - 1))  # 基于 χ² 分布的近似容差
    if (var_diff > tol_var) {
      warning(sprintf(
        paste0("性状 '%s' 的方差等价性验证失败：",
               "Var(ỹ*) = %.4f，Schur 补 = %.4f，差异 = %.4f。"),
        trait_names[i], var_obs, schur[i], var_diff
      ))
    }

    diagnostics[[i]] <- data.frame(
      trait    = trait_names[i],
      max_cov  = max_cov,      # 正交性指标（越小越好）
      var_obs  = var_obs,      # 观测方差
      schur    = schur[i],     # 理论条件方差（Schur 补）
      var_diff = var_diff,     # 方差偏差（越小越好）
      stringsAsFactors = FALSE
    )
  }

  diagnostics_df <- do.call(rbind, diagnostics)

  # ---------- 返回 ------------------------------------------------------------
  list(
    Y_cond      = Y_cond,         # n×m 条件表型矩阵（核心输出）
    gamma       = gamma_list,     # 各性状的投影系数
    V           = V,              # 使用的协方差矩阵
    schur       = schur,          # 各性状的理论条件方差
    diagnostics = diagnostics_df  # 验证诊断结果
  )
}


# ------------------------------------------------------------------------------
# 主函数 2：对 Layer 1 显著位点估计条件遗传效应（对应论文 Eq.8）
# ------------------------------------------------------------------------------

#' 估计候选位点的条件遗传效应
#'
#' 对每个性状的条件表型拟合线性模型，提取条件效应 θ^cond 及其显著性。
#' 对应论文 Eq.8 的简化实现（在模拟中用 OLS 近似混合模型）。
#'
#' @param Y_cond n × m 条件表型矩阵（由 \code{compute_conditional_phenotype} 返回的 \code{$Y_cond}）。
#' @param X_loci n × k 候选位点基因型矩阵，**仅包含 Layer 1 显著位点**。
#'   列名应与 \code{qtl_data$QTL} 中的位点 ID 一致。
#' @param alpha2 Bonferroni 校正显著性阈值。
#'   若为 \code{NULL}（默认），自动计算为 \code{0.05 / (k × m)}。
#'
#' @return 数据框，每行对应一个（性状, 位点）对，包含：
#'   \item{trait}{性状名}
#'   \item{locus}{位点 ID}
#'   \item{theta_cond}{条件效应点估计}
#'   \item{se_cond}{标准误}
#'   \item{pval_cond}{p 值}
#'   \item{sig_cond}{是否达到 Bonferroni 校正显著性（逻辑值）}
#'
#' @details
#' 显著性判定标准：\code{pval_cond < alpha2}（Bonferroni 校正，控制 FWER）。
#' 条件效应不显著不意味着直接效应不存在，仅反映当前样本量下未检测到。
#'
#' @seealso \code{\link{compute_conditional_phenotype}}
#' @export
fit_conditional_model <- function(Y_cond, X_loci, alpha2 = NULL) {

  Y_cond  <- as.matrix(Y_cond)
  X_loci  <- as.matrix(X_loci)
  m       <- ncol(Y_cond)
  k       <- ncol(X_loci)

  if (is.null(colnames(Y_cond)))  colnames(Y_cond)  <- paste0("Trait", seq_len(m))
  if (is.null(colnames(X_loci)))  colnames(X_loci)  <- paste0("SNP",   seq_len(k))

  if (is.null(alpha2)) alpha2 <- 0.05 / (k * m)

  out <- vector("list", m * k)
  idx <- 0L

  for (i in seq_len(m)) {
    for (l in seq_len(k)) {
      fit <- lm(Y_cond[, i] ~ X_loci[, l])
      cf  <- summary(fit)$coefficients

      if (nrow(cf) < 2L) {
        # 位点无变异（单态），跳过
        next
      }

      idx <- idx + 1L
      out[[idx]] <- data.frame(
        trait      = colnames(Y_cond)[i],
        locus      = colnames(X_loci)[l],
        theta_cond = cf[2L, 1L],
        se_cond    = cf[2L, 2L],
        pval_cond  = cf[2L, 4L],
        sig_cond   = cf[2L, 4L] < alpha2,
        stringsAsFactors = FALSE
      )
    }
  }

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
    row  <- result$diagnostics[i, ]
    flag <- if (row$max_cov > 0.01) "⚠️ " else "✅ "
    cat(sprintf("  %s %s：%.6f\n", flag, row$trait, row$max_cov))
  }

  cat("\n【方差等价性】Var(ỹ*_i) vs Schur 补 V_{i|-i}：\n")
  for (i in seq_len(nrow(result$diagnostics))) {
    row  <- result$diagnostics[i, ]
    flag <- if (row$var_diff > 0.05) "⚠️ " else "✅ "
    cat(sprintf("  %s %s：观测 = %.4f，理论 = %.4f，差异 = %.6f\n",
                flag, row$trait, row$var_obs, row$schur, row$var_diff))
  }

  cat("\n【投影系数 γ_{i,-i}】：\n")
  for (nm in names(result$gamma)) {
    gamma_str <- paste(sprintf("%.4f", result$gamma[[nm]]), collapse = ", ")
    other     <- setdiff(names(result$gamma), nm)
    cat(sprintf("  γ_{%s|%s} = [%s]\n", nm, paste(other, collapse = ","), gamma_str))
  }
}