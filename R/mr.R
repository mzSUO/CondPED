# ==============================================================================
# mr.R — Layer 3：孟德尔随机化因果推断
# ==============================================================================
# 实现论文 Section 2.3.1c + 2.3.2 的三步流程：
#
#   Step 1  select_instruments()     工具变量筛选（四准则）
#   Step 2  estimate_causal_gls()    GLS 因果效应估计（Eq.7）+ Bootstrap SE
#   Step 3  bidirectional_mr()       双向 MR + 显著性判定
#   Step 4  test_mediation_consistency()  反向条件一致性验证（Section 2.3.2）
#
# ── 输入/输出接口 ──────────────────────────────────────────────────────────────
#
# select_instruments(marginal_effects, cond_effects, exposure, outcome, ...)
#   marginal_effects : qtxnetwork.output.trans()$qtl_data
#                      data.frame(TRAIT, QTL, A[char], SE[char], P_Value[char])
#   cond_effects     : fit_conditional_model()
#                      data.frame(trait, locus, theta_cond, se_cond, pval_cond, sig_cond)
#   exposure/outcome : 性状名称字符串
#   → character vector of valid IV locus IDs
#
# estimate_causal_gls(theta_exp, theta_out, ld_matrix, Y, X_iv, ...)
#   theta_exp/theta_out : named numeric vectors (names = locus IDs)
#   ld_matrix           : k×k numeric matrix of |r| values; NULL → identity
#   Y                   : n×m phenotype matrix (for bootstrap SE)
#   X_iv                : n×k genotype matrix of IVs (for bootstrap SE)
#   → list(gamma, se, z, pval, n_iv, method)
#
# bidirectional_mr(marginal_effects, cond_effects, traits, Y, X_all, ...)
#   traits  : character(2) c("TraitA", "TraitB")
#   X_all   : full n×p genotype matrix; columns subsetted to IVs internally
#   → list(AB = list(gamma, se, z, pval, sig, n_iv),
#           BA = list(...))
#
# test_mediation_consistency(theta_B_marg, theta_B_cond)
#   theta_B_marg : named numeric: IV marginal effects on outcome B
#   theta_B_cond : named numeric: IV conditional effects on B | A
#   → list(pval, consistent, delta_median)
# ==============================================================================


# ------------------------------------------------------------------------------
# 内部辅助：F 统计量
# ------------------------------------------------------------------------------

#' @keywords internal
compute_f_stat <- function(A, SE) {
  A_num  <- suppressWarnings(as.numeric(A))
  SE_num <- suppressWarnings(as.numeric(SE))
  (A_num / SE_num)^2
}


# ------------------------------------------------------------------------------
# 内部辅助：LD 修剪（贪心算法）
# ------------------------------------------------------------------------------

#' @keywords internal
ld_prune_ivs <- function(loci, ld_mat, r2_threshold = 0.1) {
  if (length(loci) <= 1L) return(loci)
  keep <- character(0L)
  for (l in loci) {
    if (length(keep) == 0L) {
      keep <- l
      next
    }
    max_r <- max(abs(ld_mat[l, keep, drop = FALSE]))
    if (max_r < r2_threshold) keep <- c(keep, l)
  }
  keep
}


# ------------------------------------------------------------------------------
# Step 1：工具变量筛选
# ------------------------------------------------------------------------------

#' 筛选 MR 有效工具变量（论文 Section 2.3.1a 四准则）
#'
#' @param marginal_effects 数据框，来自 \code{qtxnetwork.output.trans()$qtl_data}。
#'   必须含列：\code{TRAIT, QTL, A, SE, P_Value}。
#' @param cond_effects 数据框，来自 \code{fit_conditional_model()}。
#'   必须含列：\code{trait, locus, pval_cond}。
#' @param exposure 暴露性状名称（字符串）。
#' @param outcome  结果性状名称（字符串）。
#' @param alpha1 全基因组显著性阈值（由 QTLNetwork 置换检验确定）。
#' @param alpha2 排他性检验 p 值截断（默认 0.05）。
#'   位点对结果性状的条件效应 p > alpha2 → 与排他性约束相容。
#' @param F_threshold 弱工具变量过滤 F 统计量阈值（默认 10）。
#' @param r2_threshold LD 修剪阈值（|r|，默认 0.1）。
#' @param ld_matrix 位点间 LD 相关矩阵（|r| 值，非 r²）。
#'   \code{NULL} 时跳过 LD 修剪。
#'
#' @return 有效工具变量的位点 ID 字符向量。若不足 3 个，返回空向量并给出警告。
#'
#' @details
#' 四准则：
#' \enumerate{
#'   \item 相关性：位点对暴露性状边际显著（\code{P_Value < alpha1}）
#'   \item 排他性相容：位点对结果性状的条件效应不显著（\code{pval_cond > alpha2}）
#'   \item 弱工具变量过滤：F 统计量 = (A/SE)² > F_threshold
#'   \item LD 修剪：位点间 |r| < r2_threshold（当提供 ld_matrix 时）
#' }
#' 注意：准则 2 是"排他性兼容"（empirical support），不是正式证明排他性成立。
#'
#' @export
select_instruments <- function(marginal_effects,
                               cond_effects,
                               exposure,
                               outcome,
                               alpha1,
                               alpha2       = 0.05,
                               F_threshold  = 10,
                               r2_threshold = 0.1,
                               ld_matrix    = NULL) {

  # ── 准则 1：暴露性状边际显著 ────────────────────────────────────────────────
  exp_sig <- marginal_effects[
    marginal_effects$TRAIT == exposure &
    suppressWarnings(as.numeric(marginal_effects$P_Value)) < alpha1,
    "QTL"
  ]
  if (length(exp_sig) == 0L) {
    warning(sprintf("暴露性状 '%s' 无边际显著位点，无法构建工具变量集。", exposure))
    return(character(0L))
  }

  # ── 准则 2：结果性状条件效应不显著（排他性相容） ────────────────────────────
  excl_ok <- cond_effects[
    cond_effects$trait == outcome &
    cond_effects$locus %in% exp_sig &
    cond_effects$pval_cond > alpha2,
    "locus"
  ]
  candidates <- intersect(exp_sig, excl_ok)

  # ── 准则 3：F > F_threshold（弱工具变量过滤） ────────────────────────────────
  exp_df <- marginal_effects[
    marginal_effects$TRAIT == exposure &
    marginal_effects$QTL   %in% candidates,
  ]
  f_stats <- setNames(
    compute_f_stat(exp_df$A, exp_df$SE),
    exp_df$QTL
  )
  f_pass     <- names(f_stats)[!is.na(f_stats) & f_stats > F_threshold]
  candidates <- intersect(candidates, f_pass)

  # ── 准则 4：LD 修剪 ──────────────────────────────────────────────────────────
  if (!is.null(ld_matrix) && length(candidates) > 1L) {
    shared   <- intersect(candidates, rownames(ld_matrix))
    if (length(shared) > 1L) {
      sub_ld   <- ld_matrix[shared, shared, drop = FALSE]
      candidates <- ld_prune_ivs(shared, sub_ld, r2_threshold)
    }
  }

  if (length(candidates) < 3L) {
    warning(sprintf(
      "方向 %s→%s：有效 IV 数 = %d < 3，MR 结果不可靠，返回空集。",
      exposure, outcome, length(candidates)
    ))
    return(character(0L))
  }

  candidates
}


# ------------------------------------------------------------------------------
# Step 2：GLS 因果效应估计（论文 Eq.7）
# ------------------------------------------------------------------------------

#' GLS 因果效应估计 + Bootstrap 标准误
#'
#' 实现论文 Eq.7：
#' \deqn{\hat{\gamma}_{A \to B} = \frac{\hat{\boldsymbol{\theta}}_A^\top \mathbf{R}^{-1} \hat{\boldsymbol{\theta}}_B}{\hat{\boldsymbol{\theta}}_A^\top \mathbf{R}^{-1} \hat{\boldsymbol{\theta}}_A}}
#'
#' LD 修剪后 R ≈ I，GLS 退化为加权最小二乘（WLS）。
#' 标准误通过以个体为单元的 Bootstrap 重抽样估计（n_boot 次）。
#'
#' @param theta_exp 命名数值向量：工具变量对暴露性状的边际效应。
#'   名称为位点 ID，顺序需与 \code{theta_out} 一致。
#' @param theta_out 命名数值向量：工具变量对结果性状的边际效应。
#' @param ld_matrix k×k 矩阵（|r| 值）。\code{NULL} 时用单位矩阵（适用于独立 SNP）。
#' @param Y n×m 个体水平表型矩阵（含暴露和结果性状）。提供时使用 Bootstrap SE。
#' @param X_iv n×k 工具变量基因型矩阵（列与 theta_exp 顺序一致）。
#' @param exp_col 暴露性状在 Y 中的列索引（默认 1）。
#' @param out_col 结果性状在 Y 中的列索引（默认 2）。
#' @param n_boot Bootstrap 重抽样次数（默认 1000）。
#'
#' @return 列表，包含：
#'   \item{gamma}{因果效应点估计}
#'   \item{se}{标准误（Bootstrap 或解析近似）}
#'   \item{z}{Wald 检验统计量}
#'   \item{pval}{双侧 p 值}
#'   \item{n_iv}{使用的工具变量数}
#'   \item{method}{"bootstrap" 或 "delta"}
#'
#' @export
estimate_causal_gls <- function(theta_exp,
                                theta_out,
                                ld_matrix = NULL,
                                Y         = NULL,
                                X_iv      = NULL,
                                exp_col   = 1L,
                                out_col   = 2L,
                                n_boot    = 1000L) {

  k <- length(theta_exp)
  stopifnot(length(theta_out) == k, k >= 3L)

  # ── GLS 因果效应点估计 ───────────────────────────────────────────────────────
  if (is.null(ld_matrix)) {
    R_inv <- diag(k)
  } else {
    loci  <- names(theta_exp)
    sub   <- ld_matrix[loci, loci, drop = FALSE]
    R_inv <- tryCatch(solve(sub),
                      error = function(e) { warning("LD 矩阵奇异，退化为单位矩阵。"); diag(k) })
  }

  gamma_hat <- as.numeric(
    (t(theta_exp) %*% R_inv %*% theta_out) /
    (t(theta_exp) %*% R_inv %*% theta_exp)
  )

  # ── 标准误 ──────────────────────────────────────────────────────────────────
  method <- "delta"
  if (!is.null(Y) && !is.null(X_iv)) {
    method  <- "bootstrap"
    n       <- nrow(Y)
    gammas  <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      idx  <- sample.int(n, n, replace = TRUE)
      Yb   <- Y[idx, , drop = FALSE]
      Xb   <- X_iv[idx, , drop = FALSE]
      # 逐位点 OLS 回归重估效应
      tA_b <- apply(Xb, 2, function(x) coef(lm.fit(cbind(1, x), Yb[, exp_col]))[2])
      tB_b <- apply(Xb, 2, function(x) coef(lm.fit(cbind(1, x), Yb[, out_col]))[2])
      denom <- sum(tA_b^2)
      gammas[b] <- if (abs(denom) > 1e-12) sum(tA_b * tB_b) / denom else NA_real_
    }
    se_gamma <- sd(gammas, na.rm = TRUE)
  } else {
    # 解析近似（delta 方法，Fieller 1954 简化版）
    sigma2   <- sum((theta_out - gamma_hat * theta_exp)^2) / max(k - 1L, 1L)
    denom    <- as.numeric(t(theta_exp) %*% R_inv %*% theta_exp)
    se_gamma <- if (abs(denom) > 1e-12) sqrt(sigma2 / denom) else NA_real_
  }

  z    <- gamma_hat / se_gamma
  pval <- 2 * pnorm(-abs(z))

  list(gamma = gamma_hat, se = se_gamma, z = z, pval = pval,
       n_iv = k, method = method)
}


# ------------------------------------------------------------------------------
# Step 3：双向 MR
# ------------------------------------------------------------------------------

#' 双向孟德尔随机化
#'
#' 对性状对 (traits[1], traits[2]) 执行双向 MR：A→B 和 B→A。
#' 工具变量按方向独立筛选；每个方向的因果效应独立估计。
#'
#' @param marginal_effects 数据框，来自 \code{qtxnetwork.output.trans()$qtl_data}。
#' @param cond_effects 数据框，来自 \code{fit_conditional_model()}。
#' @param traits 长度为 2 的字符向量，性状名称 c("TraitA", "TraitB")。
#' @param Y n×m 表型矩阵（用于 Bootstrap SE）。\code{NULL} 时使用解析近似。
#' @param X_all n×p 全基因型矩阵（将按 IV 列名自动截取）。
#' @param ld_matrix LD 相关矩阵（可选）。
#' @param alpha1 全基因组显著性阈值。
#' @param alpha2 排他性检验截断（默认 0.05）。
#' @param F_threshold 弱工具变量过滤阈值（默认 10）。
#' @param r2_threshold LD 修剪阈值（默认 0.1）。
#' @param alpha_mr MR 显著性判定阈值（默认 0.05）。
#' @param n_boot Bootstrap 次数（默认 1000）。
#'
#' @return 列表，包含：
#'   \item{AB}{\code{list(gamma, se, z, pval, sig, n_iv)} 方向 A→B}
#'   \item{BA}{\code{list(...)}}                            方向 B→A
#'   \item{iv_AB}{A→B 使用的工具变量 ID 向量}
#'   \item{iv_BA}{B→A 使用的工具变量 ID 向量}
#'
#' @export
bidirectional_mr <- function(marginal_effects,
                             cond_effects,
                             traits,
                             Y            = NULL,
                             X_all        = NULL,
                             ld_matrix    = NULL,
                             alpha1,
                             alpha2       = 0.05,
                             F_threshold  = 10,
                             r2_threshold = 0.1,
                             alpha_mr     = 0.05,
                             n_boot       = 1000L) {

  stopifnot(length(traits) == 2L)
  trA <- traits[1L]; trB <- traits[2L]

  # 内部函数：运行单方向 MR
  run_direction <- function(exposure, outcome, exp_col, out_col) {
    ivs <- select_instruments(marginal_effects, cond_effects,
                              exposure, outcome, alpha1, alpha2,
                              F_threshold, r2_threshold, ld_matrix)

    if (length(ivs) == 0L) {
      return(list(gamma = NA_real_, se = NA_real_, z = NA_real_,
                  pval = NA_real_, sig = FALSE, n_iv = 0L, ivs = character(0L)))
    }

    # 提取 IV 的暴露和结果边际效应
    get_eff <- function(trait) {
      rows <- marginal_effects[marginal_effects$TRAIT == trait &
                               marginal_effects$QTL   %in% ivs, ]
      setNames(suppressWarnings(as.numeric(rows$A)), rows$QTL)[ivs]
    }

    theta_exp <- get_eff(exposure)
    theta_out <- get_eff(outcome)

    # 移除缺失
    ok        <- !is.na(theta_exp) & !is.na(theta_out)
    theta_exp <- theta_exp[ok]; theta_out <- theta_out[ok]
    ivs_used  <- names(theta_exp)

    if (length(ivs_used) < 3L) {
      warning(sprintf("方向 %s→%s：可用 IV = %d < 3，跳过 MR。", exposure, outcome, length(ivs_used)))
      return(list(gamma = NA_real_, se = NA_real_, z = NA_real_,
                  pval = NA_real_, sig = FALSE, n_iv = length(ivs_used), ivs = ivs_used))
    }

    # 准备 Bootstrap 所需的个体水平数据
    X_iv_boot <- if (!is.null(X_all) && !is.null(Y)) {
      cols <- intersect(ivs_used, colnames(X_all))
      if (length(cols) == length(ivs_used)) X_all[, cols, drop = FALSE] else NULL
    } else NULL

    Y_pair    <- if (!is.null(Y)) Y[, c(exp_col, out_col), drop = FALSE] else NULL

    gls <- estimate_causal_gls(theta_exp, theta_out, ld_matrix,
                               Y = Y_pair, X_iv = X_iv_boot,
                               exp_col = 1L, out_col = 2L, n_boot = n_boot)

    c(gls, list(sig = !is.na(gls$pval) && gls$pval < alpha_mr, ivs = ivs_used))
  }

  # 确定各性状在 Y 中的列号
  trait_col <- function(tr) {
    if (!is.null(Y) && !is.null(colnames(Y))) {
      idx <- which(colnames(Y) == tr)
      if (length(idx) == 1L) return(idx)
    }
    if (tr == trA) 1L else 2L
  }

  mr_AB <- run_direction(trA, trB, trait_col(trA), trait_col(trB))
  mr_BA <- run_direction(trB, trA, trait_col(trB), trait_col(trA))

  list(AB     = mr_AB,
       BA     = mr_BA,
       iv_AB  = mr_AB$ivs,
       iv_BA  = mr_BA$ivs)
}


# ------------------------------------------------------------------------------
# Step 4：介导一致性验证（论文 Section 2.3.2）
# ------------------------------------------------------------------------------

#' 反向条件一致性验证
#'
#' 验证 MR 支持的因果通路 A→B 的介导机制：
#' 控制上游性状 A 后，工具变量对下游性状 B 的效应是否系统性减弱？
#'
#' 使用单侧 Wilcoxon 符号秩检验，检验配对差值
#' δ_l = θ_{Bl}^{marg} − θ_{B|A,l}^{cond} 是否随机大于零。
#' P < 0.05 视为通过一致性检验。
#'
#' @param theta_B_marg 命名数值向量：工具变量对结果性状 B 的**边际**效应。
#' @param theta_B_cond 命名数值向量：工具变量对结果性状 B 的**条件**效应（控制 A 后）。
#'   名称顺序不必相同；函数内部按名称对齐。
#'
#' @return 列表，包含：
#'   \item{pval}{单侧 Wilcoxon 符号秩检验 p 值}
#'   \item{consistent}{逻辑值，P < 0.05 时为 TRUE}
#'   \item{delta_median}{δ 的中位数（正值表示效应系统减弱）}
#'   \item{delta}{各工具变量的效应差 δ_l 向量}
#'
#' @export
test_mediation_consistency <- function(theta_B_marg, theta_B_cond) {

  shared <- intersect(names(theta_B_marg), names(theta_B_cond))
  if (length(shared) < 3L) {
    warning("共享工具变量 < 3，一致性检验不可靠。")
    return(list(pval = NA_real_, consistent = FALSE,
                delta_median = NA_real_, delta = numeric(0L)))
  }

  delta <- theta_B_marg[shared] - theta_B_cond[shared]

  wt <- tryCatch(
    wilcox.test(delta, mu = 0, alternative = "greater"),
    error = function(e) list(p.value = NA_real_)
  )

  list(
    pval         = wt$p.value,
    consistent   = !is.na(wt$p.value) && wt$p.value < 0.05,
    delta_median = median(delta),
    delta        = delta
  )
}