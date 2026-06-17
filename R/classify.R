# ==============================================================================
# classify.R — 五类分类决策（删除原 Class 2 后顺延编号）
# ==============================================================================
# 编号对照：
#   新 Class 1  = 原 Class 1  (Trait-specific)
#   新 Class 2  = 原 Class 3  (Horizontal pleiotropy)
#   新 Class 3  = 原 Class 4  (Vertical pleiotropy – complete mediation)
#   新 Class 4  = 原 Class 5  (Vertical pleiotropy – partial mediation)
#   新 Class 5  = 原 Class 6  (Bidirectional / confounded)
# ==============================================================================

#' 单个位点的五类分类判定
#'
#' 依据论文 Table 1（删除原 Class 2 协方差诱导后顺延编号），综合边际显著性、
#' 条件投影显著性与双向 MR 结果，对单个 SNP 进行最终分类。
#'
#' @param n_marg_sig 整数。该位点在 Layer 1 边际分析中显著关联的性状数。
#'   来源：\code{qtxnetwork.layer1.screen()$snp_count$n_sig}。
#' @param cond_sig_by_trait 命名逻辑向量。Layer 2 条件投影后，各性状是否仍显著。
#'   名称必须与 \code{traits} 对应。例：\code{c(TraitA = TRUE, TraitB = FALSE)}。
#' @param mr_results \code{bidirectional_mr()} 返回值（列表），或 \code{NULL}。
#'   若为 \code{NULL}，表示消融实验（无 Layer 3）。
#' @param traits 长度为 2 的字符向量，如 \code{c("TraitA", "TraitB")}。
#'   用于单向 MR 时判定 outcome 性状。消融实验可省略。
#'
#' @return 字符串：\code{"null"}, \code{"class1"}, \code{"class2"},
#'   \code{"class3"}, \code{"class4"}, \code{"class5"}。
#'
#' @details
#' \strong{消融实验降级逻辑}：当 \code{mr_results = NULL} 时，
#' 若 \code{n_cond_sig < 2}（投影后仅剩一个性状显著），框架无法区分
#' 真实 Class 1（单性状）与 Class 3（完全中介，原 Class 4），
#' 因此保守归为 \code{"class1"}，以此证明 Layer 3 MR 的必要性。
#'
#' @export
classify_locus <- function(n_marg_sig,
                           cond_sig_by_trait,
                           mr_results = NULL,
                           traits     = NULL) {

  n_marg_sig <- as.integer(n_marg_sig)

  # ── 0. 防御：输入校验 ─────────────────────────────────────────────────────
  if (is.null(traits) && is.null(names(cond_sig_by_trait))) {
    stop("traits 与 cond_sig_by_trait 名称不能同时为空。")
  }
  if (!is.null(traits) && length(traits) != 2L) {
    stop("traits 必须是长度为 2 的字符向量。")
  }

  # ── 1. Null：边际不显著 ───────────────────────────────────────────────────
  if (is.na(n_marg_sig) || n_marg_sig == 0L) {
    return("null")
  }

  # ── 2. Class 1：仅单一性状边际显著，不进入 Layer 2/3 ─────────────────────
  if (n_marg_sig == 1L) {
    return("class1")
  }

  # ── 3. 多效性候选（≥2 个性状边际显著）────────────────────────────────────
  n_cond_sig <- sum(cond_sig_by_trait, na.rm = TRUE)

  # ── 4. 消融实验：无 MR，仅 Layer 1+2 ──────────────────────────────────────
  if (is.null(mr_results)) {
    # n_cond >= 2 → 水平多效（Class 2，原 Class 3）
    # n_cond <  2 → 保守归为 Class 1（无法区分完全中介 vs 单性状）
    return(if (n_cond_sig >= 2L) "class2" else "class1")
  }

  # ── 5. 提取双向 MR 显著性 ──────────────────────────────────────────────────
  mr_AB_sig <- isTRUE(mr_results$AB$sig)
  mr_BA_sig <- isTRUE(mr_results$BA$sig)

  # ── 6. Class 5：双向均显著（反馈/混杂，原 Class 6）────────────────────────
  if (mr_AB_sig && mr_BA_sig) {
    return("class5")
  }

  # ── 7. 双向均不显著 ────────────────────────────────────────────────────────
  #   n_cond >= 2 → Class 2（水平多效，原 Class 3）
  #   n_cond <  2 → Class 1（保守降级；实践中协方差诱导已被 QTLNetwork 过滤）
  if (!mr_AB_sig && !mr_BA_sig) {
    return(if (n_cond_sig >= 2L) "class2" else "class1")
  }

  # ── 8. 单向 MR 显著：区分完全中介 vs 部分中介（Class 3 vs Class 4）──────
  #   防御：若投影后两个性状均不显著，但 MR 显著，数据矛盾，保守降级
  if (n_cond_sig == 0L) {
    return("class1")
  }

  # 判定 outcome 性状：
  #   AB 显著 → A 为暴露，B 为结局 → 检查 B 的条件效应
  #   BA 显著 → B 为暴露，A 为结局 → 检查 A 的条件效应
  if (mr_AB_sig) {
    outcome_trait <- if (!is.null(traits)) traits[2L] else names(cond_sig_by_trait)[2L]
  } else {
    outcome_trait <- if (!is.null(traits)) traits[1L] else names(cond_sig_by_trait)[1L]
  }

  outcome_cond_sig <- isTRUE(cond_sig_by_trait[outcome_trait])

  # outcome 条件效应仍显著 → 存在直接效应 + 间接效应 → 部分中介（Class 4，原 Class 5）
  # outcome 条件效应不显著 → 仅间接效应 → 完全中介（Class 3，原 Class 4）
  if (outcome_cond_sig) "class4" else "class3"
}