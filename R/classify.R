# ==============================================================================
# classify.R — 六类分类决策（对应论文 Table 1）
# ==============================================================================
# 实现论文 Section 2.3.3 的统一分类体系：
#
#   classify_locus()       单个位点分类（核心函数）
#   classify_all_loci()    批量分类（对 condped 输出的包装）
#
# ── 分类判定逻辑 ──────────────────────────────────────────────────────────────
#
#  n_marg_sig == 0                         → null
#  n_marg_sig == 1                         → class1（性状特异性）
#  n_marg_sig >= 2：
#    无 MR 结果（仅 Layer 1+2 消融场景）：
#      n_cond_sig >= 2                     → class3（水平多效）
#      n_cond_sig <  2                     → class2（无法与 class4 区分）
#    双向 MR 均显著                         → class6（双向或混杂）
#    双向 MR 均不显著：
#      n_cond_sig >= 2                     → class3
#      n_cond_sig <  2                     → class2
#    单向 MR 显著（dir：exposure→outcome）：
#      outcome 性状条件效应显著             → class5（垂直部分中介）
#      outcome 性状条件效应不显著           → class4（垂直完全中介）
#
# ── 输入/输出接口 ──────────────────────────────────────────────────────────────
#
# classify_locus(n_marg_sig, cond_sig_by_trait, mr_results, traits)
#   n_marg_sig       : integer，边际显著性状数（来自 layer1$snp_count）
#   cond_sig_by_trait: named logical vector，各性状条件效应是否显著
#                      例：c(TraitA=TRUE, TraitB=FALSE)
#                      名称需与 traits 参数一致
#   mr_results       : bidirectional_mr() 的返回值，或 NULL（不做 MR）
#   traits           : character(2)，c("TraitA","TraitB")（决定方向映射）
#   → character: "null"/"class1"/"class2"/"class3"/"class4"/"class5"/"class6"
#
# classify_all_loci(condped_result)
#   condped_result   : condped() 的返回数据框（需含 class 列）
#   → 带标签的汇总表
# ==============================================================================


#' 单个位点的六类分类判定
#'
#' 依据论文 Table 1，综合边际显著性模式、条件效应显著性、双向 MR 结果，
#' 对一个 SNP-性状集合进行最终分类。
#'
#' @param n_marg_sig 整数。该位点在几个性状上具有边际显著效应。
#'   来源：\code{qtxnetwork.layer1.screen()$snp_count$n_sig}。
#' @param cond_sig_by_trait 命名逻辑向量。各性状条件效应是否显著。
#'   名称为性状名称，顺序需与 \code{traits} 参数一致。
#'   来源：\code{fit_conditional_model()} 中 \code{sig_cond} 列，按位点聚合。
#'   例：\code{c(TraitA = TRUE, TraitB = FALSE)}。
#' @param mr_results \code{bidirectional_mr()} 的返回值（列表）。
#'   若为 \code{NULL}，跳过 Layer 3，仅基于条件效应分类（消融实验场景）。
#' @param traits 长度为 2 的字符向量，性状名称 \code{c("TraitA", "TraitB")}。
#'   决定"AB 方向"和"BA 方向"对应哪个性状为结果性状。
#'
#' @return 字符串，取值为：
#'   \code{"null"}, \code{"class1"}, \code{"class2"}, \code{"class3"},
#'   \code{"class4"}, \code{"class5"}, \code{"class6"}。
#'
#' @details
#' \strong{Class 2 vs Class 4 区分依赖 Layer 3 MR}：
#' 两者在 Layer 2 具有相同指纹（B 条件效应不显著），但 Class 4 有单向
#' 显著 MR，Class 2 双向均不显著。当 \code{mr_results = NULL} 时，
#' 两者均归为 \code{"class2"}，这正是消融实验展示 Layer 3 必要性的证据。
#'
#' @examples
#' # Class 4 场景（垂直完全中介，A→B）
#' classify_locus(
#'   n_marg_sig       = 2L,
#'   cond_sig_by_trait = c(TraitA = TRUE, TraitB = FALSE),
#'   mr_results        = list(AB = list(sig = TRUE),  BA = list(sig = FALSE)),
#'   traits            = c("TraitA", "TraitB")
#' )
#' # → "class4"
#'
#' # Class 5 场景（垂直部分中介）
#' classify_locus(
#'   n_marg_sig        = 2L,
#'   cond_sig_by_trait = c(TraitA = TRUE, TraitB = TRUE),
#'   mr_results        = list(AB = list(sig = TRUE), BA = list(sig = FALSE)),
#'   traits            = c("TraitA", "TraitB")
#' )
#' # → "class5"
#'
#' @export
classify_locus <- function(n_marg_sig,
                           cond_sig_by_trait,
                           mr_results = NULL,
                           traits     = NULL) {

  # ── 前置校验 ──────────────────────────────────────────────────────────────
  n_marg_sig <- as.integer(n_marg_sig)

  # ── 零假设：边际不显著 ────────────────────────────────────────────────────
  if (is.na(n_marg_sig) || n_marg_sig == 0L) return("null")

  # ── Class 1：仅单一性状显著 ─────────────────────────────────────────────
  if (n_marg_sig == 1L) return("class1")

  # ── 多效性候选（≥2 个性状边际显著）──────────────────────────────────────
  n_cond_sig <- sum(cond_sig_by_trait, na.rm = TRUE)

  # ── 消融实验：无 MR，仅依赖 Layer 2 ─────────────────────────────────────
  if (is.null(mr_results)) {
    # 无 MR 时无法区分 Class 2 和 Class 4（两者 Layer 2 指纹相同）
    return(if (n_cond_sig >= 2L) "class3" else "class2")
  }

  # ── 提取 MR 显著性 ─────────────────────────────────────────────────────
  mr_AB_sig <- isTRUE(mr_results$AB$sig)
  mr_BA_sig <- isTRUE(mr_results$BA$sig)

  # ── Class 6：双向均显著（反馈回路或混杂）────────────────────────────────
  if (mr_AB_sig && mr_BA_sig) return("class6")

  # ── 双向均不显著 ──────────────────────────────────────────────────────
  if (!mr_AB_sig && !mr_BA_sig) {
    # 多个条件效应显著 → 水平多效；否则 → 协方差介导
    return(if (n_cond_sig >= 2L) "class3" else "class2")
  }

  # ── 单向 MR 显著：区分完全/部分中介（Class 4 vs Class 5）──────────────
  # 找出显著方向的"结果性状"，检验其条件效应
  #
  # 约定：traits[1]=暴露("A"), traits[2]=结果("B")
  #   MR_AB 显著 → 结果是 traits[2]（B 的条件效应）
  #   MR_BA 显著 → 结果是 traits[1]（A 的条件效应）
  if (mr_AB_sig) {
    outcome_trait <- if (!is.null(traits)) traits[2L] else names(cond_sig_by_trait)[2L]
  } else {
    outcome_trait <- if (!is.null(traits)) traits[1L] else names(cond_sig_by_trait)[1L]
  }

  outcome_cond_sig <- if (!is.null(outcome_trait) && outcome_trait %in% names(cond_sig_by_trait)) {
    isTRUE(cond_sig_by_trait[outcome_trait])
  } else {
    # fallback：任意一个条件效应显著
    n_cond_sig >= 1L
  }

  # 结果性状条件效应显著 → 部分中介；否则 → 完全中介
  if (outcome_cond_sig) "class5" else "class4"
}


#' 汇总 condped 分类结果
#'
#' 对 \code{condped()} 返回的数据框进行汇总，统计各类别的频数和比例。
#'
#' @param result \code{condped()} 的返回数据框，需含 \code{class} 列。
#' @param verbose 逻辑值，是否打印汇总信息（默认 \code{TRUE}）。
#'
#' @return 数据框：\code{class}, \code{n}（频数）, \code{pct}（百分比）。
#'
#' @export
summarise_classification <- function(result, verbose = TRUE) {

  if (!is.data.frame(result) || !"class" %in% colnames(result)) {
    stop("result 必须是包含 'class' 列的数据框（condped() 的输出）。")
  }

  class_levels <- c("null", "class1", "class2", "class3",
                     "class4", "class5", "class6")
  class_labels <- c(
    null   = "Null (undetected)",
    class1 = "Class 1: Trait-specific",
    class2 = "Class 2: Covariance-induced",
    class3 = "Class 3: Horizontal pleiotropy",
    class4 = "Class 4: Vertical – complete mediation",
    class5 = "Class 5: Vertical – partial mediation",
    class6 = "Class 6: Bidirectional / confounded"
  )

  tbl <- table(factor(result$class, levels = class_levels))
  df  <- data.frame(
    class = names(tbl),
    label = class_labels[names(tbl)],
    n     = as.integer(tbl),
    pct   = round(100 * as.integer(tbl) / nrow(result), 1),
    stringsAsFactors = FALSE
  )
  df <- df[df$n > 0, , drop = FALSE]

  if (verbose) {
    cat("\n=== CondPED Classification Summary ===\n")
    cat(sprintf("Total loci: %d\n\n", nrow(result)))
    for (i in seq_len(nrow(df))) {
      cat(sprintf("  %-45s  n = %3d  (%5.1f%%)\n",
                  df$label[i], df$n[i], df$pct[i]))
    }
    cat("\n")
  }

  invisible(df)
}