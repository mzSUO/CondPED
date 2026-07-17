# ==============================================================================
# classify.R — Pair-level association-pattern annotation for CondPED
# ==============================================================================
# Scientific scope
# ----------------
# CondPED annotates empirically distinguishable locus-trait-pair association
# patterns. The output is not a proof of a biological mechanism.
#
# Output patterns
# ---------------
#   not_detected : neither trait reaches marginal significance for the locus
#   pattern1     : only one trait reaches marginal significance
#   pattern2     : both traits are marginally significant, both conditional
#                  associations are retained, and neither MR direction is
#                  statistically supported
#   pattern3     : one MR direction is supported and the outcome conditional
#                  association is not detected (complete-attenuation pattern)
#   pattern4     : one MR direction is supported and the outcome conditional
#                  association is retained (partial-retention pattern)
#   pattern5     : both MR directions are supported (bidirectional/ambiguous)
#   unresolved   : insufficient IVs, failed estimation, conflicting diagnostics,
#                  or a conditional pattern that cannot be assigned reliably
#
# Main functions
# --------------
#   classify_pair_pattern()         one locus x one unordered trait pair
#   classify_multitrait_patterns()  all loci x all trait pairs
# ============================================================================== 


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

.normalize_marginal_results <- function(marginal_results) {
  aliases <- list(
    TRAIT = c("TRAIT", "trait"),
    SNPID = c("SNPID", "snp_id", "locus"),
    P_Value = c("P_Value", "p_value", "pval", "p")
  )

  find_name <- function(options) {
    hit <- intersect(options, names(marginal_results))
    if (length(hit) == 0L) NA_character_ else hit[1L]
  }

  selected <- vapply(aliases, find_name, character(1L))
  if (anyNA(selected)) {
    stop(
      "marginal_results 需要性状、SNP 和 P 值列。支持列名: ",
      "TRAIT/trait, SNPID/snp_id/locus, P_Value/p_value/pval/p。",
      call. = FALSE
    )
  }

  out <- data.frame(
    TRAIT = as.character(marginal_results[[selected["TRAIT"]]]),
    SNPID = as.character(marginal_results[[selected["SNPID"]]]),
    P_Value = suppressWarnings(as.numeric(marginal_results[[selected["P_Value"]]])),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(out[c("TRAIT", "SNPID")])) {
    stop("marginal_results 存在重复 TRAIT-SNPID 记录。", call. = FALSE)
  }
  out
}

.normalize_conditional_results <- function(conditional_results) {
  aliases <- list(
    trait = c("trait", "TRAIT"),
    locus = c("locus", "SNPID", "snp_id"),
    pval_cond = c("pval_cond", "p_value", "P_Value", "pval", "p"),
    sig_cond = c("sig_cond", "significant", "sig")
  )

  find_name <- function(options, required = TRUE) {
    hit <- intersect(options, names(conditional_results))
    if (length(hit) == 0L) {
      if (required) NA_character_ else ""
    } else {
      hit[1L]
    }
  }

  selected <- c(
    trait = find_name(aliases$trait),
    locus = find_name(aliases$locus),
    pval_cond = find_name(aliases$pval_cond, required = FALSE),
    sig_cond = find_name(aliases$sig_cond, required = FALSE)
  )

  if (anyNA(selected[c("trait", "locus")])) {
    stop("conditional_results 需要 trait 和 locus/SNPID/snp_id 列。", call. = FALSE)
  }
  if (selected["pval_cond"] == "" && selected["sig_cond"] == "") {
    stop("conditional_results 至少需要条件 P 值列或显著性逻辑列。", call. = FALSE)
  }

  out <- data.frame(
    trait = as.character(conditional_results[[selected["trait"]]]),
    locus = as.character(conditional_results[[selected["locus"]]]),
    stringsAsFactors = FALSE
  )
  out$pval_cond <- if (selected["pval_cond"] == "") {
    rep(NA_real_, nrow(out))
  } else {
    suppressWarnings(as.numeric(conditional_results[[selected["pval_cond"]]]))
  }
  out$sig_cond_input <- if (selected["sig_cond"] == "") {
    rep(NA, nrow(out))
  } else {
    as.logical(conditional_results[[selected["sig_cond"]]])
  }

  if (anyDuplicated(out[c("trait", "locus")])) {
    stop("conditional_results 存在重复 trait-locus 记录。", call. = FALSE)
  }
  out
}

.mr_direction_valid <- function(mr_result, min_iv = 3L) {
  if (!is.list(mr_result) || !identical(mr_result$status, "ok")) return(FALSE)
  required <- c("pval", "gamma", "n_iv")
  if (!all(required %in% names(mr_result))) return(FALSE)
  if (length(mr_result$pval) != 1L || length(mr_result$gamma) != 1L ||
      length(mr_result$n_iv) != 1L) return(FALSE)
  is.finite(mr_result$pval) &&
    is.finite(mr_result$gamma) &&
    is.finite(mr_result$n_iv) &&
    mr_result$n_iv >= min_iv
}

.mr_direction_significant <- function(mr_result, min_iv = 3L) {
  .mr_direction_valid(mr_result, min_iv) && isTRUE(mr_result$sig)
}

.direction_diagnostic_status <- function(mr_result,
                                         require_steiger = TRUE,
                                         reject_heterogeneity = TRUE,
                                         require_robust_direction = FALSE) {
  reasons <- character(0L)

  if (isTRUE(require_steiger)) {
    steiger_ok <- is.list(mr_result$steiger) && isTRUE(mr_result$steiger$consistent)
    if (!steiger_ok) reasons <- c(reasons, "steiger_inconsistent_or_unavailable")
  }

  if (isTRUE(reject_heterogeneity)) {
    heterogeneous <- is.list(mr_result$q) && isTRUE(mr_result$q$heterogeneous)
    if (heterogeneous) reasons <- c(reasons, "heterogeneous_instruments")
  }

  if (isTRUE(require_robust_direction)) {
    if (!isTRUE(mr_result$robust_direction_consistent)) {
      reasons <- c(reasons, "gls_weighted_median_conflict")
    }
  }

  list(ok = length(reasons) == 0L, reasons = reasons)
}


# ------------------------------------------------------------------------------
# One locus x one trait pair
# ------------------------------------------------------------------------------

#' Classify one locus-trait pair into an empirical association pattern
#'
#' @param marg_sig_trait1,marg_sig_trait2 Marginal significance indicators.
#' @param cond_sig_trait1,cond_sig_trait2 Conditional significance indicators.
#' @param mr_results Output of bidirectional_mr(). AB means trait1 -> trait2.
#' @param traits Character vector c(trait1, trait2).
#' @param min_iv Minimum IV count for a direction to be estimable.
#' @param require_steiger Require Steiger consistency before assigning pattern3/4.
#' @param reject_heterogeneity Return unresolved when the supported direction has
#'   significant Cochran-Q heterogeneity.
#' @param require_robust_direction Require GLS and weighted-median sign/scale
#'   consistency before assigning pattern3/4.
#' @param return_details Return a one-row data.frame rather than only the pattern.
#'
#' @export
classify_pair_pattern <- function(marg_sig_trait1,
                                  marg_sig_trait2,
                                  cond_sig_trait1,
                                  cond_sig_trait2,
                                  mr_results = NULL,
                                  traits = c("Trait1", "Trait2"),
                                  min_iv = 3L,
                                  require_steiger = TRUE,
                                  reject_heterogeneity = TRUE,
                                  require_robust_direction = FALSE,
                                  return_details = TRUE) {
  if (length(traits) != 2L || anyDuplicated(traits)) {
    stop("traits 必须包含两个不同性状名称。", call. = FALSE)
  }

  m1 <- isTRUE(marg_sig_trait1)
  m2 <- isTRUE(marg_sig_trait2)
  c1 <- isTRUE(cond_sig_trait1)
  c2 <- isTRUE(cond_sig_trait2)

  pattern <- NA_character_
  direction <- "none"
  confidence <- "not_assessed"
  unresolved_reason <- NA_character_
  forward_status <- if (is.null(mr_results) || is.null(mr_results$AB$status)) {
    NA_character_
  } else {
    as.character(mr_results$AB$status)[1L]
  }
  reverse_status <- if (is.null(mr_results) || is.null(mr_results$BA$status)) {
    NA_character_
  } else {
    as.character(mr_results$BA$status)[1L]
  }

  if (!m1 && !m2) {
    pattern <- "not_detected"
    confidence <- "not_applicable"
  } else if (xor(m1, m2)) {
    pattern <- "pattern1"
    direction <- "none"
    confidence <- "descriptive"
  } else if (is.null(mr_results)) {
    pattern <- "unresolved"
    unresolved_reason <- "mr_not_provided"
    confidence <- "unresolved"
  } else {
    forward_valid <- .mr_direction_valid(mr_results$AB, min_iv)
    reverse_valid <- .mr_direction_valid(mr_results$BA, min_iv)
    forward_sig <- .mr_direction_significant(mr_results$AB, min_iv)
    reverse_sig <- .mr_direction_significant(mr_results$BA, min_iv)

    if (!forward_valid || !reverse_valid) {
      pattern <- "unresolved"
      confidence <- "unresolved"
      statuses <- c(
        if (!forward_valid) paste0(traits[1L], "_to_", traits[2L], ":", mr_results$AB$status),
        if (!reverse_valid) paste0(traits[2L], "_to_", traits[1L], ":", mr_results$BA$status)
      )
      unresolved_reason <- paste(statuses, collapse = ";")
    } else if (forward_sig && reverse_sig) {
      pattern <- "pattern5"
      direction <- "bidirectional_or_ambiguous"
      confidence <- "ambiguous"
    } else if (!forward_sig && !reverse_sig) {
      if (c1 && c2) {
        pattern <- "pattern2"
        direction <- "none"
        confidence <- "moderate"
      } else {
        pattern <- "unresolved"
        confidence <- "unresolved"
        unresolved_reason <- "no_directional_support_with_attenuated_conditional_pattern"
      }
    } else {
      supported <- if (forward_sig) mr_results$AB else mr_results$BA
      supported_direction <- if (forward_sig) {
        paste0(traits[1L], "->", traits[2L])
      } else {
        paste0(traits[2L], "->", traits[1L])
      }
      outcome_cond_sig <- if (forward_sig) c2 else c1
      diag_status <- .direction_diagnostic_status(
        supported,
        require_steiger = require_steiger,
        reject_heterogeneity = reject_heterogeneity,
        require_robust_direction = require_robust_direction
      )

      if (!diag_status$ok) {
        pattern <- "unresolved"
        direction <- supported_direction
        confidence <- "unresolved"
        unresolved_reason <- paste(diag_status$reasons, collapse = ";")
      } else {
        pattern <- if (outcome_cond_sig) "pattern4" else "pattern3"
        direction <- supported_direction
        confidence <- "direction_supported"
      }
    }
  }

  result <- data.frame(
    trait1 = traits[1L],
    trait2 = traits[2L],
    marginal_sig_trait1 = m1,
    marginal_sig_trait2 = m2,
    conditional_sig_trait1 = c1,
    conditional_sig_trait2 = c2,
    pattern = pattern,
    direction = direction,
    confidence = confidence,
    unresolved_reason = unresolved_reason,
    mr_forward_status = forward_status,
    mr_reverse_status = reverse_status,
    stringsAsFactors = FALSE
  )

  if (isTRUE(return_details)) result else pattern
}


# ------------------------------------------------------------------------------
# All loci x all trait pairs
# ------------------------------------------------------------------------------

#' Annotate all locus-trait pairs in a multi-trait dataset
#'
#' @param marginal_results Layer-1 result table. Supported columns include
#'   TRAIT/SNPID/P_Value or trait/snp_id/p_value.
#' @param conditional_results Layer-2 result table. Supported columns include
#'   trait/locus/pval_cond/sig_cond.
#' @param mr_by_pair Named list returned by run_all_trait_pairs_mr().
#' @param traits Trait names and ordering.
#' @param alpha1 Layer-1 significance threshold.
#' @param alpha2 Layer-2 significance threshold, used when sig_cond is absent.
#' @param loci Optional SNP subset. Default is the union of marginal-result SNPs.
#' @param include_not_detected Include pair rows where neither trait is marginally
#'   significant for the locus.
#'
#' @return Long-format data.frame with one row per locus x unordered trait pair.
#' @export
classify_multitrait_patterns <- function(marginal_results,
                                         conditional_results,
                                         mr_by_pair,
                                         traits,
                                         alpha1,
                                         alpha2,
                                         loci = NULL,
                                         include_not_detected = FALSE,
                                         min_iv = 3L,
                                         require_steiger = TRUE,
                                         reject_heterogeneity = TRUE,
                                         require_robust_direction = FALSE) {
  marg <- .normalize_marginal_results(marginal_results)
  cond <- .normalize_conditional_results(conditional_results)
  traits <- unique(as.character(traits))

  if (length(traits) < 2L) stop("至少需要两个性状。", call. = FALSE)
  if (!all(traits %in% unique(marg$TRAIT))) {
    stop("部分 traits 不在 marginal_results 中。", call. = FALSE)
  }
  if (!all(traits %in% unique(cond$trait))) {
    stop("部分 traits 不在 conditional_results 中。", call. = FALSE)
  }
  if (is.null(loci)) loci <- unique(marg$SNPID)
  loci <- unique(as.character(loci))

  cond$sig_cond <- ifelse(
    !is.na(cond$sig_cond_input),
    cond$sig_cond_input,
    is.finite(cond$pval_cond) & cond$pval_cond < alpha2
  )

  marginal_sig <- function(snp, trait) {
    row <- marg[marg$SNPID == snp & marg$TRAIT == trait, , drop = FALSE]
    nrow(row) == 1L && is.finite(row$P_Value) && row$P_Value < alpha1
  }
  conditional_sig <- function(snp, trait) {
    row <- cond[cond$locus == snp & cond$trait == trait, , drop = FALSE]
    nrow(row) == 1L && isTRUE(row$sig_cond)
  }

  pairs <- utils::combn(traits, 2L, simplify = FALSE)
  rows <- vector("list", length(loci) * length(pairs))
  index <- 0L

  for (snp in loci) {
    for (pair in pairs) {
      m1 <- marginal_sig(snp, pair[1L])
      m2 <- marginal_sig(snp, pair[2L])
      if (!include_not_detected && !m1 && !m2) next

      key <- paste(pair, collapse = "__")
      mr_result <- mr_by_pair[[key]]
      if (is.null(mr_result)) {
        reverse_key <- paste(rev(pair), collapse = "__")
        mr_result <- mr_by_pair[[reverse_key]]
        if (!is.null(mr_result)) {
          # Standard run_all_trait_pairs_mr() always uses sorted pair order.
          # A reversed external entry cannot be safely relabeled automatically.
          stop("mr_by_pair 使用了反向 pair key；请统一为 traits 顺序的无序 pair key。", call. = FALSE)
        }
      }

      ann <- classify_pair_pattern(
        marg_sig_trait1 = m1,
        marg_sig_trait2 = m2,
        cond_sig_trait1 = conditional_sig(snp, pair[1L]),
        cond_sig_trait2 = conditional_sig(snp, pair[2L]),
        mr_results = mr_result,
        traits = pair,
        min_iv = min_iv,
        require_steiger = require_steiger,
        reject_heterogeneity = reject_heterogeneity,
        require_robust_direction = require_robust_direction,
        return_details = TRUE
      )

      index <- index + 1L
      rows[[index]] <- cbind(
        data.frame(SNPID = snp, stringsAsFactors = FALSE),
        ann
      )
    }
  }

  if (index == 0L) {
    return(data.frame(
      SNPID = character(0L), trait1 = character(0L), trait2 = character(0L),
      pattern = character(0L), direction = character(0L),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows[seq_len(index)])
}


# ------------------------------------------------------------------------------
# Backward compatibility for old two-trait scripts
# ------------------------------------------------------------------------------

#' Deprecated two-trait compatibility wrapper
#'
#' New code should use classify_pair_pattern(). The old n_marg_sig interface
#' cannot represent pair-specific marginal significance when m > 2.
#' @export
classify_locus <- function(n_marg_sig,
                           cond_sig_by_trait,
                           mr_results = NULL,
                           traits = NULL) {
  warning(
    "classify_locus() 已弃用；请使用 classify_pair_pattern()。",
    call. = FALSE
  )
  if (is.null(traits)) traits <- names(cond_sig_by_trait)
  if (length(traits) != 2L) {
    stop("兼容包装器仅支持两个性状。", call. = FALSE)
  }
  if (length(cond_sig_by_trait) != 2L) {
    stop("cond_sig_by_trait 必须包含两个性状。", call. = FALSE)
  }

  # Old input only provides a count, so the exact marginally significant trait
  # is unknowable when n_marg_sig == 1. The wrapper returns pattern1 directly.
  if (is.na(n_marg_sig) || n_marg_sig <= 0L) return("not_detected")
  if (n_marg_sig == 1L) return("pattern1")

  classify_pair_pattern(
    marg_sig_trait1 = TRUE,
    marg_sig_trait2 = TRUE,
    cond_sig_trait1 = isTRUE(cond_sig_by_trait[traits[1L]]),
    cond_sig_trait2 = isTRUE(cond_sig_by_trait[traits[2L]]),
    mr_results = mr_results,
    traits = traits,
    return_details = FALSE
  )
}
