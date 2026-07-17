# ==============================================================================
# classify.R — Five-class locus classification
# ==============================================================================


#' Classify a single locus into one of five mechanistic classes
#'
#' Integrates Layer 1 marginal significance, Layer 2 conditional projection,
#' and Layer 3 bidirectional MR to assign a mechanistic class to each locus,
#' following Table 1 of the paper.
#'
#' @param n_marg_sig Integer. Number of traits with significant marginal
#'   association at Layer 1. From l1$snp_summary$n_marg.
#' @param cond_sig_by_trait Named logical vector. Whether each trait retains
#'   significant conditional association after Layer 2 projection. Names must
#'   match trait names in \code{traits}.
#'   Example: \code{c(Trait1 = TRUE, Trait2 = FALSE)}.
#' @param mr_results Return value of \code{bidirectional_mr()}, or \code{NULL}.
#'   NULL triggers ablation mode (Layer 1 + 2 only, no MR).
#' @param traits Character vector of length 2. Trait names c("TraitA", "TraitB").
#'
#' @return Character string: one of \code{"null"}, \code{"class1"},
#'   \code{"class2"}, \code{"class3"}, \code{"class4"}, \code{"class5"}.
#'
#' @details
#' Classification logic:
#'
#' \strong{n_marg = 0}: null (not detected at Layer 1).
#'
#' \strong{n_marg = 1}: The locus is marginally significant for only one trait.
#'   This pattern covers two distinct biological scenarios:
#'   \itemize{
#'     \item True Class 1 (trait-specific): the locus directly affects only one
#'       trait, and the global MR is non-significant.
#'     \item True Class 3 (complete mediation): the locus directly affects the
#'       upstream trait only; the indirect effect on the downstream trait
#'       (tau * beta) is too weak to reach genome-wide significance at Layer 1.
#'       In this case the global MR A->B is significant, and the locus shows
#'       no significant conditional effect on the outcome trait.
#'   }
#'   Both scenarios produce n_marg = 1 at Layer 1. MR is therefore applied to
#'   n_marg = 1 loci as well, to rescue Class 3 loci that were missed at Layer 1.
#'   This is why Layer 3 is indispensable: without MR, all n_marg = 1 loci
#'   would be conservatively labelled Class 1.
#'
#' \strong{n_marg >= 2}: Pleiotropy candidate; enters full Layer 2 + Layer 3
#'   pipeline.
#'
#' \strong{Ablation mode} (\code{mr_results = NULL}): MR is skipped.
#'   n_marg = 1 -> Class 1 (cannot distinguish from Class 3).
#'   n_marg >= 2, n_cond >= 2 -> Class 2 (horizontal pleiotropy).
#'   n_marg >= 2, n_cond < 2  -> Class 1 (conservative downgrade).
#'   This misclassification of Class 3 loci in ablation mode demonstrates
#'   that directional causal evidence is indispensable.
#'
#' @export
classify_locus <- function(n_marg_sig,
                           cond_sig_by_trait,
                           mr_results = NULL,
                           traits     = NULL) {
  
  n_marg_sig <- as.integer(n_marg_sig)
  
  # ── Input validation ────────────────────────────────────────────────────────
  if (is.null(traits) && is.null(names(cond_sig_by_trait))) {
    stop("At least one of 'traits' or names(cond_sig_by_trait) must be provided.")
  }
  if (!is.null(traits) && length(traits) != 2L) {
    stop("'traits' must be a character vector of length 2.")
  }
  
  # ── Null: not detected at Layer 1 ──────────────────────────────────────────
  if (is.na(n_marg_sig) || n_marg_sig == 0L) {
    return("null")
  }
  
  # ── Internal helper: resolve outcome trait name ─────────────────────────────
  get_outcome <- function(ab_significant) {
    if (ab_significant) {
      if (!is.null(traits)) traits[2L] else names(cond_sig_by_trait)[2L]
    } else {
      if (!is.null(traits)) traits[1L] else names(cond_sig_by_trait)[1L]
    }
  }
  
  # ── n_marg = 1: trait-specific or complete mediation ───────────────────────
  #
  # Both Class 1 and Class 3 can produce n_marg = 1 at Layer 1:
  #   Class 1: locus affects one trait directly; MR non-significant.
  #   Class 3: locus affects upstream trait; indirect effect on downstream
  #            trait too weak for Layer 1 detection; MR A->B significant
  #            and outcome conditional effect non-significant.
  #
  # MR is applied here to distinguish the two cases.
  if (n_marg_sig == 1L) {
    
    # Ablation mode: cannot distinguish Class 1 vs Class 3 without MR
    if (is.null(mr_results)) return("class1")
    
    mr_AB_sig <- isTRUE(mr_results$AB$sig)
    mr_BA_sig <- isTRUE(mr_results$BA$sig)
    
    # MR non-significant in either direction -> Class 1
    if (!mr_AB_sig && !mr_BA_sig) return("class1")
    
    # MR significant in one direction: check outcome conditional effect
    #   outcome cond non-significant -> complete mediation -> Class 3
    #   outcome cond significant     -> locus has direct effect on outcome
    #                                   despite n_marg=1 (edge case) -> Class 1
    if (mr_AB_sig || mr_BA_sig) {
      outcome_trait    <- get_outcome(mr_AB_sig)
      outcome_cond_sig <- isTRUE(cond_sig_by_trait[outcome_trait])
      return(if (!outcome_cond_sig) "class3" else "class1")
    }
  }
  
  # ── n_marg >= 2: pleiotropy candidate ──────────────────────────────────────
  n_cond_sig <- sum(cond_sig_by_trait, na.rm = TRUE)
  
  # Ablation mode: Layer 1 + 2 only
  if (is.null(mr_results)) {
    return(if (n_cond_sig >= 2L) "class2" else "class1")
  }
  
  mr_AB_sig <- isTRUE(mr_results$AB$sig)
  mr_BA_sig <- isTRUE(mr_results$BA$sig)
  
  # Class 5: bidirectional MR significant (feedback loop or confounding)
  if (mr_AB_sig && mr_BA_sig) return("class5")
  
  # MR non-significant in either direction
  #   n_cond >= 2 -> Class 2 (horizontal pleiotropy, independent direct effects)
  #   n_cond <  2 -> Class 1 (conservative downgrade)
  if (!mr_AB_sig && !mr_BA_sig) {
    return(if (n_cond_sig >= 2L) "class2" else "class1")
  }
  
  # MR significant in one direction: distinguish complete vs partial mediation
  # Guard: if no conditional effects remain after projection but MR is
  # significant, data are contradictory; downgrade conservatively.
  if (n_cond_sig == 0L) return("class1")
  
  outcome_trait    <- get_outcome(mr_AB_sig)
  outcome_cond_sig <- isTRUE(cond_sig_by_trait[outcome_trait])
  
  # outcome cond significant   -> direct + indirect effect -> Class 4 (partial)
  # outcome cond non-significant -> indirect effect only   -> Class 3 (complete)
  if (outcome_cond_sig) "class4" else "class3"
}