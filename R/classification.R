#' Classify QTL into five categories based on conditional analysis
#'
#' Implements the four-layer testing framework (method.md Section 0.2.3).
#'
#' @param gwas_results Output from run_bidirectional_gwas()
#' @param alpha1 Significance threshold for marginal test (default: 0.05)
#' @param alpha2 Significance threshold for forward conditional test
#' @param mr_results Optional output from bidirectional_causal_test() for vertical pleiotropy classification
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{SNP}{SNP name}
#'   \item{Trait}{Trait name}
#'   \item{classification}{Category: "not_detected", "covariance_induced", "trait_specific", "horizontal_pleiotropy", "vertical_pleiotropy"}
#'   \item{pval_marginal}{Marginal p-value}
#'   \item{pval_cond}{Conditional p-value}
#' }
#'
#' @details
#' Five-category classification (method.md Table 1):
#'
#' | Category              | Criteria                                                                      |
#' |----------------------|-------------------------------------------------------------------------------|
#' | not_detected        | P_marginal >= α1                                                               |
#' | covariance_induced   | P_marginal < α1, P_cond >= α2                                                 |
#' | trait_specific      | P_cond < α2, for all t != i: P_cond(t) >= α2                                |
#' | horizontal_pleiotropy | P_cond < α2 for >= 2 traits, no causal evidence                             |
#' | vertical_pleiotropy | P_cond < α2 for >= 2 traits, MR supports causal direction                   |
#'
#' @export
#' @importFrom dplyr filter select rename
#'
classify_qtl <- function(gwas_results,
                         alpha1 = 0.05,
                         alpha2 = NULL,
                         mr_results = NULL) {

  # Extract results from gwas_results
  marginal <- gwas_results$step1_marginal
  forward <- gwas_results$step2_forward

  # Get marginal GWAS results
  marg_df <- marginal$gwas_results

  # Get number of traits and candidate SNPs
  traits <- unique(marg_df$TRAIT)
  snps <- unique(marg_df$SNP)
  n_trait <- length(traits)
  n_snp <- length(snps)

  # Set default thresholds
  if (is.null(alpha2)) alpha2 <- 0.05 / (n_snp * n_trait)

  # Build conditional p-value matrix
  cond_pvals <- matrix(NA, nrow = n_snp, ncol = n_trait,
                       dimnames = list(snps, traits))
  for (trait in traits) {
    fwd_trait <- forward[[trait]]
    if (!is.null(fwd_trait) && nrow(fwd_trait) > 0) {
      cond_pvals[fwd_trait$SNP, trait] <- fwd_trait$P_Value
    }
  }

  # Initialize results data frame
  results <- data.frame(
    SNP = character(),
    Trait = character(),
    classification = character(),
    pval_marginal = numeric(),
    pval_cond = numeric(),
    stringsAsFactors = FALSE
  )

  # Process each SNP-Trait pair
  for (snp in snps) {
    for (trait in traits) {
      # Get marginal p-value
      marg_row <- marg_df[marg_df$SNP == snp & marg_df$TRAIT == trait, ]
      p_marg <- if (nrow(marg_row) > 0) marg_row$P_Value else NA

      # Get conditional p-value
      p_cond <- cond_pvals[snp, trait]

      # ===== Layer 1: Not detected =====
      if (is.na(p_marg) || p_marg >= alpha1) {
        classification <- "not_detected"
      }
      # ===== Layer 2: Covariance-induced =====
      else if (is.na(p_cond) || p_cond >= alpha2) {
        classification <- "covariance_induced"
      }
      # ===== Layer 3 & 4: Check if pleiotropic =====
      else {
        # Count how many traits have significant conditional effects for this SNP
        n_sig_traits <- sum(cond_pvals[snp, ] < alpha2, na.rm = TRUE)

        if (n_sig_traits == 1) {
          # Only this trait is significant -> trait_specific
          classification <- "trait_specific"
        } else {
          # Multiple traits significant -> need MR to distinguish horizontal/vertical
          if (!is.null(mr_results)) {
            # Check if MR supports causal direction involving this SNP
            is_vertical <- FALSE

            # Check A -> B direction
            for (trait_a in traits) {
              for (trait_b in traits) {
                if (trait_a != trait_b) {
                  # Check if this SNP is an IV for A -> B
                  ivs_a_to_b <- mr_results[[paste0(trait_a, "_to_", trait_b)]]$iv_A_to_B
                  if (!is.null(ivs_a_to_b) && snp %in% ivs_a_to_b) {
                    # Check if causal effect is significant
                    mr_result <- mr_results[[paste0(trait_a, "_to_", trait_b)]]$A_to_B
                    if (!is.null(mr_result) && !is.na(mr_result$p_value) && mr_result$p_value < 0.05) {
                      is_vertical <- TRUE
                      break
                    }
                  }
                }
              }
              if (is_vertical) break
            }

            classification <- if (is_vertical) "vertical_pleiotropy" else "horizontal_pleiotropy"
          } else {
            # No MR results -> default to horizontal_pleiotropy
            classification <- "horizontal_pleiotropy"
          }
        }
      }

      # Add to results
      results <- rbind(results, data.frame(
        SNP = snp,
        Trait = trait,
        classification = classification,
        pval_marginal = p_marg,
        pval_cond = p_cond,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Convert classification to factor with defined levels
  results$classification <- factor(results$classification,
                                  levels = c("not_detected", "covariance_induced",
                                            "trait_specific", "horizontal_pleiotropy",
                                            "vertical_pleiotropy"))

  return(results)
}


#' Summarize QTL classification results
#'
#' @param classification_results Output from classify_qtl()
#'
#' @return Data frame with counts for each category
#' @export
#'
summarize_classification <- function(classification_results) {

  summary_df <- table(classification_results$classification, classification_results$Trait)

  # Convert to data frame
  result <- as.data.frame(summary_df)
  colnames(result) <- c("Category", "Trait", "Count")

  return(result)
}
