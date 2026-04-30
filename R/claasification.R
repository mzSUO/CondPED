#' Classify QTL into trait-specific, independent, shared, or not detected categories
#'
#' Implements the bidirectional conditional testing framework (method.md Section 0.2.4).
#'
#' @param gwas_results Output from run_bidirectional_gwas()
#' @param alpha1 Significance threshold for marginal test
#' @param alpha2 Significance threshold for forward conditional test
#' @param alpha3 Significance threshold for reverse conditional test (default: 0.05)
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{SNP}{SNP name}
#'   \item{Trait}{Trait name}
#'   \item{classification}{Category: "not_detected", "shared", "independent", "trait_specific"}
#'   \item{pval_marginal}{Marginal p-value}
#'   \item{pval_cond}{Conditional p-value (forward)}
#' }
#'
#' @details
#' Classification logic (method.md Section 0.2.4):
#'
#' | Category       | Criteria                                                                      |
#' |----------------|-------------------------------------------------------------------------------|
#' | Not detected   | P_marginal >= α1                                                               |
#' | Shared         | P_marginal < α1, P_cond >= α2                                                 |
#' | Independent    | P_cond < α2, exists t: P_reverse < α3                                         |
#' | Trait-specific | P_cond < α2, for all t != i: P_reverse >= α3                                 |
#'
#' @export
#' @importFrom dplyr filter select rename
#'
classify_qtl <- function(gwas_results,
                         alpha1 = NULL,
                         alpha2 = NULL,
                         alpha3 = 0.05) {

  # Extract results from gwas_results
  marginal <- gwas_results$step1_marginal
  forward <- gwas_results$step2_forward
  reverse <- gwas_results$step3_reverse

  # Get marginal GWAS results
  marg_df <- marginal$gwas_results

  # Get number of traits and candidate SNPs
  traits <- unique(marg_df$TRAIT)
  snps <- unique(marg_df$SNP)
  n_trait <- length(traits)
  n_snp <- length(snps)

  # Set default thresholds
  if (is.null(alpha1)) alpha1 <- 0.05 / n_snp
  if (is.null(alpha2)) alpha2 <- 0.05 / (n_snp * n_trait)

  # Initialize results data frame
  results <- data.frame(
    SNP = character(),
    Trait = character(),
    classification = character(),
    pval_marginal = numeric(),
    pval_cond = numeric(),
    stringsAsFactors = FALSE
  )

  # Process each trait
  for (trait in traits) {
    # Get marginal p-values for this trait
    marg_trait <- marg_df[marg_df$TRAIT == trait, ]

    # Get forward conditional results for this trait
    fwd_trait <- forward[[trait]]

    # Get significant SNPs from forward step
    sig_fwd <- fwd_trait[fwd_trait$P_Value < alpha2, ]

    # Process each SNP
    for (snp in snps) {
      # Get marginal p-value
      marg_row <- marg_trait[marg_trait$SNP == snp, ]
      p_marg <- if (nrow(marg_row) > 0) marg_row$P_Value else NA

      # Get conditional p-value
      fwd_row <- fwd_trait[fwd_trait$SNP == snp, ]
      p_cond <- if (nrow(fwd_row) > 0) fwd_row$P_Value else NA

      # Classification
      if (is.na(p_marg) || p_marg >= alpha1) {
        # Not detected
        classification <- "not_detected"
      } else if (is.na(p_cond) || p_cond >= alpha2) {
        # Shared: significant in marginal but not in conditional
        classification <- "shared"
      } else {
        # Significant in conditional - need to check reverse
        # Get reverse results for this SNP-trait pair
        is_trait_specific <- TRUE

        for (test_trait in traits) {
          if (test_trait != trait) {
            key <- paste0(test_trait, "_excl_", trait)

            if (!is.null(reverse[[key]])) {
              rev_df <- reverse[[key]]
              rev_row <- rev_df[rev_df$SNP == snp, ]

              if (nrow(rev_row) > 0) {
                p_rev <- rev_row$P_Value

                if (!is.na(p_rev) && p_rev < alpha3) {
                  # Found a trait where reverse is significant
                  is_trait_specific <- FALSE
                  break
                }
              }
            }
          }
        }

        classification <- if (is_trait_specific) "trait_specific" else "independent"
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
                                  levels = c("not_detected", "shared", "independent", "trait_specific"))

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
