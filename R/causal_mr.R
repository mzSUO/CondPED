#' CondPED: Causal Inference via Mendelian Randomization
#'
#' This module implements Mendelian randomization for trait-level causal inference
#' using conditional independence as instrument selection criterion.
#'
#' @name causal-functions
#' @docType package
NULL


#' Select valid instrumental variables for MR
#'
#' Implements the instrument selection criterion from method.md Section 2.5.2.
#' An SNP is selected as a valid IV if:
#' 1. Relevance: P_marginal < alpha_rel
#' 2. Exclusion: P_cond > alpha_excl (no direct effect on outcome)
#' 3. LD pruning: |r| < r_LD
#'
#' @param gwas_results Output from run_bidirectional_gwas()
#' @param exposure_trait Name of exposure trait
#' @param outcome_trait Name of outcome trait
#' @param alpha_rel Relevance threshold (default: 5e-8)
#' @param alpha_excl Exclusion threshold (default: 0.05)
#' @param r_LD LD pruning threshold (default: 0.1)
#' @param G Genotype matrix for LD calculation
#'
#' @return List with:
#'   \item{iv_selected}{Character vector of selected IVs}
#'   \item{snp_effects}{Data frame with SNP effects on exposure and outcome}
#'
#' @export
#'
select_ivs <- function(gwas_results,
                       exposure_trait,
                       outcome_trait,
                       alpha_rel = 5e-8,
                       alpha_excl = 0.05,
                       r_LD = 0.1,
                       G = NULL) {

  # Get marginal GWAS results
  marg_df <- gwas_results$step1_marginal$gwas_results

  # Filter for exposure and outcome traits
  exposure_df <- marg_df[marg_df$TRAIT == exposure_trait, ]
  outcome_df <- marg_df[marg_df$TRAIT == outcome_trait, ]

  if (nrow(exposure_df) == 0 || nrow(outcome_df) == 0) {
    stop("Traits not found in GWAS results")
  }

  # Get forward conditional results for exclusion criterion
  forward_df <- gwas_results$step2_forward[[outcome_trait]]

  # Step 1: Relevance - select SNPs associated with exposure
  # (For now, use marginal p-value as proxy since we have full GWAS)
  relevant_snps <- exposure_df$SNP[exposure_df$P_Value < alpha_rel]

  if (length(relevant_snps) == 0) {
    warning("No SNPs meeting relevance criterion")
    return(list(
      iv_selected = character(),
      snp_effects = data.frame()
    ))
  }

  # Step 2: Exclusion - check if SNPs have direct effect on outcome
  # Use conditional p-value: if P_cond > alpha_excl, no direct effect
  selected_ivs <- character()

  for (snp in relevant_snps) {
    if (!is.null(forward_df)) {
      fwd_row <- forward_df[forward_df$SNP == snp, ]
      if (nrow(fwd_row) > 0) {
        p_cond <- fwd_row$P_Value
        # If conditional effect is not significant, satisfies exclusion
        if (p_cond >= alpha_excl) {
          selected_ivs <- c(selected_ivs, snp)
        }
      }
    } else {
      # If no conditional results, use marginal as fallback
      out_row <- outcome_df[outcome_df$SNP == snp, ]
      if (nrow(out_row) > 0 && out_row$P_Value >= alpha_excl) {
        selected_ivs <- c(selected_ivs, snp)
      }
    }
  }

  if (length(selected_ivs) == 0) {
    warning("No SNPs meeting exclusion criterion")
    return(list(
      iv_selected = character(),
      snp_effects = data.frame()
    ))
  }

  # Step 3: LD pruning (if genotype matrix provided)
  if (!is.null(G) && length(selected_ivs) > 1) {
    G_sub <- G[, selected_ivs, drop = FALSE]

    # Calculate LD matrix
    if (ncol(G_sub) > 1) {
      cor_matrix <- cor(G_sub)
      cor_matrix[is.na(cor_matrix)] <- 0

      # Prune: keep SNP with lowest p-value, remove correlated ones
      pruned <- selected_ivs[1]
      for (snp in selected_ivs[-1]) {
        # Check correlation with already kept SNPs
        keep <- TRUE
        for (kept_snp in pruned) {
          if (!is.na(cor_matrix[snp, kept_snp]) &&
              abs(cor_matrix[snp, kept_snp]) > r_LD) {
            keep <- FALSE
            break
          }
        }
        if (keep) {
          pruned <- c(pruned, snp)
        }
      }
      selected_ivs <- pruned
    }
  }

  # Build SNP effects data frame
  effects_df <- data.frame(
    SNP = selected_ivs,
    stringsAsFactors = FALSE
  )

  # Add exposure effects
  for (snp in selected_ivs) {
    exp_row <- exposure_df[exposure_df$SNP == snp, ]
    if (nrow(exp_row) > 0) {
      effects_df$Effect_exposure[effects_df$SNP == snp] <- exp_row$A
    }
  }

  # Add outcome effects
  for (snp in selected_ivs) {
    out_row <- outcome_df[outcome_df$SNP == snp, ]
    if (nrow(out_row) > 0) {
      effects_df$Effect_outcome[effects_df$SNP == snp] <- out_row$A
    }
  }

  return(list(
    iv_selected = selected_ivs,
    snp_effects = effects_df
  ))
}


#' Estimate causal effect using Wald ratio
#'
#' @param snp_effects Data frame with SNP effects on exposure and outcome
#'
#' @return List with causal estimate and standard error
#' @export
#'
estimate_wald_ratio <- function(snp_effects) {

  if (nrow(snp_effects) == 0) {
    return(list(
      causal_effect = NA,
      se = NA,
      p_value = NA,
      n_iv = 0
    ))
  }

  # For single SNP: Wald ratio = beta_outcome / beta_exposure
  if (nrow(snp_effects) == 1) {
    beta_exp <- snp_effects$Effect_exposure[1]
    beta_out <- snp_effects$Effect_outcome[1]

    if (is.na(beta_exp) || beta_exp == 0) {
      return(list(
        causal_effect = NA,
        se = NA,
        p_value = NA,
        n_iv = 1
      ))
    }

    # Delta method for SE
    gamma <- beta_out / beta_exp
    # Approximate SE (simplified)
    se <- abs(gamma) / abs(beta_exp)

    # Z-test
    z <- gamma / se
    p_value <- 2 * pnorm(-abs(z))

    return(list(
      causal_effect = gamma,
      se = se,
      p_value = p_value,
      n_iv = 1
    ))
  }

  # Multiple SNPs: use inverse-variance weighted (IVW) method
  beta_exp <- snp_effects$Effect_exposure
  beta_out <- snp_effects$Effect_outcome

  # Remove NAs
  valid <- !is.na(beta_exp) & !is.na(beta_out) & beta_exp != 0
  beta_exp <- beta_exp[valid]
  beta_out <- beta_out[valid]

  if (length(beta_exp) == 0) {
    return(list(
      causal_effect = NA,
      se = NA,
      p_value = NA,
      n_iv = 0
    ))
  }

  # IVW estimate: weighted mean of Wald ratios
  weights <- 1 / (beta_exp^2)
  gamma_ivw <- sum(weights * beta_out / beta_exp) / sum(weights)

  # SE (IVW)
  resid <- beta_out - gamma_ivw * beta_exp
  k <- length(beta_exp)
  if (k > 1) {
    se_ivw <- sqrt(sum(weights * resid^2) / ((k - 1) * sum(weights)))
  } else {
    se_ivw <- NA
  }

  # Z-test
  if (!is.na(se_ivw) && se_ivw > 0) {
    z <- gamma_ivw / se_ivw
    p_value <- 2 * pnorm(-abs(z))
  } else {
    p_value <- NA
  }

  return(list(
    causal_effect = gamma_ivw,
    se = se_ivw,
    p_value = p_value,
    n_iv = length(beta_exp)
  ))
}


#' Bidirectional causal estimation
#'
#' Tests both causal directions: A -> B and B -> A
#'
#' @param gwas_results Output from run_bidirectional_gwas()
#' @param trait_a Name of first trait
#' @param trait_b Name of second trait
#' @param G Genotype matrix for LD pruning (optional)
#'
#' @return List with causal estimates for both directions
#' @export
#'
bidirectional_causal_test <- function(gwas_results,
                                    trait_a,
                                    trait_b,
                                    G = NULL) {

  # Direction A -> B
  cat("Testing direction:", trait_a, "->", trait_b, "\n")
  iv_a_to_b <- select_ivs(gwas_results, trait_a, trait_b, G = G)

  if (length(iv_a_to_b$iv_selected) > 0) {
    result_a_to_b <- estimate_wald_ratio(iv_a_to_b$snp_effects)
  } else {
    result_a_to_b <- list(
      causal_effect = NA,
      se = NA,
      p_value = NA,
      n_iv = 0
    )
  }

  # Direction B -> A
  cat("Testing direction:", trait_b, "->", trait_a, "\n")
  iv_b_to_a <- select_ivs(gwas_results, trait_b, trait_a, G = G)

  if (length(iv_b_to_a$iv_selected) > 0) {
    result_b_to_a <- estimate_wald_ratio(iv_b_to_a$snp_effects)
  } else {
    result_b_to_a <- list(
      causal_effect = NA,
      se = NA,
      p_value = NA,
      n_iv = 0
    )
  }

  return(list(
    A_to_B = result_a_to_b,
    B_to_A = result_b_to_a,
    iv_A_to_B = iv_a_to_b$iv_selected,
    iv_B_to_A = iv_b_to_a$iv_selected
  ))
}


#' Heterogeneity test for MR
#'
#' Tests whether IVs provide consistent estimates (Cochran's Q)
#'
#' @param snp_effects Data frame with SNP effects
#'
#' @return List with Q statistic and p-value
#' @export
#'
test_heterogeneity <- function(snp_effects) {

  beta_exp <- snp_effects$Effect_exposure
  beta_out <- snp_effects$Effect_outcome

  # Remove NAs
  valid <- !is.na(beta_exp) & !is.na(beta_out) & beta_exp != 0
  beta_exp <- beta_exp[valid]
  beta_out <- beta_out[valid]

  k <- length(beta_exp)

  if (k < 2) {
    return(list(
      Q = NA,
      df = NA,
      p_value = NA
    ))
  }

  # Calculate IVW estimate
  weights <- 1 / (beta_exp^2)
  gamma_ivw <- sum(weights * beta_out / beta_exp) / sum(weights)

  # Calculate Q statistic
  Q <- sum(weights * (beta_out - gamma_ivw * beta_exp)^2)
  df <- k - 1
  p_value <- 1 - pchisq(Q, df)

  return(list(
    Q = Q,
    df = df,
    p_value = p_value
  ))
}
