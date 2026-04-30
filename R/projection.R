#' projection.R
#' Conditional Projection Framework for Multi-trait Analysis
#' 
#' This module implements orthogonal projection to construct conditional phenotypes
#' that isolate trait-specific genetic variation by removing linear dependence 
#' on other traits.

# ==============================================================================
# Core Functions
# ==============================================================================

#' Construct conditional phenotype via orthogonal projection
#'
#' Implements Equation 4 from the manuscript:
#' \deqn{\tilde{y}_{i}^* = \tilde{y}_{i} - \tilde{\mathbf{y}}_{-i}^T \gamma_{i,-i}}
#' where \eqn{\gamma_{i,-i} = V_{-i}^{-1} C_{-i,i}}
#'
#' @param pheno_data a data frame with columns: GID (character/factor) and trait columns (numeric)
#' @param target_trait character, name of the target trait to be conditioned
#' @param adjust_covariates logical, whether to adjust for covariates before projection (default: FALSE)
#' @param covariate_cols character vector, names of covariate columns if adjust_covariates = TRUE
#'
#' @return a list with:
#'   \item{pheno_cond}{data frame with GID and conditional phenotype}
#'   \item{gamma}{projection coefficient vector}
#'   \item{var_original}{original phenotypic variance}
#'   \item{var_conditional}{conditional phenotypic variance (equals V_{i|-i})}
#'   \item{var_reduction}{proportion of variance explained by other traits}
#'   \item{conditioning_traits}{character vector of traits used for conditioning}
#'
#' @export
#' @import dplyr
#'
#' @examples 
#' \dontrun{
#' # Basic usage
#' cond_result <- conditional_phenotype(pheno_data, "Height")
#' 
#' # With covariate adjustment
#' cond_result <- conditional_phenotype(pheno_data, "Height", 
#'                                     adjust_covariates = TRUE,
#'                                     covariate_cols = c("Sex", "Age"))
#' }
conditional_phenotype <- function(pheno_data, 
                                 target_trait,
                                 adjust_covariates = FALSE,
                                 covariate_cols = NULL)
{
  ## Validate input
  if(!target_trait %in% colnames(pheno_data)) {
    stop(sprintf("Error: target trait '%s' not found in pheno_data", target_trait))
  }
  
  if(!"GID" %in% colnames(pheno_data)) {
    stop("Error: pheno_data must contain a 'GID' column")
  }
  
  ## Extract trait columns
  exclude_cols <- c("GID", covariate_cols)
  trait_cols <- setdiff(colnames(pheno_data), exclude_cols)
  
  if(!target_trait %in% trait_cols) {
    stop("Error: target_trait must be a trait column, not GID or covariate")
  }
  
  ## Step 1: Adjust for covariates if requested
  pheno_resid <- pheno_data
  
  if(adjust_covariates) {
    if(is.null(covariate_cols)) {
      stop("Error: covariate_cols must be provided when adjust_covariates = TRUE")
    }
    
    for(trait in trait_cols) {
      formula_str <- paste(trait, "~", paste(covariate_cols, collapse = " + "))
      lm_fit <- lm(as.formula(formula_str), data = pheno_data)
      pheno_resid[[trait]] <- residuals(lm_fit)
    }
  }
  
  ## Step 2: Prepare matrices
  other_traits <- setdiff(trait_cols, target_trait)
  
  if(length(other_traits) == 0) {
    stop("Error: no other traits available for conditioning")
  }
  
  y_target <- pheno_resid[[target_trait]]
  Y_others <- as.matrix(pheno_resid[, other_traits, drop = FALSE])
  
  ## Step 3: Calculate projection coefficients (Equation 4)
  # γ_{i,-i} = V_{-i}^{-1} C_{-i,i}
  V_others <- cov(Y_others, use = "complete.obs")
  C_others_target <- cov(Y_others, y_target, use = "complete.obs")
  
  ## Handle singular covariance matrix
  gamma <- tryCatch({
    solve(V_others) %*% C_others_target
  }, error = function(e) {
    warning("Covariance matrix V_{-i} is singular. Using Moore-Penrose pseudo-inverse.")
    MASS::ginv(V_others) %*% C_others_target
  })
  
  ## Step 4: Construct conditional phenotype
  # ỹ*_i = ỹ_i - ỹ_{-i}^T γ_{i,-i}
  y_cond <- y_target - Y_others %*% gamma
  y_cond <- as.vector(y_cond)
  
  ## Step 5: Compute variance components
  var_original <- var(y_target, na.rm = TRUE)
  var_conditional <- var(y_cond, na.rm = TRUE)
  var_reduction <- 1 - var_conditional / var_original
  
  ## Step 6: Prepare output
  pheno_cond <- data.frame(
    GID = pheno_data$GID,
    trait = y_cond
  )
  colnames(pheno_cond)[2] <- paste0(target_trait, "_cond")
  
  ## Return results
  return(list(
    pheno_cond = pheno_cond,
    gamma = gamma,
    var_original = var_original,
    var_conditional = var_conditional,
    var_reduction = var_reduction,
    conditioning_traits = other_traits
  ))
}


#' Construct reverse conditional phenotype
#'
#' Condition target trait on all traits EXCEPT the exclusion trait.
#' Used in reverse conditional GWAS to test mediation hypothesis.
#'
#' @param pheno_data a data frame with GID and trait columns
#' @param target_trait character, name of the trait to be conditioned
#' @param exclude_trait character, name of the trait to exclude from conditioning set
#' @param adjust_covariates logical, whether to adjust for covariates (default: FALSE)
#' @param covariate_cols character vector, names of covariate columns
#'
#' @return same structure as conditional_phenotype()
#' @export
#' @import dplyr
#'
#' @examples 
#' \dontrun{
#' # Test if QTL effect on Yield is mediated through Height
#' reverse_result <- reverse_conditional(pheno_data, 
#'                                      target_trait = "Yield",
#'                                      exclude_trait = "Height")
#' }
reverse_conditional <- function(pheno_data, 
                               target_trait,
                               exclude_trait,
                               adjust_covariates = FALSE,
                               covariate_cols = NULL)
{
  ## Validate input
  exclude_cols <- c("GID", covariate_cols)
  trait_cols <- setdiff(colnames(pheno_data), exclude_cols)
  
  if(!target_trait %in% trait_cols) {
    stop(sprintf("Error: target trait '%s' not found", target_trait))
  }
  
  if(!exclude_trait %in% trait_cols) {
    stop(sprintf("Error: exclude trait '%s' not found", exclude_trait))
  }
  
  if(target_trait == exclude_trait) {
    stop("Error: target_trait and exclude_trait cannot be the same")
  }
  
  ## Construct conditioning set: all traits except target and exclude
  conditioning_set <- setdiff(trait_cols, c(target_trait, exclude_trait))
  
  if(length(conditioning_set) == 0) {
    stop("Error: no traits remaining for conditioning after exclusion")
  }
  
  ## Create subset data
  pheno_subset <- pheno_data[, c("GID", target_trait, conditioning_set, covariate_cols)]
  
  ## Call conditional_phenotype on the subset
  result <- conditional_phenotype(pheno_subset, 
                                 target_trait,
                                 adjust_covariates = adjust_covariates,
                                 covariate_cols = covariate_cols)
  
  ## Update output name to reflect reverse conditioning
  colnames(result$pheno_cond)[2] <- paste0(target_trait, "_cond_excl_", exclude_trait)
  
  return(result)
}


# ==============================================================================
# Diagnostic Functions
# ==============================================================================

#' Check orthogonality property of conditional phenotype
#'
#' Verifies Equation 6: Cov(ỹ*_i, ỹ_{-i}) = 0
#' and Equation 7: Var(ỹ*_i) = V_{i|-i}
#'
#' @param pheno_original original phenotype data frame
#' @param pheno_conditional conditional phenotype (output from conditional_phenotype)
#' @param target_trait character, name of the target trait
#' @param tolerance numeric, tolerance for near-zero covariance (default: 1e-10)
#'
#' @return a list with:
#'   \item{orthogonality_check}{logical vector, TRUE if covariance < tolerance}
#'   \item{covariances}{numeric vector of covariances between y*_i and each y_{-i}}
#'   \item{variance_check}{logical, whether Var(y*_i) matches theoretical V_{i|-i}}
#'   \item{var_empirical}{empirical variance of conditional phenotype}
#'   \item{var_theoretical}{theoretical conditional variance}
#'   \item{diagnostic_message}{character, summary message}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cond_result <- conditional_phenotype(pheno_data, "Height")
#' diag <- check_orthogonality(pheno_data, cond_result, "Height")
#' print(diag$diagnostic_message)
#' }
check_orthogonality <- function(pheno_original, 
                               pheno_conditional,
                               target_trait,
                               tolerance = 1e-10)
{
  ## Extract conditioning traits
  conditioning_traits <- pheno_conditional$conditioning_traits
  
  ## Get conditional phenotype values
  y_cond <- pheno_conditional$pheno_cond[, 2]
  
  ## Get original other traits
  Y_others <- as.matrix(pheno_original[, conditioning_traits, drop = FALSE])
  
  ## Compute covariances
  covs <- apply(Y_others, 2, function(y_other) {
    cov(y_cond, y_other, use = "complete.obs")
  })
  
  ## Check orthogonality
  ortho_check <- abs(covs) < tolerance
  
  ## Check variance equivalence
  var_empirical <- pheno_conditional$var_conditional
  var_theoretical <- pheno_conditional$var_original - 
                    sum(pheno_conditional$gamma * cov(Y_others, pheno_original[[target_trait]]))
  
  var_check <- abs(var_empirical - var_theoretical) < tolerance
  
  ## Generate diagnostic message
  if(all(ortho_check) & var_check) {
    msg <- "✓ Projection is valid: orthogonality and variance equivalence satisfied"
  } else {
    msg <- sprintf("⚠ Warning: projection may be inaccurate\n  - Orthogonality: %d/%d traits\n  - Variance equivalence: %s",
                   sum(ortho_check), length(ortho_check),
                   ifelse(var_check, "passed", "failed"))
  }
  
  ## Return results
  return(list(
    orthogonality_check = ortho_check,
    covariances = covs,
    variance_check = var_check,
    var_empirical = var_empirical,
    var_theoretical = var_theoretical,
    diagnostic_message = msg
  ))
}


#' Compute variance reduction matrix across all trait pairs
#'
#' Calculates the proportion of variance explained when conditioning 
#' each trait on every other single trait.
#'
#' @param pheno_data a data frame with GID and trait columns
#' @param adjust_covariates logical, whether to adjust for covariates
#' @param covariate_cols character vector, covariate column names
#'
#' @return a matrix where element [i,j] is the variance reduction when 
#'         conditioning trait i on trait j
#' @export
#'
#' @examples
#' \dontrun{
#' vr_matrix <- compute_variance_reduction_matrix(pheno_data)
#' heatmap(vr_matrix, main = "Variance Reduction Matrix")
#' }
compute_variance_reduction_matrix <- function(pheno_data,
                                              adjust_covariates = FALSE,
                                              covariate_cols = NULL)
{
  ## Get trait columns
  exclude_cols <- c("GID", covariate_cols)
  trait_cols <- setdiff(colnames(pheno_data), exclude_cols)
  n_traits <- length(trait_cols)
  
  ## Initialize matrix
  vr_matrix <- matrix(NA, nrow = n_traits, ncol = n_traits,
                     dimnames = list(trait_cols, trait_cols))
  
  ## Fill diagonal with 0 (no reduction when no conditioning)
  diag(vr_matrix) <- 0
  
  ## Compute variance reduction for each pair
  for(i in 1:n_traits) {
    for(j in 1:n_traits) {
      if(i != j) {
        ## Create temporary data with only two traits
        pheno_temp <- pheno_data[, c("GID", trait_cols[i], trait_cols[j], covariate_cols)]
        
        ## Condition trait i on trait j
        cond_result <- conditional_phenotype(pheno_temp, 
                                            trait_cols[i],
                                            adjust_covariates = adjust_covariates,
                                            covariate_cols = covariate_cols)
        
        vr_matrix[i, j] <- cond_result$var_reduction
      }
    }
  }
  
  return(vr_matrix)
}


# ==============================================================================
# Utility Functions
# ==============================================================================

#' Print summary of conditional projection result
#'
#' @param cond_result output from conditional_phenotype() or reverse_conditional()
#'
#' @return NULL (prints to console)
#' @export
print_projection_summary <- function(cond_result)
{
  cat("\n")
  cat("═══════════════════════════════════════════════════════\n")
  cat("  Conditional Phenotype Projection Summary\n")
  cat("═══════════════════════════════════════════════════════\n")
  cat(sprintf("Target trait:        %s\n", 
              gsub("_cond.*", "", colnames(cond_result$pheno_cond)[2])))
  cat(sprintf("Conditioning traits: %s\n", 
              paste(cond_result$conditioning_traits, collapse = ", ")))
  cat(sprintf("Original variance:   %.4f\n", cond_result$var_original))
  cat(sprintf("Conditional variance: %.4f\n", cond_result$var_conditional))
  cat(sprintf("Variance reduction:  %.2f%%\n", cond_result$var_reduction * 100))
  cat("═══════════════════════════════════════════════════════\n")
  cat("\n")
}


#' Batch construction of conditional phenotypes for all traits
#'
#' Convenience wrapper to construct conditional phenotypes for multiple traits
#'
#' @param pheno_data a data frame with GID and trait columns
#' @param adjust_covariates logical, whether to adjust for covariates
#' @param covariate_cols character vector, covariate column names
#' @param verbose logical, whether to print progress (default: TRUE)
#'
#' @return a named list where each element is the output of conditional_phenotype()
#' @export
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' all_cond <- batch_conditional_phenotypes(pheno_data)
#' 
#' # Access conditional phenotype for Height
#' height_cond <- all_cond$Height$pheno_cond
#' }
batch_conditional_phenotypes <- function(pheno_data,
                                        adjust_covariates = FALSE,
                                        covariate_cols = NULL,
                                        verbose = TRUE)
{
  ## Get trait columns
  exclude_cols <- c("GID", covariate_cols)
  trait_cols <- setdiff(colnames(pheno_data), exclude_cols)
  
  ## Initialize result list
  results <- list()
  
  ## Loop through traits
  for(trait in trait_cols) {
    if(verbose) {
      cat(sprintf("→ Constructing conditional phenotype for: %s\n", trait))
    }
    
    results[[trait]] <- conditional_phenotype(pheno_data, 
                                             trait,
                                             adjust_covariates = adjust_covariates,
                                             covariate_cols = covariate_cols)
    
    if(verbose) {
      cat(sprintf("  ✓ Variance reduction: %.2f%%\n\n", 
                  results[[trait]]$var_reduction * 100))
    }
  }
  
  return(results)
}