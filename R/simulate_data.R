#' Simulate Multi-Trait Phenotype and Genotype Data
#'
#' Generates synthetic genotype and phenotype data with known genetic architectures
#' for testing the conditional projection framework.
#'
#' This is the main user-facing function that combines:
#' \itemize{
#'   \item Generate genotype matrix
#'   \item Generate true genetic effects (with 5-class SNP classification)
#'   \item Generate phenotypes with causal pathway structure
#' }
#'
#' @param n_ind Integer. Number of individuals (default: 500)
#' @param n_trait Integer. Number of traits (default: 3)
#' @param n_snp Integer. Number of SNPs (default: 200)
#' @param n_env Integer. Number of environments (default: 1)
#' @param maf_range Numeric vector of length 2. Range for minor allele frequency (default: c(0.05, 0.5))
#' @param trait_cor Numeric. Target trait correlation for AR(1) structure (default: 0.5)
#' @param heritability Numeric. Heritability for each trait (default: 0.5)
#' @param effect_config List. Configuration for genetic effect patterns. See Details.
#' @param causal_matrix Matrix. Causal effect matrix (n_trait x n_trait). 
#'   causal_matrix[i,j] represents causal effect of trait j on trait i.
#'   If NULL, uses default_causal_matrix().
#'
#' @details
#' The \code{effect_config} parameter controls the genetic architecture with 5 SNP classes:
#' \describe{
#'   \item{covariance_induced}{SNPs with effects only on trait 1, creating apparent 
#'     associations on trait 2 through covariance structure}
#'   \item{trait_specific}{SNPs with effects only on trait 2 (trait-specific)}
#'   \item{horizontal}{SNPs with independent effects on multiple traits simultaneously}
#'   \item{vertical}{SNPs with effects only on upstream trait 1, propagating to trait 2 
#'     through causal_matrix}
#'   \item{bidirectional}{SNPs with effects on multiple traits including potential feedback}
#' }
#'
#' The \code{causal_matrix} parameter specifies causal pathways between traits:
#' Y_i = G * beta_i + sum_j(causal_matrix[i,j] * Y_j) + E_i
#'
#' @return List containing:
#' \describe{
#'   \item{G}{Genotype matrix (n_ind x n_snp), standardized {-1, 0, 1}}
#'   \item{Y}{Phenotype matrix (n_ind x n_trait)}
#'   \item{Y_base}{Base phenotype without causal propagation (n_ind x n_trait)}
#'   \item{E}{Environmental residuals (n_ind x n_trait)}
#'   \item{beta_true}{True total effect matrix (n_snp x n_trait)}
#'   \item{snp_class}{Character vector indicating the design class of each SNP}
#'   \item{V_true}{True trait covariance matrix (n_trait x n_trait)}
#'   \item{causal_matrix}{Causal effect matrix between traits}
#'   \item{A}{Additive relationship matrix (n_ind x n_ind)}
#'   \item{snp_info}{Data frame with SNP information (CHR, SNP, BP)}
#'   \item{config}{Input configuration parameters}
#' }
#'
#' @examples
#' # Basic usage with default configuration
#' data <- simulate_multitrait_data(n_ind = 200, n_trait = 3, n_snp = 100)
#'
#' # Custom genetic architecture with causal structure
#' causal_mat <- default_causal_matrix(3)
#' causal_mat[2, 1] <- 0.5  # Y1 -> Y2 effect
#' data <- simulate_multitrait_data(
#'   n_ind = 500,
#'   n_trait = 3,
#'   n_snp = 200,
#'   causal_matrix = causal_mat
#' )
#'
#' @export
#' @importFrom stats rbinom rnorm runif
#' @import Matrix MASS

# ============================================================================
# Internal Helper Functions
# ============================================================================

#' Compute Additive Relationship Matrix
#'
#' Calculates the additive relationship matrix from genotype data.
#' Implementation based on rrBLUP A.mat function.
#'
#' @param G Genotype matrix (n x m), already standardized
#'
#' @return Relationship matrix (n x n)
#'
#' @keywords internal
#'
compute_A_mat <- function(G) {
  m <- ncol(G)
  A <- tcrossprod(G) / m
  return(A)
}


#' Default Causal Matrix
#'
#' Creates a default causal effect matrix between traits.
#' Implements a simple DAG: trait 1 -> trait 2 -> trait 3
#'
#' @param n_trait Number of traits
#'
#' @return Causal effect matrix (n_trait x n_trait)
#'   causal_matrix[i,j] represents the causal effect of trait j on trait i
#'
#' @keywords internal
#'
#' @examples
#' default_causal_matrix(3)
#' # Returns:
#' #       trait_1 trait_2 trait_3
#' # trait_1   0.0    0.0    0.0
#' # trait_2   0.6    0.0    0.0
#' # trait_3   0.0    0.4    0.0
#'
default_causal_matrix <- function(n_trait) {
  C <- matrix(0, nrow = n_trait, ncol = n_trait)
  rownames(C) <- colnames(C) <- paste0("trait_", 1:n_trait)
  
  # Trait 1 -> Trait 2 (coefficient 0.6)
  if (n_trait >= 2) {
    C[2, 1] <- 0.6
  }
  
  # Trait 2 -> Trait 3 (coefficient 0.4)
  if (n_trait >= 3) {
    C[3, 2] <- 0.4
  }
  
  return(C)
}


#' Decompose Effects for Analysis (Not for Ground Truth)
#'
#' Mathematical decomposition of effects into independent and shared components.
#' This is an auxiliary analysis function, NOT used for generating ground truth labels.
#' Ground truth labels come from the simulation design (snp_class).
#'
#' @param beta Total effect matrix (n_snp x n_trait)
#' @param V Trait covariance matrix (n_trait x n_trait)
#'
#' @return List with:
#'   \item{beta_ind}{Independent effects (portion unique to each trait)}
#'   \item{beta_shared}{Shared effects (portion explained by covariance)}
#'
#' @keywords internal
#'
decompose_true_effects <- function(beta, V) {
  
  n_snp <- nrow(beta)
  m <- ncol(beta)
  
  beta_ind <- matrix(0, nrow = n_snp, ncol = m)
  beta_shared <- matrix(0, nrow = n_snp, ncol = m)
  
  colnames(beta_ind) <- colnames(beta_shared) <- colnames(beta)
  rownames(beta_ind) <- rownames(beta_shared) <- rownames(beta)
  
  for (i in 1:m) {
    # Compute conditional regression coefficients
    # b_{i,-i} = V_{-i}^{-1} * C_{-i,i}
    V_minus_i <- V[-i, -i, drop = FALSE]
    C_minus_i <- V[-i, i, drop = FALSE]
    
    b <- tryCatch(
      solve(V_minus_i, C_minus_i),
      error = function(e) {
        MASS::ginv(V_minus_i) %*% C_minus_i
      }
    )
    
    # Independent effect = beta_i - beta_{-i} %*% b
    beta_ind[, i] <- beta[, i] - beta[, -i, drop = FALSE] %*% b
  }
  
  # Shared effect = total - independent
  beta_shared <- beta - beta_ind
  
  return(list(
    beta_ind = beta_ind,
    beta_shared = beta_shared
  ))
}


# ============================================================================
# Main Simulation Functions
# ============================================================================

#' @export
#' @importFrom stats rbinom rnorm runif
#' @import Matrix MASS
#'
simulate_multitrait_data <- function(
    n_ind = 500,
    n_trait = 3,
    n_snp = 200,
    n_env = 1,
    maf_range = c(0.05, 0.5),
    trait_cor = 0.5,
    heritability = 0.5,
    effect_config = NULL,
    causal_matrix = NULL
) {
  # ============================================================================
  # Input Validation
  # ============================================================================
  stopifnot(
    "n_ind must be positive" = n_ind > 0,
    "n_trait must be at least 2" = n_trait >= 2,
    "n_snp must be positive" = n_snp > 0,
    "maf_range must have length 2" = length(maf_range) == 2,
    "maf_range[1] < maf_range[2]" = maf_range[1] < maf_range[2],
    "heritability must be between 0 and 1" = heritability > 0 && heritability < 1
  )
  
  # ============================================================================
  # Step 1: Generate Genotype
  # ============================================================================
  n_chr <- max(3, ceiling(n_snp / 30))
  chr_snpNum <- rep(floor(n_snp / n_chr), n_chr)
  diff <- n_snp - sum(chr_snpNum)
  if (diff != 0) {
    add_per_chr <- diff %/% n_chr
    remainder <- diff %% n_chr
    chr_snpNum <- chr_snpNum + add_per_chr
    if (remainder != 0) {
      chr_snpNum[seq_len(remainder)] <- chr_snpNum[seq_len(remainder)] + 1
    }
  }
  n_chr <- length(chr_snpNum)
  
  # Generate MAF
  maf <- runif(n_snp, maf_range[1], maf_range[2])
  
  # Generate genotype matrix (additive coding: -1, 0, 1)
  G_raw <- matrix(0L, nrow = n_ind, ncol = n_snp)
  for (j in seq_len(n_snp)) {
    G_raw[, j] <- rbinom(n_ind, size = 2, prob = maf[j]) - 1L
  }
  
  # Create SNP info data frame
  snp_info <- data.frame(
    CHR = rep(paste0("chr", seq_len(n_chr)), times = chr_snpNum),
    SNP = paste0("snp", seq_len(n_snp)),
    BP = unlist(lapply(chr_snpNum, seq_len)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # Standardize genotype matrix for analysis
  G <- scale(G_raw)
  G[is.na(G)] <- 0
  rownames(G) <- paste0("ind_", seq_len(n_ind))
  colnames(G) <- snp_info$SNP
  
  # ============================================================================
  # Step 2: Generate Trait Covariance Matrix
  # ============================================================================
  V_true <- compute_trait_covariance(n_trait = n_trait, rho = trait_cor)
  
  # ============================================================================
  # Step 3: Set Up Causal Matrix (if not provided)
  # ============================================================================
  if (is.null(causal_matrix)) {
    causal_matrix <- default_causal_matrix(n_trait)
  }
  
  # Validate causal_matrix dimensions
  if (!identical(dim(causal_matrix), c(n_trait, n_trait))) {
    stop("causal_matrix dimensions must match n_trait x n_trait")
  }
  
  # ============================================================================
  # Step 4: Generate Genetic Effects with 5-Class Classification
  # ============================================================================
  if (is.null(effect_config)) {
    effect_config <- default_effect_config(n_trait = n_trait, n_snp = n_snp)
  }
  
  effect_result <- simulate_genetic_effects(
    n_snp = n_snp,
    n_trait = n_trait,
    effect_config = effect_config,
    V = V_true
  )
  beta_true <- effect_result$beta
  snp_class <- effect_result$snp_class  # Design-based ground truth labels
  
  # ============================================================================
  # Step 5: Generate Phenotypes with Causal Structure
  # ============================================================================
  pheno_result <- simulate_phenotype(
    G = G,
    beta = beta_true,
    V = V_true,
    heritability = heritability,
    n_env = n_env,
    causal_matrix = causal_matrix
  )
  Y <- pheno_result$Y
  Y_base <- pheno_result$Y_base
  E <- pheno_result$E
  
  # ============================================================================
  # Step 6: Compute Additive Relationship Matrix
  # ============================================================================
  A <- compute_A_mat(G)
  rownames(A) <- colnames(A) <- rownames(G)
  
  # ============================================================================
  # Return Results
  # ============================================================================
  result <- list(
    G = G,
    G_raw = G_raw,
    Y = Y,
    Y_base = Y_base,
    E = E,
    beta_true = beta_true,
    snp_class = snp_class,  # Design-based ground truth labels
    V_true = V_true,
    causal_matrix = causal_matrix,
    A = A,
    snp_info = snp_info,
    config = list(
      n_ind = n_ind,
      n_trait = n_trait,
      n_snp = n_snp,
      n_env = n_env,
      maf_range = maf_range,
      trait_cor = trait_cor,
      heritability = heritability,
      effect_config = effect_config
    )
  )
  
  class(result) <- "condped_simdata"
  return(result)
}


#' Default Effect Configuration for 5-Class SNP Architecture
#'
#' Creates a default genetic architecture with 5 SNP classes:
#' \enumerate{
#'   \item covariance_induced: Effects only on trait 1, apparent association on trait 2 via covariance
#'   \item trait_specific: Effects only on trait 2
#'   \item horizontal: Independent effects on multiple traits (traits 1 and 3)
#'   \item vertical: Effects only on upstream trait 1, propagating via causal_matrix to trait 2
#'   \item bidirectional: Effects on multiple traits
#' }
#'
#' Each class contains 5 SNPs by default (25 functional SNPs total).
#'
#' @param n_trait Number of traits (minimum 2)
#' @param n_snp Number of SNPs
#'
#' @return List with 5 effect configurations
#'
#' @keywords internal
#'
#' @examples
#' cfg <- default_effect_config(n_trait = 3, n_snp = 200)
#' names(cfg)
#' # [1] "covariance_induced" "trait_specific" "horizontal" 
#' # [4] "vertical" "bidirectional"
#'
default_effect_config <- function(n_trait, n_snp) {
  
  n_per_class <- 5
  
  # Ensure we have enough SNPs
  total_functional <- 5 * n_per_class
  if (total_functional > n_snp) {
    n_per_class <- floor(n_snp / 5)
    total_functional <- 5 * n_per_class
  }
  
  list(
    # -------------------------------------------------------------------------
    # 1. Covariance-induced pleiotropy
    # Real genetic effect only on trait 1
    # Apparent association on trait 2 appears through V correlation
    # Condition on trait 1, effect on trait 2 disappears
    # -------------------------------------------------------------------------
    covariance_induced = list(
      snp_ids = 1:n_per_class,
      target_trait = 1,
      effect_size = 0.5
    ),
    
    # -------------------------------------------------------------------------
    # 2. Trait-specific effect
    # Real genetic effect only on trait 2
    # Assumes trait 2 is conditionally independent (or weak V correlation)
    # -------------------------------------------------------------------------
    trait_specific = list(
      snp_ids = (n_per_class + 1):(2 * n_per_class),
      target_trait = 2,
      effect_size = 0.5
    ),
    
    # -------------------------------------------------------------------------
    # 3. Horizontal pleiotropy
    # Independent real effects on multiple traits (traits 1 and 3)
    # SNP directly affects traits without mediation
    # -------------------------------------------------------------------------
    horizontal = list(
      snp_ids = (2 * n_per_class + 1):(3 * n_per_class),
      effect_matrix = matrix(
        c(0.5, 0, 0.4),  # Effects on trait 1 and 3, zero on trait 2
        nrow = n_per_class,
        ncol = n_trait,
        byrow = TRUE
      )
    ),
    
    # -------------------------------------------------------------------------
    # 4. Vertical pleiotropy
    # Real genetic effect only on upstream trait 1
    # Effect on trait 2 is entirely mediated through causal pathway Y1 -> Y2
    # IMPORTANT: beta must be strictly 0 on downstream traits
    # -------------------------------------------------------------------------
    vertical = list(
      snp_ids = (3 * n_per_class + 1):(4 * n_per_class),
      target_trait = 1,  # Only upstream trait
      effect_size = 0.5
      # Note: beta on other traits must be 0 by design
    ),
    
    # -------------------------------------------------------------------------
    # 5. Bidirectional / Confounded
    # Real effects on multiple traits with potential feedback
    # Used to test fallback classification
    # -------------------------------------------------------------------------
    bidirectional = list(
      snp_ids = (4 * n_per_class + 1):(5 * n_per_class),
      effect_matrix = matrix(
        c(0.3, 0.3, 0.3),
        nrow = n_per_class,
        ncol = n_trait,
        byrow = TRUE
      )
    )
  )
}


#' Simulate Genetic Effects with 5-Class SNP Classification
#'
#' Generates true genetic effect matrices according to the 5-class architecture.
#' Returns snp_class vector for design-based ground truth labeling.
#'
#' @param n_snp Number of SNPs
#' @param n_trait Number of traits
#' @param effect_config List with 5-class effect configuration
#' @param V Trait covariance matrix
#'
#' @return List with:
#'   \item{beta}{Total genetic effect matrix (n_snp x n_trait)}
#'   \item{snp_class}{Character vector: class label for each SNP}
#'   \item{config}{Input effect configuration}
#'
#' @keywords internal
#'
simulate_genetic_effects <- function(n_snp, n_trait, effect_config, V) {
  
  beta <- matrix(0, nrow = n_snp, ncol = n_trait)
  rownames(beta) <- paste0("snp_", 1:n_snp)
  colnames(beta) <- paste0("trait_", 1:n_trait)
  
  # Record design class for each SNP (ground truth labels)
  snp_class <- rep("null", n_snp)
  
  # ============================================================================
  # 1. Covariance-induced pleiotropy
  # ============================================================================
  if (!is.null(effect_config$covariance_induced)) {
    cfg <- effect_config$covariance_induced
    s <- cfg$snp_ids
    t <- cfg$target_trait
    
    # Validate SNP indices
    if (any(s > n_snp)) {
      warning("covariance_induced: some snp_ids exceed n_snp, skipping")
    } else {
      beta[s, t] <- cfg$effect_size
      snp_class[s] <- "covariance_induced"
    }
  }
  
  # ============================================================================
  # 2. Trait-specific effect
  # ============================================================================
  if (!is.null(effect_config$trait_specific)) {
    cfg <- effect_config$trait_specific
    s <- cfg$snp_ids
    t <- cfg$target_trait
    
    if (any(s > n_snp)) {
      warning("trait_specific: some snp_ids exceed n_snp, skipping")
    } else {
      beta[s, t] <- cfg$effect_size
      snp_class[s] <- "trait_specific"
    }
  }
  
  # ============================================================================
  # 3. Horizontal pleiotropy
  # ============================================================================
  if (!is.null(effect_config$horizontal)) {
    cfg <- effect_config$horizontal
    s <- cfg$snp_ids
    
    if (any(s > n_snp)) {
      warning("horizontal: some snp_ids exceed n_snp, skipping")
    } else {
      if (!is.null(cfg$effect_matrix)) {
        beta[s, ] <- cfg$effect_matrix
      } else if (!is.null(cfg$effect_size)) {
        beta[s, ] <- matrix(
          rnorm(length(s) * n_trait, mean = 0, sd = cfg$effect_size),
          nrow = length(s),
          ncol = n_trait
        )
      }
      snp_class[s] <- "horizontal"
    }
  }
  
  # ============================================================================
  # 4. Vertical pleiotropy
  # CRITICAL: beta must be strictly 0 on downstream traits
  # This is enforced by only setting beta on target_trait
  # ============================================================================
  if (!is.null(effect_config$vertical)) {
    cfg <- effect_config$vertical
    s <- cfg$snp_ids
    t <- cfg$target_trait
    
    if (any(s > n_snp)) {
      warning("vertical: some snp_ids exceed n_snp, skipping")
    } else {
      beta[s, t] <- cfg$effect_size
      # Key invariant: beta[s, -t] must be 0
      # This is ensured by only modifying beta[s, t]
      snp_class[s] <- "vertical"
    }
  }
  
  # ============================================================================
  # 5. Bidirectional / Confounded
  # ============================================================================
  if (!is.null(effect_config$bidirectional)) {
    cfg <- effect_config$bidirectional
    s <- cfg$snp_ids
    
    if (any(s > n_snp)) {
      warning("bidirectional: some snp_ids exceed n_snp, skipping")
    } else {
      if (!is.null(cfg$effect_matrix)) {
        beta[s, ] <- cfg$effect_matrix
      } else if (!is.null(cfg$effect_size)) {
        beta[s, ] <- matrix(
          rnorm(length(s) * n_trait, mean = 0, sd = cfg$effect_size),
          nrow = length(s),
          ncol = n_trait
        )
      }
      snp_class[s] <- "bidirectional"
    }
  }
  
  return(list(
    beta = beta,
    snp_class = snp_class,
    config = effect_config
  ))
}


#' Generate Phenotypes with Causal Pathway Structure
#'
#' Computes phenotypes as:
#' Y_base = G * beta + E
#' Y_i = Y_base_i + sum_j(causal_matrix[i,j] * Y_j)
#'
#' The causal_matrix implements cascading phenotype generation where
#' downstream traits are affected by upstream traits.
#'
#' @param G Genotype matrix (n_ind x n_snp)
#' @param beta True effect matrix (n_snp x n_trait)
#' @param V Trait covariance matrix (n_trait x n_trait)
#' @param heritability Heritability value
#' @param n_env Number of environments
#' @param causal_matrix Causal effect matrix (n_trait x n_trait).
#'   causal_matrix[i,j] represents causal effect of trait j on trait i.
#'   NULL for no causal propagation.
#'
#' @return List with:
#'   \item{Y}{Final phenotype matrix (n_ind x n_trait)}
#'   \item{Y_base}{Base phenotype without causal propagation (n_ind x n_trait)}
#'   \item{E}{Environmental effects (n_ind x n_trait)}
#'
#' @export
#' @import MASS
#'
simulate_phenotype <- function(G, beta, V, heritability, n_env = 1,
                                causal_matrix = NULL) {
  
  n_ind <- nrow(G)
  n_trait <- ncol(V)
  
  # ============================================================================
  # Genetic component: G %*% beta
  # ============================================================================
  G_effect <- G %*% beta
  
  # ============================================================================
  # Environmental component
  # Scale V to match desired heritability
  # ============================================================================
  var_g <- apply(G_effect, 2, var)
  
  # Environmental variance to achieve target heritability
  ve_factor <- (1 - heritability) / heritability
  Sigma_e <- ve_factor * outer(sqrt(var_g), sqrt(var_g)) * V
  
  # Add small jitter to ensure positive definiteness
  jitter <- 1e-6
  Sigma_e <- Sigma_e + jitter * diag(n_trait)
  
  # Generate environmental effects
  E <- MASS::mvrnorm(n = n_ind, mu = rep(0, n_trait), Sigma = Sigma_e)
  
  # ============================================================================
  # Base phenotype (genetic + environmental)
  # ============================================================================
  Y_base <- G_effect + E
  
  # ============================================================================
  # Apply causal propagation if causal_matrix is provided
  # Process in topological order (assumes upstream traits have smaller indices)
  # ============================================================================
  Y <- Y_base
  
  if (!is.null(causal_matrix)) {
    # Simple topological processing: iterate from trait 2 to n_trait
    # More robust: use DAG topological sort
    processed <- rep(FALSE, n_trait)
    
    # Find topological order (simple greedy approach)
    # A trait can be processed when all its parents have been processed
    for (iteration in 1:n_trait) {
      for (i in 1:n_trait) {
        if (!processed[i]) {
          parents <- which(causal_matrix[i, ] != 0)
          if (length(parents) == 0 || all(processed[parents])) {
            # All parents processed, now we can add effects from parents
            if (length(parents) > 0) {
              for (p in parents) {
                Y[, i] <- Y[, i] + causal_matrix[i, p] * Y[, p]
              }
            }
            processed[i] <- TRUE
          }
        }
      }
    }
  }
  
  # Set names
  rownames(Y) <- rownames(Y_base) <- rownames(G)
  colnames(Y) <- colnames(Y_base) <- paste0("trait_", 1:n_trait)
  rownames(E) <- rownames(G)
  colnames(E) <- paste0("trait_", 1:n_trait)
  
  return(list(Y = Y, Y_base = Y_base, E = E))
}


#' Compute Trait Covariance Matrix
#'
#' Creates an AR(1) structured covariance matrix.
#'
#' @param n_trait Number of traits
#' @param rho Correlation parameter
#'
#' @return Covariance matrix (n_trait x n_trait)
#'
#' @keywords internal
#'
compute_trait_covariance <- function(n_trait, rho = 0.5) {
  # AR(1) structure
  V <- rho^abs(outer(1:n_trait, 1:n_trait, "-"))
  
  # Scale to have variance 1
  V <- V / sqrt(diag(V) %o% diag(V))
  
  rownames(V) <- colnames(V) <- paste0("trait_", 1:n_trait)
  return(V)
}


# ============================================================================
# Print and Summary Methods
# ============================================================================

#' Print Method for Simulated Data
#'
#' @param x Object of class condped_simdata
#' @param ... Additional arguments passed to print
#'
#' @export
#'
print.condped_simdata <- function(x, ...) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("         CondPED Simulated Data\n")
  cat("═══════════════════════════════════════════════════════════\n")
  
  cat(sprintf("  Individuals:      %d\n", nrow(x$Y)))
  cat(sprintf("  Traits:           %d\n", ncol(x$Y)))
  cat(sprintf("  SNPs:             %d\n", ncol(x$G)))
  cat(sprintf("  Heritability:     %.2f\n", x$config$heritability))
  
  cat("\n  ─────────────────────────────────────────────────────────\n")
  cat("  SNP Classification (Ground Truth):\n")
  cat("  ─────────────────────────────────────────────────────────\n")
  
  # Count SNPs by class
  snp_class <- x$snp_class
  class_counts <- table(snp_class)
  for (nm in names(class_counts)) {
    cat(sprintf("    %-20s: %d SNPs\n", nm, class_counts[nm]))
  }
  
  cat("\n  ─────────────────────────────────────────────────────────\n")
  cat("  Causal Structure (causal_matrix):\n")
  cat("  ─────────────────────────────────────────────────────────\n")
  if (!is.null(x$causal_matrix)) {
    print(round(x$causal_matrix, 3))
  } else {
    cat("    (No causal structure)\n")
  }
  
  cat("\n  ─────────────────────────────────────────────────────────\n")
  cat("  Trait Correlations (V_true):\n")
  cat("  ─────────────────────────────────────────────────────────\n")
  print(round(x$V_true, 3))
  
  cat("\n  ─────────────────────────────────────────────────────────\n")
  cat("  Sample Phenotype Correlations:\n")
  cat("  ─────────────────────────────────────────────────────────\n")
  print(round(cor(x$Y), 3))
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("\n")
  
  invisible(x)
}


#' Summary Method for Simulated Data
#'
#' @param object Object of class condped_simdata
#' @param ... Additional arguments
#'
#' @export
#'
summary.condped_simdata <- function(object, ...) {
  
  cat("\n")
  cat("CondPED Simulated Data Summary\n")
  cat("==============================\n\n")
  
  # Dimensions
  cat("Dimensions:\n")
  cat(sprintf("  - Individuals: %d\n", nrow(object$Y)))
  cat(sprintf("  - Traits: %d\n", ncol(object$Y)))
  cat(sprintf("  - SNPs: %d\n", ncol(object$G)))
  
  # SNP classification summary
  cat("\nSNP Classification (Ground Truth):\n")
  snp_class <- object$snp_class
  class_counts <- table(snp_class)
  for (nm in names(class_counts)) {
    cat(sprintf("  - %-20s: %d SNPs\n", nm, class_counts[nm]))
  }
  
  # Effect summary
  cat("\nTrue Effect Summary:\n")
  cat(sprintf("  - Total genetic effects: %d non-zero\n", sum(object$beta_true != 0)))
  
  # Effects by class
  cat("\nEffects by SNP Class:\n")
  for (cls in names(table(snp_class))) {
    idx <- which(snp_class == cls)
    if (length(idx) > 0) {
      non_zero <- sum(object$beta_true[idx, ] != 0)
      cat(sprintf("  - %-20s: %d non-zero effects\n", cls, non_zero))
    }
  }
  
  # Variance components
  var_g <- apply(object$G %*% object$beta_true, 2, var)
  var_e <- apply(object$E, 2, var)
  
  cat("\nVariance Components:\n")
  cat("  Trait    Var_G    Var_E    h2 (sample)\n")
  cat("  ----------------------------------------\n")
  for (i in 1:ncol(object$Y)) {
    h2_sample <- var_g[i] / (var_g[i] + var_e[i])
    cat(sprintf("  %-7s %7.3f %7.3f %7.3f\n",
                colnames(object$Y)[i], var_g[i], var_e[i], h2_sample))
  }
  
  # Causal structure
  cat("\nCausal Matrix Summary:\n")
  if (!is.null(object$causal_matrix)) {
    causal_elements <- which(object$causal_matrix != 0, arr.ind = TRUE)
    if (nrow(causal_elements) > 0) {
      for (r in 1:nrow(causal_elements)) {
        i <- causal_elements[r, 1]
        j <- causal_elements[r, 2]
        cat(sprintf("  - trait_%d <- %.2f * trait_%d\n", i, object$causal_matrix[i, j], j))
      }
    } else {
      cat("  (No causal structure)\n")
    }
  }
  
  cat("\n")
  invisible(object)
}
