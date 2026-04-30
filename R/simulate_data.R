#' Simulate Multi-Trait Phenotype and Genotype Data
#'
#' Generates synthetic genotype and phenotype data with known genetic architectures
#' for testing the conditional projection framework.
#'
#' This is the main user-facing function that combines:
#' \itemize{
#'   \item \code{\link{simulate_genotype}} - Generate genotype matrix
#'   \item \code{\link{simulate_genetic_effects}} - Generate true genetic effects
#'   \item \code{\link{simulate_phenotype}} - Generate phenotypes from genotype and effects
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
#'
#' @details
#' The \code{effect_config} parameter controls the genetic architecture:
#' \describe{
#'   \item{trait_specific}{List with elements:
#'     \describe{
#'       \item{snp_ids}{Integer vector of SNP indices}
#'       \item{trait_ids}{Integer vector of trait indices}
#'       \item{effect_size}{Numeric, effect size (standardized)}
#'     }
#'   }
#'   \item{shared}{List with elements:
#'     \describe{
#'       \item{snp_ids}{Integer vector of SNP indices}
#'       \item{effect_size}{Numeric, effect size on the primary trait}
#'     }
#'   }
#'   \item{independent}{List with elements:
#'     \describe{
#'       \item{snp_ids}{Integer vector of SNP indices}
#'       \item{effect_matrix}{Matrix of effect sizes (n_snp x n_trait)}
#'     }
#'   }
#' }
#'
#' @return List containing:
#' \describe{
#'   \item{G}{Genotype matrix (n_ind x n_snp), standardized {-1, 0, 1}}
#'   \item{Y}{Phenotype matrix (n_ind x n_trait)}
#'   \item{E}{Environmental residuals (n_ind x n_trait)}
#'   \item{beta_true}{True total effect matrix (n_snp x n_trait)}
#'   \item{beta_ind_true}{True independent effects (n_snp x n_trait)}
#'   \item{beta_shared_true}{True shared effects (n_snp x n_trait)}
#'   \item{V_true}{True trait covariance matrix (n_trait x n_trait)}
#'   \item{A}{Additive relationship matrix (n_ind x n_ind)}
#'   \item{snp_info}{Data frame with SNP information (CHR, SNP, BP)}
#'   \item{config}{Input configuration parameters}
#' }
#'
#' @examples
#' # Basic usage with default configuration
#' data <- simulate_multitrait_data(n_ind = 200, n_trait = 3, n_snp = 100)
#'
#' # Custom genetic architecture
#' data <- simulate_multitrait_data(
#'   n_ind = 500,
#'   n_trait = 3,
#'   n_snp = 200,
#'   effect_config = list(
#'     trait_specific = list(
#'       snp_ids = 1:10,
#'       trait_ids = 1,
#'       effect_size = 0.5
#'     ),
#'     shared = list(
#'       snp_ids = 11:20,
#'       effect_size = 0.3
#'     ),
#'     independent = list(
#'       snp_ids = 21:30,
#'       effect_matrix = matrix(rnorm(10 * 3, 0, 0.2), nrow = 10, ncol = 3)
#'     )
#'   )
#' )
#'
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
  # Relationship matrix: A = G %*% t(G) / m
  # G is already standardized (mean=0, var=1)
  m <- ncol(G)
  A <- tcrossprod(G) / m

  return(A)
}

#' @export
#' @importFrom stats rbinom rnorm runif cor
#' @import Matrix
#'
simulate_multitrait_data <- function(
    n_ind = 500,
    n_trait = 3,
    n_snp = 200,
    n_env = 1,
    maf_range = c(0.05, 0.5),
    trait_cor = 0.5,
    heritability = 0.5,
    effect_config = NULL
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
  geno_result <- simulate_genotype(
    n_ind = n_ind,
    n_snp = n_snp,
    maf_range = maf_range
  )
  G <- geno_result$G          # Standardized genotype (-1, 0, 1 after scaling)
  G_raw <- geno_result$G_raw  # Original genotype (-1, 0, 1 before scaling)
  snp_info <- geno_result$snp_info

  # ============================================================================
  # Step 2: Generate Trait Covariance Matrix
  # ============================================================================
  V_true <- compute_trait_covariance(n_trait = n_trait, rho = trait_cor)

  # ============================================================================
  # Step 3: Generate Genetic Effects
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
  beta_ind_true <- effect_result$beta_ind
  beta_shared_true <- effect_result$beta_shared

  # ============================================================================
  # Step 4: Generate Phenotypes
  # ============================================================================
  pheno_result <- simulate_phenotype(
    G = G,
    beta = beta_true,
    V = V_true,
    heritability = heritability,
    n_env = n_env
  )
  Y <- pheno_result$Y
  E <- pheno_result$E

  # ============================================================================
  # Step 5: Compute Additive Relationship Matrix
  # ============================================================================
  A <- compute_A_mat(G)
  rownames(A) <- colnames(A) <- rownames(G)

  # ============================================================================
  # Return Results
  # ============================================================================
  result <- list(
    G = G,
    G_raw = G_raw,  # Original genotype for QTLNetwork
    Y = Y,
    E = E,
    beta_true = beta_true,
    beta_ind_true = beta_ind_true,
    beta_shared_true = beta_shared_true,
    V_true = V_true,
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


#' Generate Genotype Matrix
#'
#' Simulates genotype data with Hardy-Weinberg equilibrium.
#' Genotypes are coded as -1, 0, 1 (additive coding).
#'
#' @param n_ind Number of individuals
#' @param n_snp Number of SNPs
#' @param maf_range Range for minor allele frequency
#' @param chr_snpNum Vector. Number of SNPs per chromosome (optional)
#'
#' @return List with G (standardized genotype), G_raw (original -1,0,1), and snp_info
#'
#' @export
#' @importFrom stats runif rbinom
#'
simulate_genotype <- function(n_ind, n_snp, maf_range = c(0.05, 0.5), chr_snpNum = NULL) {

  # Generate MAF for each SNP
  maf <- runif(n_snp, maf_range[1], maf_range[2])

  # Generate genotype matrix (additive coding: -1, 0, 1)
  G_raw <- matrix(0, nrow = n_ind, ncol = n_snp)
  for (j in 1:n_snp) {
    # rbinom with size=2 gives 0,1,2; subtract 1 to get -1,0,1
    G_raw[, j] <- rbinom(n_ind, size = 2, prob = maf[j]) - 1
  }

  # Save original genotype before standardization
  Ga <- G_raw

  # Standardize genotype matrix
  Ga <- scale(Ga)
  Ga[is.na(Ga)] <- 0  # Handle constant columns

  # Set row and column names
  rownames(Ga) <- rownames(G_raw) <- paste0("ind_", 1:n_ind)
  colnames(Ga) <- colnames(G_raw) <- paste0("snp_", 1:n_snp)

  # Create SNP info
  if (is.null(chr_snpNum)) {
    # Default: equal distribution across chromosomes
    n_chr <- max(5, floor(n_snp / 20))
    chr_snpNum <- rep(floor(n_snp / n_chr), n_chr)

    # Adjust to match exact n_snp
    diff <- n_snp - sum(chr_snpNum)
    if (diff != 0) {
      # Distribute the difference across first chromosomes
      for (i in 1:abs(diff)) {
        chr_snpNum[i] <- chr_snpNum[i] + sign(diff)
      }
    }
  }

  chr_names <- paste0("chr", 1:length(chr_snpNum))
  snp_info <- data.frame(
    CHR = rep(1:length(chr_snpNum), times = chr_snpNum),
    SNP = paste0("snp_", 1:n_snp),
    BP = unlist(lapply(chr_snpNum, function(x) 1:x)),
    MAF = maf,
    stringsAsFactors = FALSE
  )

  return(list(G = Ga, G_raw = G_raw, snp_info = snp_info))
}


#' Generate Genetic Effects (Additive + Epistasis)
#'
#' Creates true genetic effect matrices with specified architecture.
#' Includes both additive effects and additive-by-additive epistatic effects.
#' Decomposes total effects into independent and shared components.
#'
#' @param n_snp Number of SNPs
#' @param n_trait Number of traits
#' @param effect_config List specifying effect patterns (additive and epistatic)
#' @param V Trait covariance matrix
#'
#' @return List with:
#'   \item{beta_add}{Additive effect matrix (n_snp x n_trait)}
#'   \item{beta_epi}{Epistatic effect matrix (n_epairs x n_trait)}
#'   \item{epairs}{Matrix of SNP pairs for epistasis (n_epairs x 2)}
#'   \item{beta_add_ind}{Additive independent effects}
#'   \item{beta_add_shared}{Additive shared effects}
#'   \item{beta_epi_ind}{Epistatic independent effects}
#'   \item{beta_epi_shared}{Epistatic shared effects}
#'
#' @export
#'
simulate_genetic_effects <- function(n_snp, n_trait, effect_config, V) {

  # Initialize additive effect matrix
  beta_add <- matrix(0, nrow = n_snp, ncol = n_trait)
  rownames(beta_add) <- paste0("snp_", 1:n_snp)
  colnames(beta_add) <- paste0("trait_", 1:n_trait)

  # ============================================================================
  # Additive effects
  # ============================================================================

  # ---- Trait-specific additive effects ----
  if (!is.null(effect_config$add_trait_specific)) {
    cfg <- effect_config$add_trait_specific
    snp_ids <- cfg$snp_ids
    trait_ids <- cfg$trait_ids

    if (any(snp_ids > n_snp)) {
      stop("snp_ids in add_trait_specific exceeds n_snp")
    }

    for (i in seq_along(snp_ids)) {
      beta_add[snp_ids[i], trait_ids[i]] <- cfg$effect_size[i %% length(cfg$effect_size) + 1]
    }
  }

  # ---- Shared additive effects (pleiotropic via covariance) ----
  if (!is.null(effect_config$add_shared)) {
    cfg <- effect_config$add_shared
    snp_ids <- cfg$snp_ids

    if (any(snp_ids > n_snp)) {
      stop("snp_ids in add_shared exceeds n_snp")
    }

    # Set effect on first trait
    beta_add[snp_ids, 1] <- cfg$effect_size

    # Effects on other traits follow covariance-mediated pattern
    for (j in 2:n_trait) {
      beta_add[snp_ids, j] <- (V[j, 1] / V[1, 1]) * beta_add[snp_ids, 1]
    }
  }

  # ---- Independent additive pleiotropy ----
  if (!is.null(effect_config$add_independent)) {
    cfg <- effect_config$add_independent
    snp_ids <- cfg$snp_ids

    if (any(snp_ids > n_snp)) {
      stop("snp_ids in add_independent exceeds n_snp")
    }

    if (!is.null(cfg$effect_matrix)) {
      beta_add[snp_ids, ] <- cfg$effect_matrix
    } else if (!is.null(cfg$effect_size)) {
      beta_add[snp_ids, ] <- matrix(
        rnorm(length(snp_ids) * n_trait, mean = 0, sd = cfg$effect_size),
        nrow = length(snp_ids),
        ncol = n_trait
      )
    }
  }

  # ============================================================================
  # Epistatic effects (Additive x Additive)
  # ============================================================================

  beta_epi <- matrix(0, nrow = 0, ncol = n_trait)
  epairs <- matrix(nrow = 0, ncol = 2)

  # ---- Trait-specific epistatic effects ----
  if (!is.null(effect_config$epi_trait_specific)) {
    cfg <- effect_config$epi_trait_specific

    # Generate SNP pairs for epistasis
    n_epairs <- length(cfg$snp_pairs)
    if (n_epairs > 0) {
      new_epairs <- matrix(unlist(cfg$snp_pairs), nrow = n_epairs, byrow = TRUE)

      # Validate
      if (any(new_epairs > n_snp)) {
        stop("SNP indices in epi_trait_specific exceed n_snp")
      }

      # Add to epairs
      epairs <- rbind(epairs, new_epairs)

      # Create effect matrix
      epi_effects <- matrix(0, nrow = n_epairs, ncol = n_trait)
      for (i in 1:n_epairs) {
        epi_effects[i, cfg$trait_id] <- cfg$effect_size[i %% length(cfg$effect_size) + 1]
      }
      beta_epi <- rbind(beta_epi, epi_effects)
    }
  }

  # ---- Shared epistatic effects ----
  if (!is.null(effect_config$epi_shared)) {
    cfg <- effect_config$epi_shared
    n_epairs <- length(cfg$snp_pairs)

    if (n_epairs > 0) {
      new_epairs <- matrix(unlist(cfg$snp_pairs), nrow = n_epairs, byrow = TRUE)

      if (any(new_epairs > n_snp)) {
        stop("SNP indices in epi_shared exceed n_snp")
      }

      epairs <- rbind(epairs, new_epairs)

      # Shared effects across all traits
      epi_effects <- matrix(rep(cfg$effect_size, n_epairs * n_trait),
                           nrow = n_epairs, ncol = n_trait)
      beta_epi <- rbind(beta_epi, epi_effects)
    }
  }

  # ---- Independent epistatic effects ----
  if (!is.null(effect_config$epi_independent)) {
    cfg <- effect_config$epi_independent
    n_epairs <- length(cfg$snp_pairs)

    if (n_epairs > 0) {
      new_epairs <- matrix(unlist(cfg$snp_pairs), nrow = n_epairs, byrow = TRUE)

      if (any(new_epairs > n_snp)) {
        stop("SNP indices in epi_independent exceed n_snp")
      }

      epairs <- rbind(epairs, new_epairs)

      # Random independent epistatic effects
      epi_effects <- matrix(
        rnorm(n_epairs * n_trait, mean = 0, sd = cfg$effect_size),
        nrow = n_epairs,
        ncol = n_trait
      )
      beta_epi <- rbind(beta_epi, epi_effects)
    }
  }

  # Add row names to epairs
  if (nrow(epairs) > 0) {
    rownames(epairs) <- paste0("epi_", 1:nrow(epairs))
    rownames(beta_epi) <- rownames(epairs)
  }

  # ============================================================================
  # Return results (legacy compatibility: also compute total beta)
  # ============================================================================

  # For backward compatibility, compute total beta (additive only for now)
  # Epistatic effects are stored separately
  beta <- beta_add

  # Decompose additive effects
  decomp_add <- decompose_true_effects(beta_add, V)
  beta_add_ind <- decomp_add$beta_ind
  beta_add_shared <- decomp_add$beta_shared

  # Decompose epistatic effects (same V structure)
  if (nrow(beta_epi) > 0) {
    decomp_epi <- decompose_true_effects(beta_epi, V)
    beta_epi_ind <- decomp_epi$beta_ind
    beta_epi_shared <- decomp_epi$beta_shared
  } else {
    beta_epi_ind <- matrix(0, nrow = 0, ncol = n_trait)
    beta_epi_shared <- matrix(0, nrow = 0, ncol = n_trait)
  }

  return(list(
    beta_add = beta_add,
    beta_epi = beta_epi,
    epairs = epairs,
    beta_add_ind = beta_add_ind,
    beta_add_shared = beta_add_shared,
    beta_epi_ind = beta_epi_ind,
    beta_epi_shared = beta_epi_shared,
    # Legacy compatibility
    beta = beta,
    beta_ind = beta_add_ind,
    beta_shared = beta_add_shared
  ))
}


#' Decompose True Effects
#'
#' Decomposes total genetic effects into independent and shared components
#' based on the conditional projection framework.
#'
#' @param beta Total effect matrix (n_snp x n_trait)
#' @param V Trait covariance matrix (n_trait x n_trait)
#'
#' @return List with beta_ind and beta_shared matrices
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


#' Generate Phenotypes from Genotype and Effects
#'
#' Computes phenotypes as: Y = G * beta + E
#' where E ~ MVN(0, V * (1-h2)/h2 * sigma_g^2)
#'
#' @param G Genotype matrix (n_ind x n_snp)
#' @param beta True effect matrix (n_snp x n_trait)
#' @param V Trait covariance matrix
#' @param heritability Heritability value
#' @param n_env Number of environments
#'
#' @return List with Y (phenotype) and E (environmental effects)
#'
#' @export
#' @import MASS
#'
simulate_phenotype <- function(G, beta, V, heritability, n_env = 1) {

  n_ind <- nrow(G)
  n_trait <- ncol(V)

  # Genetic component: G %*% beta
  G_effect <- G %*% beta

  # Compute genetic variance for each trait
  var_g <- apply(G_effect, 2, var)

  # Environmental variance: V_e = V * (1-h2)/h2 * var_g
  # Scale V to match the desired heritability
  # Add small jitter to ensure positive definiteness
  jitter <- 1e-6
  V_scaled <- V * outer(sqrt(var_g), sqrt(var_g))
  diag(V_scaled) <- var_g + jitter

  # Ensure positive definite
  V_scaled <- V_scaled + jitter * diag(nrow(V_scaled))

  # Environmental variance to achieve target heritability
  ve_factor <- (1 - heritability) / heritability
  Sigma_e <- ve_factor * V_scaled

  # Generate environmental effects
  E <- mvrnorm(n = n_ind, mu = rep(0, n_trait), Sigma = Sigma_e)

  # Total phenotype
  Y <- G_effect + E

  # Set names
  rownames(Y) <- rownames(G)
  colnames(Y) <- paste0("trait_", 1:n_trait)
  rownames(E) <- rownames(G)
  colnames(E) <- paste0("trait_", 1:n_trait)

  return(list(Y = Y, E = E))
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


#' Default Effect Configuration
#'
#' Creates a default genetic architecture for testing.
#' Includes both additive and epistatic effects.
#'
#' @param n_trait Number of traits
#' @param n_snp Number of SNPs
#'
#' @return List with default effect configuration
#'
#' @keywords internal
#'
default_effect_config <- function(n_trait, n_snp) {

  # ---- Additive Effects ----

  # Trait-specific additive: first 5 SNPs, each affects one trait
  add_trait_specific <- list(
    snp_ids = 1:min(5, floor(n_snp * 0.05)),
    trait_ids = rep(1:n_trait, length.out = min(5, floor(n_snp * 0.05))),
    effect_size = rep(0.5, min(5, floor(n_snp * 0.05)))
  )

  # Shared additive: next 5 SNPs, affect all traits via covariance
  n_shared <- min(5, floor(n_snp * 0.05))
  add_shared <- list(
    snp_ids = (max(add_trait_specific$snp_ids) + 1):(max(add_trait_specific$snp_ids) + n_shared),
    effect_size = 0.3
  )

  # Independent additive: next 5 SNPs, independent effects on multiple traits
  n_indep <- min(5, floor(n_snp * 0.05))
  add_independent <- list(
    snp_ids = (max(add_shared$snp_ids) + 1):(max(add_shared$snp_ids) + n_indep),
    effect_size = 0.2
  )

  # ---- Epistatic Effects ----

  # Generate some SNP pairs for epistasis (if enough SNPs available)
  if (n_snp >= 10) {
    # Trait-specific epistasis: one pair affecting only one trait
    epi_trait_specific <- list(
      snp_pairs = list(c(n_snp - 4, n_snp - 3)),
      trait_id = 1,
      effect_size = 0.3
    )

    # Shared epistasis: one pair affecting all traits
    epi_shared <- list(
      snp_pairs = list(c(n_snp - 2, n_snp - 1)),
      effect_size = 0.2
    )
  } else {
    epi_trait_specific <- NULL
    epi_shared <- NULL
  }

  return(list(
    add_trait_specific = add_trait_specific,
    add_shared = add_shared,
    add_independent = add_independent,
    epi_trait_specific = epi_trait_specific,
    epi_shared = epi_shared
  ))
}


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
  cat("  Genetic Architecture:\n")
  cat("  ─────────────────────────────────────────────────────────\n")

  # Count non-zero effects by category
  n_trait_specific <- sum(x$beta_true != 0 &
                          outer(1:nrow(x$beta_true), 1:ncol(x$beta_true),
                                function(i, j) {
                                  decomp <- decompose_true_effects(x$beta_true, x$V_true)
                                  i <= length(x$config$effect_config$trait_specific$snp_ids) &
                                  abs(decomp$beta_shared[i, j]) < 1e-10
                                }))

  cat(sprintf("  Non-zero effects: %d\n", sum(x$beta_true != 0)))

  cat("\n  ─────────────────────────────────────────────────────────\n")
  cat("  Trait Correlations (Target):\n")
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

  # Effect summary
  cat("\nTrue Effect Summary:\n")
  cat(sprintf("  - Total genetic effects: %d non-zero\n", sum(object$beta_true != 0)))
  cat(sprintf("  - Independent effects: %d non-zero\n", sum(object$beta_ind_true != 0)))
  cat(sprintf("  - Shared effects: %d non-zero\n", sum(object$beta_shared_true != 0)))

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

  cat("\n")
  invisible(object)
}
