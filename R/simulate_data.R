#' Simulate multi-trait phenotype data with known effect patterns
#'
#' Generates synthetic phenotype and genotype data under specified genetic 
#' architectures for simulation studies.
#'
#' @param n_ind Integer. Number of individuals (default: 500)
#' @param n_trait Integer. Number of traits (default: 3)
#' @param n_snp Integer. Number of SNPs (default: 100)
#' @param n_env Integer. Number of environments (default: 1)
#' @param maf Numeric. Minor allele frequency (default: 0.3)
#' @param V Matrix. Trait covariance matrix. If NULL, uses AR(1) with rho=0.5
#' @param beta_config List. Specifies effect patterns. See Details.
#'
#' @details
#' The `beta_config` parameter controls genetic effect patterns:
#' \itemize{
#'   \item `trait_specific`: List with `snp_ids`, `trait_ids`, `effect_size`
#'   \item `shared`: List with `snp_ids`, `pattern = "covariance_mediated"`
#'   \item `independent`: List with `snp_ids`, `effects` (matrix)
#' }
#'
#' @return List containing:
#' \describe{
#'   \item{Y}{Phenotype matrix (n_ind x n_trait)}
#'   \item{G}{Genotype matrix (n_ind x n_snp), standardized}
#'   \item{beta_true}{True effect matrix (n_snp x n_trait)}
#'   \item{beta_ind_true}{True independent effects}
#'   \item{beta_shared_true}{True shared effects}
#'   \item{V_true}{True trait covariance matrix}
#'   \item{config}{Input configuration}
#' }
#'
#' @examples
#' # Trait-specific QTL scenario
#' data <- simulate_multitrait_data(
#'   n_ind = 200,
#'   n_trait = 3,
#'   n_snp = 50,
#'   beta_config = list(
#'     trait_specific = list(
#'       snp_ids = 1:10,
#'       trait_ids = 1,
#'       effect_size = 0.3
#'     )
#'   )
#' )
#'
#' @export
simulate_multitrait_data <- function(
    n_ind = 500,
    n_trait = 3,
    n_snp = 100,
    n_env = 1,
    maf = 0.3,
    V = NULL,
    beta_config = NULL
) {
  
  # Input validation
  stopifnot(
    n_ind > 0, n_trait > 1, n_snp > 0,
    maf > 0, maf < 0.5
  )
  
  # 1. Generate genotype matrix (Hardy-Weinberg equilibrium)
  p <- maf
  G <- matrix(
    rbinom(n_ind * n_snp, size = 2, prob = p),
    nrow = n_ind, 
    ncol = n_snp
  )
  
  # Standardize genotypes
  G <- scale(G)
  colnames(G) <- paste0("SNP", 1:n_snp)
  
  # 2. Construct trait covariance matrix
  if (is.null(V)) {
    # Default: AR(1) structure with rho = 0.5
    rho <- 0.5
    V <- rho^abs(outer(1:n_trait, 1:n_trait, "-"))
    rownames(V) <- colnames(V) <- paste0("Trait", 1:n_trait)
  }
  
  # 3. Initialize effect matrix
  beta_true <- matrix(0, nrow = n_snp, ncol = n_trait)
  rownames(beta_true) <- paste0("SNP", 1:n_snp)
  colnames(beta_true) <- paste0("Trait", 1:n_trait)
  
  # 4. Configure effects based on scenarios
  if (is.null(beta_config)) {
    # Default: mixed scenario
    beta_config <- list(
      trait_specific = list(
        snp_ids = 1:10,
        trait_ids = 1,
        effect_size = 0.3
      ),
      shared = list(
        snp_ids = 11:20,
        pattern = "covariance_mediated"
      ),
      independent = list(
        snp_ids = 21:30,
        effects = matrix(rnorm(10 * n_trait, mean = 0, sd = 0.2), 
                        ncol = n_trait)
      )
    )
  }
  
  # 4.1 Trait-specific QTL
  if (!is.null(beta_config$trait_specific)) {
    cfg <- beta_config$trait_specific
    beta_true[cfg$snp_ids, cfg$trait_ids] <- cfg$effect_size
  }
  
  # 4.2 Shared QTL (covariance-mediated pattern)
  if (!is.null(beta_config$shared)) {
    cfg <- beta_config$shared
    if (cfg$pattern == "covariance_mediated") {
      # Set effect on trait 1
      beta_true[cfg$snp_ids, 1] <- 0.3
      
      # Effects on other traits = V[i,1]/V[1,1] * beta[1]
      for (i in 2:n_trait) {
        beta_true[cfg$snp_ids, i] <- (V[i, 1] / V[1, 1]) * 
                                       beta_true[cfg$snp_ids, 1]
      }
    }
  }
  
  # 4.3 Independent pleiotropy
  if (!is.null(beta_config$independent)) {
    cfg <- beta_config$independent
    beta_true[cfg$snp_ids, ] <- cfg$effects
  }
  
  # 5. Generate phenotype matrix
  Y <- matrix(0, nrow = n_ind, ncol = n_trait)
  colnames(Y) <- paste0("Trait", 1:n_trait)
  
  # Genetic component
  G_effect <- G %*% beta_true
  
  # Environmental error (multivariate normal)
  E <- MASS::mvrnorm(n = n_ind, mu = rep(0, n_trait), Sigma = V)
  
  # Total phenotype
  Y <- G_effect + E
  
  # 6. Decompose true effects into independent and shared components
  beta_decomp <- decompose_true_effects(beta_true, V)
  
  # Return
  structure(
    list(
      Y = Y,
      G = G,
      beta_true = beta_true,
      beta_ind_true = beta_decomp$beta_ind,
      beta_shared_true = beta_decomp$beta_shared,
      V_true = V,
      config = beta_config
    ),
    class = "condped_simdata"
  )
}


#' Decompose true effects into independent and shared components
#'
#' @param beta Effect matrix (n_snp x n_trait)
#' @param V Trait covariance matrix (n_trait x n_trait)
#'
#' @return List with beta_ind and beta_shared matrices
#' @keywords internal
decompose_true_effects <- function(beta, V) {
  
  n_snp <- nrow(beta)
  m <- ncol(beta)
  
  beta_ind <- matrix(0, nrow = n_snp, ncol = m)
  colnames(beta_ind) <- colnames(beta)
  rownames(beta_ind) <- rownames(beta)
  
  for (i in 1:m) {
    # Compute conditional regression coefficients b_{i,-i}
    V_minus_i <- V[-i, -i, drop = FALSE]
    C_minus_i <- V[-i, i]
    
    b <- solve(V_minus_i, C_minus_i)
    
    # Independent effect = beta_i - beta_{-i} %*% b
    beta_ind[, i] <- beta[, i] - beta[, -i, drop = FALSE] %*% b
  }
  
  beta_shared <- beta - beta_ind
  
  list(
    beta_ind = beta_ind,
    beta_shared = beta_shared
  )
}


#' Print method for simulated data
#' @export
print.condped_simdata <- function(x, ...) {
  cat("CondPED Simulated Data\n")
  cat("======================\n")
  cat(sprintf("Individuals: %d\n", nrow(x$Y)))
  cat(sprintf("Traits: %d\n", ncol(x$Y)))
  cat(sprintf("SNPs: %d\n", ncol(x$G)))
  cat(sprintf("Non-zero effects: %d\n", sum(x$beta_true != 0)))
  cat("\nTrait correlations:\n")
  print(round(cor(x$Y), 3))
  invisible(x)
}