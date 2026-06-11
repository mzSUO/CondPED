# ==============================================================================
# CondPED Simulation System (RIL + Additive + Optional Epistasis)
# ==============================================================================
# Changes:
#   1. sim_genotype: supports RIL (-1/1) and F2 (-1/0/1) via population parameter
#   2. sim_phenotype_M1: supports additive + additive epistasis (AA)
#   3. All generate_* functions: pass population and epistasis parameters
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Basic Helper Functions
# ------------------------------------------------------------------------------

#' Simulate SNP Genotype Matrix
#'
#' Supports RIL (homozygous, -1/1) and F2 (with heterozygotes, -1/0/1).
#'
#' @param n Sample size.
#' @param p Number of SNPs.
#' @param maf Minor allele frequency. Default 0.3.
#' @param population Population type: "RIL" (default, homozygous -1/1) or "F2" (-1/0/1).
#' @param seed Random seed.
#' @return Integer matrix, n x p. RIL: {-1, 1}; F2: {-1, 0, 1}.
#' @export
sim_genotype <- function(n, p, maf = 0.3, population = c("RIL", "F2"), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(n > 0, p > 0, maf > 0, maf < 1)
  population <- match.arg(population)

  if (population == "RIL") {
    # RIL: homozygous only. Sample -1 or 1 with prob (1-maf) and maf
    X <- matrix(
      sample(c(-1, 1), n * p, replace = TRUE, prob = c(1 - maf, maf)),
      nrow = n, ncol = p
    )
  } else {
    # F2: -1/0/1 dosage coding
    X <- matrix(stats::rbinom(n * p, size = 2, prob = maf), nrow = n, ncol = p)
    X[X == 0] <- -1
    X[X == 1] <- 0
    X[X == 2] <- 1
  }

  colnames(X) <- paste0("SNP", seq_len(p))
  rownames(X) <- paste0("Ind", seq_len(n))
  X
}


#' Calculate Effect Size SD from PVE
#'
#' Formula: sd = sqrt(PVE / Var(X)), where Var(X) = 2 * maf * (1 - maf) for F2,
#' and Var(X) = 4 * maf * (1 - maf) for RIL (since X is -1/1, E[X^2]=1, Var=1-(2*maf-1)^2 = 4*maf*(1-maf)).
#' Wait: for RIL with P(X=1)=maf, P(X=-1)=1-maf:
#'   E[X] = maf - (1-maf) = 2*maf - 1
#'   E[X^2] = 1
#'   Var(X) = 1 - (2*maf - 1)^2 = 4*maf*(1-maf)
#' For F2 with dosage 0/1/2 converted to -1/0/1:
#'   This is just a linear transform of 0/1/2: X_f2 = X_raw - 1
#'   Var(X_f2) = Var(X_raw) = 2*maf*(1-maf)
#'
#' @param pve Target phenotypic variance explained.
#' @param maf Minor allele frequency. Default 0.3.
#' @param population Population type. Default "RIL".
#' @return Standard deviation of effect size.
#' @keywords internal
sd_from_pve <- function(pve, maf = 0.3, population = c("RIL", "F2")) {
  population <- match.arg(population)
  if (population == "RIL") {
    var_x <- 4 * maf * (1 - maf)
  } else {
    var_x <- 2 * maf * (1 - maf)
  }
  sqrt(pve / var_x)
}


# ------------------------------------------------------------------------------
# 2. Core Phenotype Generation Models
# ------------------------------------------------------------------------------

#' Base Model M1 (Additive + Optional AA Epistasis)
#'
#' Generative model:
#'   Y = X %*% B + sum_{(j,k)} aa_{jk} * (X[,j] * X[,k]) + E
#'
#' Uses deterministic expected variance to control residual, avoiding small-sample fluctuations.
#'
#' @param X Genotype matrix, n x p. Coding: -1/1 (RIL) or -1/0/1 (F2).
#' @param B SNP effect matrix, p x m. Additive effects only.
#' @param Sigma_E m x m residual covariance matrix.
#' @param epi_pairs List of epistatic SNP pairs. Each element: c(snp_j, snp_k) (1-based indices).
#'   Default NULL (no epistasis).
#' @param epi_effects Vector of epistatic effect sizes, same length as epi_pairs. Default NULL.
#' @param maf Minor allele frequency (for expected variance calculation). Default 0.3.
#' @param population Population type: "RIL" or "F2". Default "RIL".
#' @param seed Random seed.
#' @export
sim_phenotype_M1 <- function(X, B, Sigma_E,
                              epi_pairs = NULL, epi_effects = NULL,
                              maf = 0.3, population = c("RIL", "F2"),
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  m <- ncol(B)
  population <- match.arg(population)

  stopifnot(ncol(X) == nrow(B))
  stopifnot(dim(Sigma_E) == c(m, m))

  # 1. Additive genetic effects
  G_add <- X %*% B

  # 2. Epistatic effects (AA: additive x additive)
  G_epi <- matrix(0, n, m)
  if (!is.null(epi_pairs) && length(epi_pairs) > 0) {
    if (is.null(epi_effects)) {
      stop("epi_effects must be provided when epi_pairs is non-NULL")
    }
    if (length(epi_pairs) != length(epi_effects)) {
      stop("epi_pairs and epi_effects must have the same length")
    }
    for (i in seq_along(epi_pairs)) {
      pair <- epi_pairs[[i]]
      eff <- epi_effects[i]
      # Hadamard product of the two SNP columns
      interaction <- X[, pair[1]] * X[, pair[2]]
      G_epi <- G_epi + eff * interaction
    }
  }

  G <- G_add + G_epi

  # 3. Expected variance control (deterministic)
  # Additive variance
  expected_var_G_add <- colSums(B^2) * ifelse(population == "RIL", 4 * maf * (1 - maf), 2 * maf * (1 - maf))

  # Epistatic variance (approximate: assuming independence between pairs)
  expected_var_G_epi <- rep(0, m)
  if (!is.null(epi_pairs) && length(epi_pairs) > 0) {
    for (i in seq_along(epi_pairs)) {
      pair <- epi_pairs[[i]]
      eff <- epi_effects[i]
      # Var(X_j * X_k) for independent SNPs
      if (population == "RIL") {
        var_interaction <- 1  # (X_j*X_k) is -1 or 1, E=0, Var=1
      } else {
        # F2: X_j, X_k in {-1,0,1}. For independent SNPs:
        # E[X_j*X_k] = E[X_j]*E[X_k] = 0
        # E[(X_j*X_k)^2] = E[X_j^2]*E[X_k^2] = (2*maf*(1-maf)) * (2*maf*(1-maf))... wait
        # Actually for -1/0/1 with p(1)=maf^2, p(0)=2*maf*(1-maf), p(-1)=(1-maf)^2:
        # E[X^2] = 1*maf^2 + 0 + 1*(1-maf)^2 = maf^2 + (1-maf)^2
        # This is NOT 2*maf*(1-maf). Let me recalculate.
        # For dosage coding 0/1/2: E[X_raw]=2*maf, Var=2*maf*(1-maf)
        # X_f2 = X_raw - 1, so E[X_f2] = 2*maf - 1, Var = 2*maf*(1-maf)
        # E[X_f2^2] = Var + E^2 = 2*maf*(1-maf) + (2*maf-1)^2 = 2maf-2maf^2 + 4maf^2-4maf+1 = 2maf^2-2maf+1
        # Hmm, simpler: E[X_f2^2] = 1*p(1 or -1) + 0 = 1*(maf^2 + (1-maf)^2) = 2maf^2 - 2maf + 1
        # For independent X_j, X_k: E[(X_j*X_k)^2] = E[X_j^2]*E[X_k^2] = (2maf^2-2maf+1)^2
        var_interaction <- (2 * maf^2 - 2 * maf + 1)^2
      }
      expected_var_G_epi <- expected_var_G_epi + eff^2 * var_interaction
    }
  }

  expected_var_G <- expected_var_G_add + expected_var_G_epi
  sigma_eps <- pmax(sqrt(1 - expected_var_G), 0.1)

  D <- diag(sigma_eps)
  Cor_E <- stats::cov2cor(Sigma_E)
  Sigma_E_actual <- D %*% Cor_E %*% D

  E <- mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = Sigma_E_actual)

  # 4. Phenotype
  Y <- G + E

  # 5. Validation
  vars <- apply(Y, 2, stats::var)
  if (any(abs(vars - 1) > 0.15)) {
    warning(sprintf("M1: Trait variance deviates from 1: %s", paste(round(vars, 3), collapse = ", ")))
  }

  list(
    Y = Y,
    G = G,
    G_add = G_add,
    G_epi = G_epi,
    E = E,
    B = B,
    epi_pairs = epi_pairs,
    epi_effects = epi_effects,
    Sigma_E = Sigma_E_actual,
    expected_var_G = expected_var_G,
    params = list(model = "M1", target_var = 1, actual_var = as.numeric(vars))
  )
}


#' Causal Chain Model M2 (Additive + Optional AA Epistasis + Causal Path)
#'
#' Generative model:
#'   Y_base = X %*% B + sum_{(j,k)} aa_{jk} * (X[,j] * X[,k]) + E
#'   Y = Y_base %*% t(solve(I_m - Tau))
#'
#' @param X Genotype matrix, n x p.
#' @param B SNP effect matrix, p x m.
#' @param Tau m x m causal path matrix.
#' @param Sigma_E m x m residual covariance matrix.
#' @param epi_pairs List of epistatic SNP pairs. Default NULL.
#' @param epi_effects Vector of epistatic effect sizes. Default NULL.
#' @param maf Minor allele frequency. Default 0.3.
#' @param population Population type: "RIL" or "F2". Default "RIL".
#' @param seed Random seed.
#' @export
sim_phenotype_M2 <- function(X, B, Tau, Sigma_E,
                              epi_pairs = NULL, epi_effects = NULL,
                              maf = 0.3, population = c("RIL", "F2"),
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  m <- ncol(B)
  population <- match.arg(population)

  stopifnot(ncol(X) == nrow(B))
  stopifnot(dim(Tau) == c(m, m))
  stopifnot(dim(Sigma_E) == c(m, m))

  # 1. Additive effects
  G_add <- X %*% B

  # 2. Epistatic effects
  G_epi <- matrix(0, n, m)
  if (!is.null(epi_pairs) && length(epi_pairs) > 0) {
    if (is.null(epi_effects)) stop("epi_effects must be provided when epi_pairs is non-NULL")
    if (length(epi_pairs) != length(epi_effects)) stop("epi_pairs and epi_effects must have the same length")
    for (i in seq_along(epi_pairs)) {
      pair <- epi_pairs[[i]]
      eff <- epi_effects[i]
      interaction <- X[, pair[1]] * X[, pair[2]]
      G_epi <- G_epi + eff * interaction
    }
  }

  G <- G_add + G_epi

  # 3. Expected variance control
  expected_var_G_add <- colSums(B^2) * ifelse(population == "RIL", 4 * maf * (1 - maf), 2 * maf * (1 - maf))

  expected_var_G_epi <- rep(0, m)
  if (!is.null(epi_pairs) && length(epi_pairs) > 0) {
    for (i in seq_along(epi_pairs)) {
      pair <- epi_pairs[[i]]
      eff <- epi_effects[i]
      if (population == "RIL") {
        var_interaction <- 1
      } else {
        var_interaction <- (2 * maf^2 - 2 * maf + 1)^2
      }
      expected_var_G_epi <- expected_var_G_epi + eff^2 * var_interaction
    }
  }

  expected_var_G <- expected_var_G_add + expected_var_G_epi
  sigma_eps <- pmax(sqrt(1 - expected_var_G), 0.1)

  D <- diag(sigma_eps)
  Cor_E <- stats::cov2cor(Sigma_E)
  Sigma_E_actual <- D %*% Cor_E %*% D

  E <- mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = Sigma_E_actual)

  # 4. Base phenotype
  Y_base <- G + E

  # 5. Causal transformation
  I_T <- diag(m) - Tau
  I_T_inv <- solve(I_T)
  Y <- Y_base %*% t(I_T_inv)

  # 6. Validation
  vars <- apply(Y, 2, stats::var)
  if (any(abs(vars - 1) > 0.5)) {
    warning(sprintf("M2: Post-causal trait variance: %s. Tau-induced scale change is normal.",
                    paste(round(vars, 3), collapse = ", ")))
  }

  tau_observed <- NULL
  if (m == 2 && abs(Tau[2, 1]) > 1e-10 && abs(Tau[1, 2]) < 1e-10) {
    fit_check <- stats::lm(Y[, 2] ~ Y[, 1])
    tau_obs <- unname(stats::coef(fit_check)[2])
    tau_true <- Tau[2, 1]
    tau_observed <- tau_obs
    if (abs(tau_obs - tau_true) > 0.15) {
      warning(sprintf("M2: Observed tau (%.3f) deviates from true tau (%.3f).", tau_obs, tau_true))
    }
  }

  list(
    Y = Y,
    Y_base = Y_base,
    G = G,
    G_add = G_add,
    G_epi = G_epi,
    E = E,
    B = B,
    Tau = Tau,
    epi_pairs = epi_pairs,
    epi_effects = epi_effects,
    Sigma_E = Sigma_E_actual,
    expected_var_G = expected_var_G,
    params = list(model = "M2", target_var = 1, actual_var = as.numeric(vars), tau_observed = tau_observed)
  )
}


# ------------------------------------------------------------------------------
# 3. Five-Class Genetic Mechanism Generators
# ------------------------------------------------------------------------------

#' @export
generate_class1 <- function(n, p = 200, pve = 0.01, maf = 0.3,
                             population = "RIL", seed = NULL,
                             epi_pairs = NULL, epi_effects = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf, population = population)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf, population))
  out <- sim_phenotype_M1(X, B, Sigma_E = diag(2),
                           epi_pairs = epi_pairs, epi_effects = epi_effects,
                           maf = maf, population = population)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class1", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class1", seed = seed))
}

#' @export
generate_class3 <- function(n, p = 200, pve = 0.01, maf = 0.3,
                             population = "RIL", seed = NULL,
                             epi_pairs = NULL, epi_effects = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf, population = population)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf, population))
  B[1:30, 2] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf, population))
  out <- sim_phenotype_M1(X, B, Sigma_E = diag(2),
                           epi_pairs = epi_pairs, epi_effects = epi_effects,
                           maf = maf, population = population)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class3", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class3", seed = seed))
}

#' @export
generate_class4 <- function(n, p = 200, pve_A = 0.02, tau = 0.3, maf = 0.3,
                             population = "RIL", seed = NULL,
                             epi_pairs = NULL, epi_effects = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf, population = population)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve_A, maf, population))
  Tau <- matrix(0, 2, 2); Tau[2, 1] <- tau
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E = diag(2),
                             epi_pairs = epi_pairs, epi_effects = epi_effects,
                             maf = maf, population = population)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class4", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], tau = tau, stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class4", seed = seed))
}

#' @export
generate_class5 <- function(n, p = 200, pve_A = 0.02, pve_B = 0.003,
                             tau = 0.3, maf = 0.3,
                             population = "RIL", seed = NULL,
                             epi_pairs = NULL, epi_effects = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf, population = population)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve_A, maf, population))
  B[1:30, 2] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve_B, maf, population))
  Tau <- matrix(0, 2, 2); Tau[2, 1] <- tau
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E = diag(2),
                             epi_pairs = epi_pairs, epi_effects = epi_effects,
                             maf = maf, population = population)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class5", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], tau = tau, stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class5", seed = seed))
}

#' @export
generate_class6 <- function(n, p = 200, pve = 0.01,
                             tau_AB = 0.2, tau_BA = 0.2, maf = 0.3,
                             population = "RIL", seed = NULL,
                             epi_pairs = NULL, epi_effects = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf, population = population)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf, population))
  Tau <- matrix(0, 2, 2)
  Tau[2, 1] <- tau_AB
  Tau[1, 2] <- tau_BA
  out <- sim_phenotype_M2(X, B, Tau, Sigma_E = diag(2),
                             epi_pairs = epi_pairs, epi_effects = epi_effects,
                             maf = maf, population = population)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class6", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2],
    tau_AB = tau_AB, tau_BA = tau_BA, stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "class6", seed = seed))
}

# ------------------------------------------------------------------------------
# 4. Class 2 Stress Test Generator
# ------------------------------------------------------------------------------

#' @export
generate_covariance_stress <- function(n, p = 200, pve = 0.01,
                                        rho = 0.5, maf = 0.3,
                                        population = "RIL", seed = NULL,
                                        epi_pairs = NULL, epi_effects = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- sim_genotype(n, p, maf = maf, population = population)
  B <- matrix(0, p, 2)
  B[1:30, 1] <- stats::rnorm(30, mean = 0, sd = sd_from_pve(pve, maf, population))
  Sigma_E <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
  out <- sim_phenotype_M1(X, B, Sigma_E,
                           epi_pairs = epi_pairs, epi_effects = epi_effects,
                           maf = maf, population = population)
  truth <- data.frame(
    SNP = colnames(X), class = c(rep("class1", 30), rep("null", p - 30)),
    beta_A = B[, 1], beta_B = B[, 2], rho = rho, stringsAsFactors = FALSE
  )
  c(out, list(X = X, truth = truth, class_label = "stress_test", seed = seed))
}

# ------------------------------------------------------------------------------
# 5. Unified Entry Function
# ------------------------------------------------------------------------------

#' @export
generate_condped <- function(scenario = c("class1", "class3", "class4", "class5", "class6", "stress_test"),
                              n, p = 200, ...) {
  scenario <- match.arg(scenario)
  switch(scenario,
    class1 = generate_class1(n = n, p = p, ...),
    class3 = generate_class3(n = n, p = p, ...),
    class4 = generate_class4(n = n, p = p, ...),
    class5 = generate_class5(n = n, p = p, ...),
    class6 = generate_class6(n = n, p = p, ...),
    stress_test = generate_covariance_stress(n = n, p = p, ...)
  )
}