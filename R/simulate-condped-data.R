#' Simulate multi-trait data with known ground truth
#'
#' Generates simulated data with a known multi-trait effect matrix, true
#' associated-trait sets, true conditional-deviation sets, theoretical PVE
#' and complete covariance ground truth, following the interface contract
#' (section 4.1) and the Methods Simulation Study (sections 5.2-5.4).
#'
#' @details
#' Markers are generated independently (no LD blocks). With
#' `structured = TRUE`, individuals belong to `n_groups` families or
#' subpopulations whose allele frequencies are mildly differentiated
#' (Balding-Nichols model with parameter `fst`), so that the genomic
#' relationship matrix keeps an identifiable eigenvalue dispersion instead
#' of degenerating to the identity. The background relationship matrix is
#' built from the standardised background genotypes and scaled to
#' `mean(diag(K_bg)) = 1`. The polygenic background is drawn through the
#' eigen decomposition of `K_bg`.
#'
#' The focal QTL covariance
#' \eqn{\Sigma_Q = \sum_l v_l \beta_l \beta_l^\top} (with
#' \eqn{v_l = 2 p_l (1 - p_l)}) is subtracted from the target total genetic
#' covariance, \eqn{\Sigma_{G,bg} = \Sigma_{G,total} - \Sigma_Q}. If
#' \eqn{\Sigma_{G,bg}} is not positive semi-definite, the focal loci and
#' their effects are re-drawn (up to `max_attempts` times); as a last
#' resort the focal effects are shrunk by a recorded scalar factor so that
#' PSD holds, and every component of `truth` is updated accordingly.
#'
#' True conditional projection coefficients are computed from the total
#' phenotypic covariance \eqn{\Sigma_{P,total} = \Sigma_{G,total} +
#' \Sigma_E} once, and the true conditional effects satisfy
#' \eqn{\eta_{il} = c_i^\top \beta_l}. In the `covariance_aligned`
#' architecture the target-trait \eqn{\eta} is exactly zero by
#' construction; in `projection_induced` the target-trait \eqn{\beta = 0}
#' while \eqn{\eta \neq 0} under correlated traits.
#'
#' @param n Positive integer; number of individuals. Main simulations use
#'   500, 750, 1000 or 1250.
#' @param m Integer (>= 2); number of traits (default 4).
#' @param p Integer; total number of background markers.
#' @param maf_range Numeric vector of length 2 giving the range of minor
#'   allele frequencies.
#' @param h2 Numeric vector of length `m` (or a recycled scalar); target
#'   SNP heritabilities. Total phenotypic variance per trait is 1.
#' @param R_G Optional `m x m` genetic correlation matrix; `NULL` uses an
#'   exchangeable correlation of 0.4.
#' @param R_E Optional `m x m` residual correlation matrix; `NULL` uses an
#'   exchangeable correlation of 0.4.
#' @param architecture Genetic architecture of the focal locus/loci; one of
#'   `"null"`, `"single_trait"`, `"shared_same"`, `"shared_opposite"`,
#'   `"dense"`, `"covariance_aligned"`, `"conditional_deviation"`,
#'   `"projection_induced"`.
#' @param locus_pve Scalar or per-active-trait vector of target locus PVE.
#'   A single scaling factor per locus preserves the effect-template
#'   ratios, so the target is matched on average over the active traits
#'   (a vector is summarised by its mean).
#' @param delta True deviation effect in the conditional-deviation
#'   architecture; after scaling, the target-trait conditional effect is
#'   exactly `s_l * delta`.
#' @param target_trait Integer index of the target trait.
#' @param n_groups Number of families/subpopulations used when
#'   `structured = TRUE`.
#' @param fst Mild between-group allele frequency differentiation
#'   (Balding-Nichols parameter, must be > 0 when `structured = TRUE`).
#' @param structured Logical; when `TRUE` individuals carry family/subgroup
#'   structure so that the genomic relationship matrix does not degenerate
#'   to the identity.
#' @param n_qtl Integer; number of focal QTL (sampled without replacement).
#' @param max_attempts Maximum number of re-draws when the background
#'   genetic covariance is not positive semi-definite, before falling back
#'   to a recorded scaling of the focal effects.
#' @param seed Optional random seed.
#' @param return_latent Logical; whether to return latent genetic and
#'   residual values.
#'
#' @return A list with components:
#'   \describe{
#'     \item{Y}{`n x m` phenotype matrix.}
#'     \item{W}{`n x 1` intercept design matrix.}
#'     \item{G}{`n x p` genotype dosage matrix (0/1/2).}
#'     \item{K_bg}{`n x n` background genomic relationship matrix,
#'       `mean(diag(K_bg)) = 1`.}
#'     \item{qtl_index}{Integer vector of focal QTL positions in `G`.}
#'     \item{truth}{List with `B_Q` and `beta` (`n_qtl x m` effect
#'       matrices), `A` (per-locus associated-trait sets), `D` (per-locus
#'       conditional-deviation sets, subsets of `A`), `eta` (true
#'       conditional effects), `locus_pve`, `conditional_pve`,
#'       `Sigma_Q`, `Sigma_G_bg`, `Sigma_G_total`, `Sigma_E`,
#'       `Sigma_P_total`, `gamma` (`m x (m-1)`, row `i` is
#'       \eqn{\gamma_{i,-i}}) and `C` (`m x m`, column `i` is the contrast
#'       vector \eqn{c_i}).}
#'     \item{latent}{`list(U, E)` when `return_latent = TRUE`, else
#'       `NULL`.}
#'     \item{generator}{List with `seed`, `attempts`, `scaled`,
#'       `scale_factor` and `settings`.}
#'     \item{status}{Standard CondPED status list.}
#'     \item{diagnostics}{List with `mean_diag_K`, `K_eigen_sd`, `K_rank`
#'       and `min_eigen_Sigma_G_bg`.}
#'   }
#' @export
simulate_condped_data <- function(
  n = 1000L,
  m = 4L,
  p = 2000L,
  maf_range = c(0.05, 0.50),
  h2 = rep(0.50, m),
  R_G = NULL,
  R_E = NULL,
  architecture = c(
    "null",
    "single_trait",
    "shared_same",
    "shared_opposite",
    "dense",
    "covariance_aligned",
    "conditional_deviation",
    "projection_induced"
  ),
  locus_pve = 0.01,
  delta = 0,
  target_trait = 1L,
  n_groups = 8L,
  fst = 0.05,
  structured = TRUE,
  n_qtl = 1L,
  max_attempts = 100L,
  seed = NULL,
  return_latent = TRUE
) {
  architecture <- match.arg(architecture)

  # ---- input validation ----------------------------------------------------
  .check_count(n, "n", min = 2L)
  .check_count(m, "m", min = 2L)
  .check_count(p, "p", min = 1L)
  .check_count(n_qtl, "n_qtl", min = 1L)
  .check_count(max_attempts, "max_attempts", min = 1L)
  .check_count(n_groups, "n_groups", min = 1L)
  .check_count(target_trait, "target_trait", min = 1L)
  if (n_qtl > p) stop("n_qtl must not exceed p.", call. = FALSE)
  if (target_trait > m) stop("target_trait must be between 1 and m.", call. = FALSE)
  if (!is.numeric(maf_range) || length(maf_range) != 2L ||
      any(!is.finite(maf_range)) || maf_range[1L] <= 0 ||
      maf_range[1L] >= maf_range[2L] || maf_range[2L] > 0.5) {
    stop("maf_range must satisfy 0 < maf_range[1] < maf_range[2] <= 0.5.",
         call. = FALSE)
  }
  if (!is.numeric(h2) || !(length(h2) %in% c(1L, m)) ||
      any(!is.finite(h2)) || any(h2 <= 0) || any(h2 >= 1)) {
    stop("h2 must have length 1 or m, with all entries in (0, 1).", call. = FALSE)
  }
  h2 <- rep(h2, length.out = m)
  if (!is.numeric(locus_pve) || any(!is.finite(locus_pve)) || any(locus_pve < 0)) {
    stop("locus_pve must be a non-negative numeric scalar or vector.", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta)) {
    stop("delta must be a single finite number.", call. = FALSE)
  }
  if (structured) {
    if (n_groups < 2L) {
      stop("n_groups must be at least 2 when structured = TRUE.", call. = FALSE)
    }
    if (!is.numeric(fst) || length(fst) != 1L || !is.finite(fst) ||
        fst <= 0 || fst >= 1) {
      stop("fst must be a single number in (0, 1) when structured = TRUE.",
           call. = FALSE)
    }
  }

  if (!is.null(seed)) set.seed(seed)

  trait_names <- paste0("Trait", seq_len(m))
  ind_names <- paste0("Ind", seq_len(n))
  marker_names <- paste0("M", seq_len(p))

  # ---- target covariance structure ------------------------------------------
  if (is.null(R_G)) R_G <- .exchangeable_cor(m, 0.4)
  if (is.null(R_E)) R_E <- .exchangeable_cor(m, 0.4)
  .check_correlation_matrix(R_G, m, "R_G")
  .check_correlation_matrix(R_E, m, "R_E")

  Sigma_G_total <- .cov_from_h2(h2, R_G)
  Sigma_E <- .cov_from_h2(1 - h2, R_E)
  Sigma_P_total <- Sigma_G_total + Sigma_E
  dimnames(Sigma_G_total) <- dimnames(Sigma_E) <- dimnames(Sigma_P_total) <-
    list(trait_names, trait_names)

  # True conditional contrasts from the total phenotypic covariance (once).
  ctr <- .contrasts_from_cov(Sigma_P_total, trait_names)

  # ---- genotypes: independent markers, optional group structure -------------
  maf <- stats::runif(p, min = maf_range[1L], max = maf_range[2L])
  G <- .sim_genotypes(n, maf, structured = structured,
                      n_groups = n_groups, fst = fst)
  dimnames(G) <- list(ind_names, marker_names)

  # Background relationship matrix (focal QTL excluded) via eigen-based GRM.
  psd_tol <- 1e-10

  # ---- focal QTL effects with PSD enforcement --------------------------------
  attempt <- 0L
  scaled <- FALSE
  scale_factor <- 1
  repeat {
    attempt <- attempt + 1L
    qtl_index <- sort(sample.int(p, n_qtl))
    v_qtl <- 2 * maf[qtl_index] * (1 - maf[qtl_index])
    eff <- .focal_effects(
      architecture = architecture, m = m, target_trait = target_trait,
      gamma = ctr$gamma, delta = delta, v_qtl = v_qtl,
      locus_pve = locus_pve, scale = scale_factor
    )
    Sigma_Q <- crossprod(eff$beta, v_qtl * eff$beta)
    Sigma_G_bg <- Sigma_G_total - Sigma_Q
    min_eigen_bg <- .min_eigen_sym(Sigma_G_bg)
    if (min_eigen_bg >= -psd_tol) break
    if (attempt >= max_attempts) {
      # Last resort: deterministic shrinkage of the focal effects, fully
      # recorded; all downstream truth components are recomputed from the
      # scaled effects below.
      lambda <- .rescale_to_psd(Sigma_G_total, Sigma_Q, tol = psd_tol)
      scale_factor <- scale_factor * sqrt(lambda)
      scaled <- TRUE
      eff <- .focal_effects(
        architecture = architecture, m = m, target_trait = target_trait,
        gamma = ctr$gamma, delta = delta, v_qtl = v_qtl,
        locus_pve = locus_pve, scale = scale_factor
      )
      Sigma_Q <- crossprod(eff$beta, v_qtl * eff$beta)
      Sigma_G_bg <- Sigma_G_total - Sigma_Q
      min_eigen_bg <- .min_eigen_sym(Sigma_G_bg)
      break
    }
  }

  beta <- eff$beta
  dimnames(beta) <- list(marker_names[qtl_index], trait_names)

  # ---- background GRM and latent components ----------------------------------
  bg <- setdiff(seq_len(p), qtl_index)
  maf_bg <- maf[bg]
  Z_bg <- sweep(G[, bg, drop = FALSE], 2L, 2 * maf_bg, `-`)
  Z_bg <- sweep(Z_bg, 2L, sqrt(2 * maf_bg * (1 - maf_bg)), `/`)
  K_bg <- .make_grm(Z_bg)
  dimnames(K_bg) <- list(ind_names, ind_names)

  eig_K <- eigen(K_bg, symmetric = TRUE)
  lambda_K <- eig_K$values
  K_rank <- sum(lambda_K > sqrt(.Machine$double.eps) * max(lambda_K))

  # Polygenic background U ~ MN(0, K_bg, Sigma_G_bg) via the eigen
  # decomposition of K_bg; residual E ~ MN(0, I_n, Sigma_E).
  S_G <- .symmetric_psd_sqrt(Sigma_G_bg)
  S_E <- .symmetric_psd_sqrt(Sigma_E)
  Z_u <- matrix(stats::rnorm(n * m), nrow = n, ncol = m)
  U <- sweep(eig_K$vectors, 2L, sqrt(pmax(lambda_K, 0)), `*`) %*% Z_u %*% S_G
  E <- matrix(stats::rnorm(n * m), nrow = n, ncol = m) %*% S_E
  dimnames(U) <- dimnames(E) <- list(ind_names, trait_names)

  # ---- phenotypes -------------------------------------------------------------
  W <- matrix(1, nrow = n, ncol = 1,
              dimnames = list(ind_names, "Intercept"))
  Y <- G[, qtl_index, drop = FALSE] %*% beta + U + E

  # ---- ground truth -------------------------------------------------------------
  A <- lapply(seq_len(n_qtl), function(l) which(beta[l, ] != 0))
  eta <- beta %*% ctr$C
  D <- lapply(seq_len(n_qtl), function(l) {
    intersect(A[[l]], which(abs(eta[l, ]) > 1e-10))
  })
  sigma2_P <- diag(Sigma_P_total)
  locus_pve_true <- sweep(beta^2, 2L, sigma2_P, `/`) * v_qtl
  cond_pve_true <- sweep(eta^2, 2L, ctr$conditional_variance, `/`) * v_qtl
  dimnames(eta) <- dimnames(locus_pve_true) <- dimnames(cond_pve_true) <-
    dimnames(beta)

  truth <- list(
    B_Q = beta,
    beta = beta,
    A = A,
    D = D,
    eta = eta,
    locus_pve = locus_pve_true,
    conditional_pve = cond_pve_true,
    Sigma_Q = Sigma_Q,
    Sigma_G_bg = Sigma_G_bg,
    Sigma_G_total = Sigma_G_total,
    Sigma_E = Sigma_E,
    Sigma_P_total = Sigma_P_total,
    gamma = ctr$gamma,
    C = ctr$C
  )

  ok <- min_eigen_bg >= -psd_tol
  status <- .new_status(
    ok = ok,
    code = if (ok) "ok" else "unstable",
    message = if (ok) "" else "Sigma_G_bg is not PSD even after focal-effect scaling.",
    warnings = if (scaled) {
      sprintf(
        paste0(
          "Sigma_G_bg was not PSD after %d re-draw(s); focal effects were ",
          "shrunk by a cumulative factor of %.4g and all truth components ",
          "were updated accordingly."
        ),
        attempt, scale_factor
      )
    } else {
      character()
    }
  )

  list(
    Y = Y,
    W = W,
    G = G,
    K_bg = K_bg,
    qtl_index = qtl_index,
    truth = truth,
    latent = if (isTRUE(return_latent)) list(U = U, E = E) else NULL,
    generator = list(
      seed = seed,
      attempts = attempt,
      scaled = scaled,
      scale_factor = scale_factor,
      settings = list(
        n = n, m = m, p = p, maf_range = maf_range, h2 = h2,
        architecture = architecture, locus_pve = locus_pve, delta = delta,
        target_trait = target_trait, n_groups = n_groups, fst = fst,
        structured = structured, n_qtl = n_qtl, max_attempts = max_attempts,
        return_latent = return_latent
      )
    ),
    status = status,
    diagnostics = list(
      mean_diag_K = mean(diag(K_bg)),
      K_eigen_sd = stats::sd(lambda_K),
      K_rank = K_rank,
      min_eigen_Sigma_G_bg = min_eigen_bg
    )
  )
}

# ---- internal helpers ---------------------------------------------------------

#' Check a positive integer-like scalar argument
#' @param x Value to check.
#' @param name Argument name used in the error message.
#' @param min Minimum allowed value.
#' @return Invisible `TRUE`; stops otherwise.
#' @keywords internal
.check_count <- function(x, name, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != round(x) || x < min) {
    stop(sprintf("%s must be a single integer >= %d.", name, min), call. = FALSE)
  }
  invisible(TRUE)
}

#' Exchangeable correlation matrix
#' @param m Dimension.
#' @param r Common correlation.
#' @return An `m x m` correlation matrix.
#' @keywords internal
.exchangeable_cor <- function(m, r) {
  matrix(r, m, m) + diag(1 - r, m)
}

#' Validate a correlation matrix argument
#' @param R Candidate matrix.
#' @param m Expected dimension.
#' @param name Argument name used in error messages.
#' @return Invisible `TRUE`; stops otherwise.
#' @keywords internal
.check_correlation_matrix <- function(R, m, name) {
  if (!is.matrix(R) || !is.numeric(R) || nrow(R) != m || ncol(R) != m) {
    stop(sprintf("%s must be a numeric %d x %d matrix.", name, m, m),
         call. = FALSE)
  }
  if (anyNA(R) || max(abs(R - t(R))) > 1e-10) {
    stop(sprintf("%s must be symmetric and free of missing values.", name),
         call. = FALSE)
  }
  if (max(abs(diag(R) - 1)) > 1e-10) {
    stop(sprintf("%s must have unit diagonal.", name), call. = FALSE)
  }
  if (.min_eigen_sym(R) <= 0) {
    stop(sprintf("%s must be positive definite.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Build a covariance matrix from marginal variances and a correlation matrix
#' @param var Marginal variances (length m).
#' @param R Correlation matrix (m x m).
#' @return The covariance matrix `diag(sqrt(var)) R diag(sqrt(var))`.
#' @keywords internal
.cov_from_h2 <- function(var, R) {
  R * outer(sqrt(var), sqrt(var))
}

#' Smallest eigenvalue of a symmetric matrix
#' @param A Symmetric matrix.
#' @return Numeric scalar.
#' @keywords internal
.min_eigen_sym <- function(A) {
  min(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
}

#' Symmetric square root of a PSD matrix
#' @param Sigma Symmetric positive semi-definite matrix.
#' @return Symmetric matrix `S` with `S %*% S = Sigma`; tiny negative
#'   eigenvalues are floored at zero.
#' @keywords internal
.symmetric_psd_sqrt <- function(Sigma) {
  eig <- eigen(Sigma, symmetric = TRUE)
  sweep(eig$vectors, 2L, sqrt(pmax(eig$values, 0)), `*`) %*% t(eig$vectors)
}

#' True conditional contrasts from a phenotypic covariance matrix
#'
#' Computes, for every trait `i`, the projection coefficients
#' \eqn{\gamma_{i,-i} = \Sigma_{P,-i,-i}^{-1} \Sigma_{P,-i,i}}, the contrast
#' vectors \eqn{c_i = e_i - E_{-i} \gamma_{i,-i}} and the conditional
#' variances (Schur complements). All inversions go through
#' [`.safe_inverse()`].
#'
#' @param Sigma_P `m x m` phenotypic covariance matrix.
#' @param trait_names Character vector of trait names.
#' @return List with `gamma` (`m x (m-1)`, row `i` is
#'   \eqn{\gamma_{i,-i}}), `C` (`m x m`, column `i` is \eqn{c_i}) and
#'   `conditional_variance` (length `m`).
#' @keywords internal
.contrasts_from_cov <- function(Sigma_P, trait_names = colnames(Sigma_P)) {
  m <- ncol(Sigma_P)
  gamma <- matrix(0, nrow = m, ncol = m - 1L,
                  dimnames = list(trait_names, NULL))
  C <- diag(m)
  dimnames(C) <- list(trait_names, trait_names)
  conditional_variance <- numeric(m)
  for (i in seq_len(m)) {
    inv <- .safe_inverse(Sigma_P[-i, -i, drop = FALSE])
    if (inv$status == "failed") {
      stop(sprintf("Failed to invert the phenotypic covariance submatrix for trait %d.", i),
           call. = FALSE)
    }
    g <- drop(inv$inverse %*% Sigma_P[-i, i])
    gamma[i, ] <- g
    C[-i, i] <- -g
    conditional_variance[i] <- Sigma_P[i, i] - sum(Sigma_P[i, -i] * g)
  }
  list(gamma = gamma, C = C, conditional_variance = conditional_variance)
}

#' Generate genotype dosages with optional group structure
#'
#' Markers are independent (no LD blocks). With `structured = TRUE`,
#' individuals are split into `n_groups` groups whose allele frequencies
#' follow a Balding-Nichols model around the base MAF with parameter
#' `fst`, giving the relationship matrix a non-degenerate eigenvalue
#' dispersion.
#'
#' @param n Number of individuals.
#' @param maf Base minor allele frequencies (length p).
#' @param structured Logical; use group structure.
#' @param n_groups Number of groups.
#' @param fst Balding-Nichols differentiation parameter.
#' @return An `n x p` dosage matrix (0/1/2).
#' @keywords internal
.sim_genotypes <- function(n, maf, structured, n_groups, fst) {
  p <- length(maf)
  if (!structured) {
    return(matrix(
      stats::rbinom(n * p, size = 2, prob = rep(maf, each = n)),
      nrow = n, ncol = p
    ))
  }
  sizes <- rep(n %/% n_groups, n_groups)
  remainder <- n %% n_groups
  if (remainder > 0L) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
  group <- rep(seq_len(n_groups), times = sizes)
  shape1 <- maf * (1 - fst) / fst
  shape2 <- (1 - maf) * (1 - fst) / fst
  p_group <- matrix(
    stats::rbeta(p * n_groups, rep(shape1, times = n_groups),
                 rep(shape2, times = n_groups)),
    nrow = n_groups, byrow = TRUE
  )
  matrix(
    stats::rbinom(n * p, size = 2, prob = p_group[group, ]),
    nrow = n, ncol = p
  )
}

#' Raw focal-QTL effect template for one architecture
#'
#' @param architecture Architecture label.
#' @param m Number of traits.
#' @param target_trait Target trait index.
#' @param gamma True projection coefficients (`m x (m-1)`).
#' @return Numeric effect template of length `m`.
#' @keywords internal
.effect_template <- function(architecture, m, target_trait, gamma) {
  neighbour <- if (target_trait == m) 1L else target_trait + 1L
  switch(
    architecture,
    null = numeric(m),
    single_trait = {
      a <- numeric(m)
      a[target_trait] <- 1
      a
    },
    shared_same = {
      a <- numeric(m)
      a[c(target_trait, neighbour)] <- 1
      a
    },
    shared_opposite = {
      a <- numeric(m)
      a[target_trait] <- 1
      a[neighbour] <- -1
      a
    },
    dense = {
      order_from_target <- (seq_len(m) - target_trait) %% m
      1 - 0.2 * order_from_target
    },
    covariance_aligned = ,
    conditional_deviation = {
      # Template of the Methods (Q7/Q8): the target-trait component equals
      # the projection of the remaining components, so that the true
      # conditional effect on the target trait is exactly zero before any
      # deviation is added.
      others <- setdiff(seq_len(m), target_trait)
      a <- numeric(m)
      a[others] <- 1 - 0.2 * seq(0, length.out = m - 1L)
      a[target_trait] <- sum(gamma[target_trait, ] * a[others])
      a
    },
    projection_induced = {
      a <- numeric(m)
      a[neighbour] <- 1
      a
    }
  )
}

#' Scaled focal-QTL effects with exact projection identities
#'
#' Scales the effect template so that the average PVE over the
#' template-active traits equals the target (Methods eq. 69), and enforces
#' the exact conditional identities for the projection-based
#' architectures: `covariance_aligned` sets the target-trait effect to
#' exactly \eqn{\gamma^\top \beta_{-t}}; `conditional_deviation` adds
#' exactly `s_l * delta` on top.
#'
#' @param architecture Architecture label.
#' @param m Number of traits.
#' @param target_trait Target trait index.
#' @param gamma True projection coefficients (`m x (m-1)`).
#' @param delta Deviation effect (template units).
#' @param v_qtl Expected genotype variances of the focal loci.
#' @param locus_pve Target average PVE (scalar or vector; the mean is used).
#' @param scale Additional scalar factor applied to all loci (used by the
#'   PSD shrinkage fallback).
#' @return List with `beta` (`n_qtl x m` effect matrix) and `s` (per-locus
#'   cumulative scaling factors, including `scale`).
#' @keywords internal
.focal_effects <- function(architecture, m, target_trait, gamma, delta,
                           v_qtl, locus_pve, scale = 1) {
  a <- .effect_template(architecture, m, target_trait, gamma)
  active <- which(a != 0)
  n_qtl <- length(v_qtl)
  beta <- matrix(0, nrow = n_qtl, ncol = m)
  s <- numeric(n_qtl)
  pi_target <- mean(locus_pve)
  if (length(active) > 0L && pi_target > 0) {
    for (l in seq_len(n_qtl)) {
      s[l] <- scale * sqrt(pi_target / (v_qtl[l] * mean(a[active]^2)))
      beta[l, ] <- s[l] * a
    }
  }
  if (architecture %in% c("covariance_aligned", "conditional_deviation")) {
    g_t <- gamma[target_trait, ]
    for (l in seq_len(n_qtl)) {
      deviation <- if (architecture == "conditional_deviation") s[l] * delta else 0
      beta[l, target_trait] <- sum(g_t * beta[l, -target_trait]) + deviation
    }
  }
  list(beta = beta, s = s)
}

#' Largest shrinkage factor keeping the background covariance PSD
#'
#' Finds, by bisection, the largest lambda between 0 and 1 such that
#' `Sigma_G_total - lambda * Sigma_Q` is positive semi-definite. The
#' smallest eigenvalue is concave in lambda (minimum of linear functions),
#' so bisection is valid.
#'
#' @param Sigma_G_total Target total genetic covariance (PD).
#' @param Sigma_Q Focal QTL covariance (PSD).
#' @param tol Numerical PSD tolerance.
#' @return Scalar lambda in (0, 1].
#' @keywords internal
.rescale_to_psd <- function(Sigma_G_total, Sigma_Q, tol = 1e-10) {
  f <- function(lambda) .min_eigen_sym(Sigma_G_total - lambda * Sigma_Q)
  if (f(1) >= -tol) return(1)
  lo <- 0
  hi <- 1
  for (i in seq_len(60L)) {
    mid <- (lo + hi) / 2
    if (f(mid) >= 0) lo <- mid else hi <- mid
  }
  lo
}
