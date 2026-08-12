#' Simulate CondPED data with population structural truth (v1.0)
#'
#' Generates individual-level multi-trait data with one simulated
#' causal variant and the full population-level structural truth: the
#' candidate trait set, the complete subset loss table, the minimum
#' representative trait sets and the irreducible trait modules. All
#' truth is computed with the SAME deterministic loss path used by
#' estimation ([`.compute_subset_loss()`] and the Stage-2 extractors),
#' evaluated on population parameters.
#'
#' @details
#' The data-generating model is
#' \deqn{Y = WB + \tilde{x}_l\beta_l^{\top} + U + E}
#' with `U ~ MN(0, K_bg, Sigma_G_bg)` and `E ~ MN(0, I_n, Sigma_E)`.
#' The causal variant never enters the background GRM, and
#' `Sigma_G_bg = Sigma_G_total - v_l beta beta'`; only generations
#' with positive semi-definite `Sigma_G_bg` are accepted. Effect
#' vectors are scaled as a whole to the target PVE; the representation
#' loss is invariant to this scaling, and the truth tables are
#' computed AFTER scaling from the final beta and Sigma.
#'
#' `representative_singleton`/`representative_pair` use the closed-form
#' construction of Methods eqs. (39)-(41): with
#' `R = A \ S*`,
#' \deqn{\beta_R = \Gamma_{R|S*}^{\top}\beta_{S*} + \delta d,\qquad
#'   d^{\top}\Omega^{-1}d = 1,\qquad
#'   \delta = \sqrt{\frac{r\,\beta_{S*}^{\top}\Sigma_{S*S*}^{-1}\beta_{S*}}
#'   {1-r}}}
#' so that the target set has population loss exactly r = target_loss.
#' All acceptance rules read population beta, Sigma and the truth
#' subset table only -- never sample p-values or fitted results.
#'
#' @param n Number of individuals.
#' @param m Number of traits.
#' @param p Number of markers.
#' @param maf_range Range for base minor allele frequencies.
#' @param h2 Trait heritabilities (length 1 or `m`).
#' @param R_G,R_E Optional `m x m` correlation matrices; default from
#'   `correlation`.
#' @param architecture One of the ten frozen scenarios; see the
#'   interface contract.
#' @param effect_direction Sign pattern of the candidate effect
#'   template: `"concordant"` (all same direction), `"discordant"`
#'   (all non-leading effects flipped), `"mixed"` (alternating signs,
#'   identical to `"discordant"` for a pair). Only meaningful for
#'   `candidate_pair` and `candidate_dense`; must be `"default"` for
#'   every other architecture.
#' @param locus_pve Target mean per-candidate-trait PVE of the causal
#'   variant.
#' @param tolerance Primary loss tolerance tau.
#' @param target_loss Closed-form target loss `r`; only used by the
#'   representative scenarios.
#' @param target_representative_set Optional target set for the
#'   representative scenarios (defaults: trait 1 / traits 1-2).
#' @param target_irreducible_modules Optional list of module sets for
#'   the irreducible scenarios (defaults: trait 1 / traits 1-2 /
#'   `{1}` and `{3..m}`).
#' @param correlation `"block"` (consecutive trait pairs form modules
#'   with within-module correlation 0.6 and between-module 0.2) or
#'   `"independent"`. Note that joint irreducible
#'   modules (`irreducible_pair`, `multiple_modules`) require correlated
#'   traits; under `"independent"` their acceptance rule correctly
#'   never passes and the function returns status `unstable`.
#' @param structured Logical; group structure for GRM eigenvalue
#'   dispersion.
#' @param n_groups,fst Balding-Nichols group structure parameters.
#' @param truth_margin Minimum allowed distance between any subset
#'   loss and `tolerance` in the structural scenarios.
#' @param max_attempts Maximum number of generation attempts.
#' @param seed Optional random seed.
#' @param return_latent Logical; keep `U` and `E` (default `FALSE`).
#'
#' @return A list with components `Y`, `W`, `G`, `K_bg`,
#'   `causal_index`, `truth` (list with `beta`, `candidate_traits`,
#'   `subset_table`, `minimum_representative_sets`,
#'   `irreducible_modules`, `tolerance_path`, `locus_pve`,
#'   `target_loss`, `delta`, `Sigma_causal`, `Sigma_G_total`,
#'   `Sigma_G_bg`, `Sigma_E`, `Sigma_P_total`), `generator`, `status`
#'   and `diagnostics`.
#' @export
simulate_condped_data <- function(
  n = 1000L,
  m = 4L,
  p = 1000L,
  maf_range = c(0.05, 0.50),
  h2 = rep(0.50, m),
  R_G = NULL,
  R_E = NULL,
  architecture = c(
    "null",
    "candidate_single",
    "candidate_pair",
    "candidate_dense",
    "representative_singleton",
    "representative_pair",
    "representative_full",
    "irreducible_singleton",
    "irreducible_pair",
    "multiple_modules"
  ),
  effect_direction = c("default", "concordant", "discordant", "mixed"),
  locus_pve = 0.01,
  tolerance = 0.10,
  target_loss = 0.05,
  target_representative_set = NULL,
  target_irreducible_modules = NULL,
  correlation = c("block", "independent"),
  structured = TRUE,
  n_groups = 8L,
  fst = 0.05,
  truth_margin = 0.02,
  max_attempts = 200L,
  seed = NULL,
  return_latent = FALSE
) {
  architecture <- match.arg(architecture)
  correlation <- match.arg(correlation)
  effect_direction <- match.arg(effect_direction)
  if (effect_direction != "default" &&
      !architecture %in% c("candidate_pair", "candidate_dense")) {
    .stop_invalid_input(
      paste0(
        "effect_direction is only meaningful for candidate_pair and ",
        "candidate_dense; use \"default\" for architecture \"%s\"."
      ),
      architecture
    )
  }
  .check_count(n, "n", min = 2L)
  .check_count(m, "m", min = 1L)
  .check_count(p, "p", min = 2L)
  .check_count(n_groups, "n_groups", min = 2L)
  .check_count(max_attempts, "max_attempts", min = 1L)
  if (!is.numeric(maf_range) || length(maf_range) != 2L ||
      any(!is.finite(maf_range)) || maf_range[1L] <= 0 ||
      maf_range[1L] >= maf_range[2L] || maf_range[2L] > 0.5) {
    .stop_invalid_input(
      "maf_range must satisfy 0 < maf_range[1] < maf_range[2] <= 0.5."
    )
  }
  if (!is.numeric(h2) || !(length(h2) %in% c(1L, m)) ||
      any(!is.finite(h2)) || any(h2 <= 0) || any(h2 >= 1)) {
    .stop_invalid_input(
      "h2 must have length 1 or m, with all entries in (0, 1)."
    )
  }
  h2 <- rep(h2, length.out = m)
  if (!is.numeric(locus_pve) || length(locus_pve) != 1L ||
      !is.finite(locus_pve) || locus_pve < 0 || locus_pve >= 1) {
    .stop_invalid_input("locus_pve must be a single number in [0, 1).")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0 || tolerance >= 1) {
    .stop_invalid_input("tolerance must be a single number in [0, 1).")
  }
  if (!is.numeric(target_loss) || length(target_loss) != 1L ||
      !is.finite(target_loss) || target_loss < 0 || target_loss >= 1) {
    .stop_invalid_input("target_loss must be a single number in [0, 1).")
  }
  if (!is.numeric(truth_margin) || length(truth_margin) != 1L ||
      !is.finite(truth_margin) || truth_margin < 0) {
    .stop_invalid_input("truth_margin must be a single non-negative number.")
  }
  if (structured && (!is.numeric(fst) || length(fst) != 1L ||
                     !is.finite(fst) || fst <= 0 || fst >= 1)) {
    .stop_invalid_input("fst must be a single number in (0, 1).")
  }
  if (!is.logical(return_latent) || length(return_latent) != 1L ||
      is.na(return_latent)) {
    .stop_invalid_input("return_latent must be TRUE or FALSE.")
  }

  trait_names <- paste0("Trait", seq_len(m))
  ind_names <- paste0("Ind", seq_len(n))
  marker_names <- paste0("M", seq_len(p))

  # ---- covariance structures -------------------------------------------------
  if (is.null(R_G)) {
    R_G <- if (correlation == "block") .block_cor(m) else diag(m)
  }
  if (is.null(R_E)) {
    R_E <- if (correlation == "block") .block_cor(m) else diag(m)
  }
  .check_correlation_matrix(R_G, m, "R_G")
  .check_correlation_matrix(R_E, m, "R_E")
  Sigma_G_total <- .cov_from_h2(h2, R_G)
  Sigma_E <- .cov_from_h2(1 - h2, R_E)
  Sigma_P_total <- Sigma_G_total + Sigma_E
  dimnames(Sigma_G_total) <- dimnames(Sigma_E) <- dimnames(Sigma_P_total) <-
    list(trait_names, trait_names)

  # ---- target sets -------------------------------------------------------------
  targets <- .resolve_scene_targets(
    architecture, m, trait_names,
    target_representative_set, target_irreducible_modules
  )
  structural <- architecture %in% c(
    "representative_singleton", "representative_pair",
    "representative_full", "irreducible_singleton",
    "irreducible_pair", "multiple_modules"
  )
  if (architecture %in% c("representative_singleton",
                          "representative_pair") &&
      abs(target_loss - tolerance) < truth_margin) {
    .stop_invalid_input(
      paste0(
        "|target_loss - tolerance| must be at least truth_margin ",
        "(%.3g); got |%.3g - %.3g| = %.3g."
      ),
      truth_margin, target_loss, tolerance, abs(target_loss - tolerance)
    )
  }

  if (!is.null(seed)) set.seed(seed)

  # Genotypes are drawn once; the acceptance loop only re-draws the
  # causal marker and the effect proposal (both population-level).
  maf <- stats::runif(p, min = maf_range[1L], max = maf_range[2L])
  G <- .sim_genotypes(n, maf, structured = structured,
                      n_groups = n_groups, fst = fst)
  dimnames(G) <- list(ind_names, marker_names)

  # ---- acceptance loop over population truth -----------------------------------
  psd_tol <- 1e-10
  tolerances <- unique(c(tolerance, 0.05, 0.10, 0.20))
  attempt <- 0L
  psd_failures <- 0L
  accepted <- FALSE
  beta <- NULL
  delta <- NA_real_
  causal_index <- NA_integer_
  v_l <- NA_real_
  truth_tabs <- NULL

  repeat {
    attempt <- attempt + 1L
    causal_index <- sample.int(p, 1L)
    v_l <- 2 * maf[causal_index] * (1 - maf[causal_index])
    proposal <- .propose_effects(
      architecture, m, Sigma_P_total, targets, target_loss,
      effect_direction
    )
    beta <- proposal$beta
    delta <- proposal$delta
    A <- which(beta != 0)
    # overall scaling to the target PVE (loss-invariant)
    if (length(A) > 0L && locus_pve > 0) {
      c_scale <- sqrt(locus_pve / (v_l * mean(beta[A]^2)))
      beta <- c_scale * beta
    }
    Sigma_causal <- v_l * tcrossprod(beta)
    Sigma_G_bg <- Sigma_G_total - Sigma_causal
    min_eigen_bg <- .min_eigen_sym(Sigma_G_bg)
    if (min_eigen_bg < -psd_tol) {
      psd_failures <- psd_failures + 1L
      if (attempt >= max_attempts) break
      next
    }
    truth_tabs <- .population_truth(
      beta, Sigma_P_total, tolerance, tolerances,
      marker_names[causal_index], trait_names
    )
    accepted <- !structural ||
      .accept_scene(architecture, truth_tabs, targets, tolerance,
                    truth_margin)
    if (accepted || attempt >= max_attempts) break
  }

  # ---- data generation ----------------------------------------------------------
  bg <- setdiff(seq_len(p), causal_index)
  maf_bg <- maf[bg]
  Z_bg <- sweep(G[, bg, drop = FALSE], 2L, 2 * maf_bg, `-`)
  Z_bg <- sweep(Z_bg, 2L, sqrt(2 * maf_bg * (1 - maf_bg)), `/`)
  K_bg <- .make_grm(Z_bg)
  dimnames(K_bg) <- list(ind_names, ind_names)

  eig_K <- eigen(K_bg, symmetric = TRUE)
  lambda_K <- eig_K$values
  K_rank <- sum(lambda_K > sqrt(.Machine$double.eps) * max(lambda_K))

  S_G <- .symmetric_psd_sqrt(Sigma_G_bg)
  S_E <- .symmetric_psd_sqrt(Sigma_E)
  Z_u <- matrix(stats::rnorm(n * m), nrow = n, ncol = m)
  U <- sweep(eig_K$vectors, 2L, sqrt(pmax(lambda_K, 0)), `*`) %*%
    Z_u %*% S_G
  E <- matrix(stats::rnorm(n * m), nrow = n, ncol = m) %*% S_E
  dimnames(U) <- dimnames(E) <- list(ind_names, trait_names)

  W <- matrix(1, nrow = n, ncol = 1,
              dimnames = list(ind_names, "Intercept"))
  Y <- matrix(G[, causal_index], nrow = n) %*% beta + U + E
  dimnames(Y) <- list(ind_names, trait_names)

  # ---- truth ----------------------------------------------------------------------
  beta_mat <- matrix(beta, nrow = 1L,
                     dimnames = list(marker_names[causal_index],
                                     trait_names))
  sigma2_P <- diag(Sigma_P_total)
  locus_pve_traits <- v_l * beta^2 / sigma2_P
  names(locus_pve_traits) <- trait_names

  truth <- list(
    beta = beta_mat,
    candidate_traits = truth_tabs$candidate_traits,
    subset_table = truth_tabs$subset_table,
    minimum_representative_sets = truth_tabs$minimum_representative_sets,
    irreducible_modules = truth_tabs$irreducible_modules,
    tolerance_path = truth_tabs$tolerance_path,
    locus_pve = locus_pve_traits,
    target_loss = if (architecture %in% c("representative_singleton",
                                          "representative_pair")) {
      target_loss
    } else {
      NA_real_
    },
    delta = delta,
    Sigma_causal = Sigma_causal,
    Sigma_G_total = Sigma_G_total,
    Sigma_G_bg = Sigma_G_bg,
    Sigma_E = Sigma_E,
    Sigma_P_total = Sigma_P_total
  )

  status <- .new_status(
    ok = accepted,
    code = if (accepted) "ok" else "unstable",
    message = if (accepted) "" else sprintf(
      paste0(
        "No generation satisfied the population acceptance rule for ",
        "architecture \"%s\" within %d attempts; the returned object ",
        "contains the last (unaccepted) generation."
      ),
      architecture, max_attempts
    )
  )

  list(
    Y = Y,
    W = W,
    G = G,
    K_bg = K_bg,
    causal_index = causal_index,
    truth = truth,
    latent = if (isTRUE(return_latent)) list(U = U, E = E) else NULL,
    generator = list(
      seed = seed,
      attempts = attempt,
      settings = list(
        n = n, m = m, p = p, maf_range = maf_range, h2 = h2,
        architecture = architecture, locus_pve = locus_pve,
        effect_direction = effect_direction,
        tolerance = tolerance, target_loss = target_loss,
        target_representative_set = targets$representative_set,
        target_irreducible_modules = targets$irreducible_modules,
        correlation = correlation, structured = structured,
        n_groups = n_groups, fst = fst, truth_margin = truth_margin,
        max_attempts = max_attempts, return_latent = return_latent
      )
    ),
    status = status,
    diagnostics = list(
      attempts = attempt,
      psd_failures = psd_failures,
      accepted = accepted,
      mean_diag_K = mean(diag(K_bg)),
      K_eigen_sd = stats::sd(lambda_K),
      K_rank = K_rank,
      min_eigen_Sigma_G_bg = min_eigen_bg,
      pve_realised = if (length(A) > 0L) {
        mean(locus_pve_traits[beta != 0])
      } else {
        0
      },
      truth_margin_achieved = truth_tabs$margin_achieved
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
    .stop_invalid_input("%s must be a single integer >= %d.", name, min)
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

#' Block correlation matrix (within-block 0.6, between-block 0.2)
#'
#' Blocks are consecutive trait pairs (traits 1-2, 3-4, ...) as in the
#' Methods simulation design (traits A,B and C,D form modules).
#' @param m Dimension.
#' @param within Within-block correlation.
#' @param between Between-block correlation.
#' @return An `m x m` correlation matrix.
#' @keywords internal
.block_cor <- function(m, within = 0.6, between = 0.2) {
  R <- matrix(between, m, m)
  block_id <- (seq_len(m) + 1L) %/% 2L
  same <- outer(block_id, block_id, `==`)
  R[same] <- within
  diag(R) <- 1
  R
}

#' Validate a correlation matrix argument
#' @param R Candidate matrix.
#' @param m Expected dimension.
#' @param name Argument name used in error messages.
#' @return Invisible `TRUE`; stops otherwise.
#' @keywords internal
.check_correlation_matrix <- function(R, m, name) {
  if (!is.matrix(R) || !is.numeric(R) || nrow(R) != m || ncol(R) != m) {
    .stop_invalid_input("%s must be a numeric %d x %d matrix.", name, m, m)
  }
  if (anyNA(R) || max(abs(R - t(R))) > 1e-10) {
    .stop_invalid_input(
      "%s must be symmetric and free of missing values.", name
    )
  }
  if (max(abs(diag(R) - 1)) > 1e-10) {
    .stop_invalid_input("%s must have unit diagonal.", name)
  }
  if (.min_eigen_sym(R) <= 0) {
    .stop_invalid_input("%s must be positive definite.", name)
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
#' Legacy helper still used by the `fit_mt_null()` compatibility
#' fields (`compute_gamma`). Scheduled for removal when the null-model
#' interface is migrated.
#'
#' @param Sigma_P `m x m` phenotypic covariance matrix.
#' @param trait_names Character vector of trait names.
#' @return List with `gamma`, `C` and `conditional_variance`.
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

#' Resolve target sets for the structural scenarios
#'
#' @param architecture Architecture label.
#' @param m Number of traits.
#' @param trait_names Trait names.
#' @param target_representative_set,target_irreducible_modules User
#'   overrides or `NULL`.
#' @return List with `representative_set` (character or NULL) and
#'   `irreducible_modules` (list of character vectors or NULL).
#' @keywords internal
.resolve_scene_targets <- function(architecture, m, trait_names,
                                   target_representative_set,
                                   target_irreducible_modules) {
  rep_set <- target_representative_set
  mod_sets <- target_irreducible_modules
  if (architecture == "representative_singleton" && is.null(rep_set)) {
    rep_set <- trait_names[1L]
  }
  if (architecture == "representative_pair" && is.null(rep_set)) {
    if (m < 2L) {
      .stop_invalid_input("representative_pair requires m >= 2.")
    }
    rep_set <- trait_names[1:2]
  }
  if (architecture == "irreducible_singleton" && is.null(mod_sets)) {
    mod_sets <- list(trait_names[1L])
  }
  if (architecture == "irreducible_pair" && is.null(mod_sets)) {
    if (m < 2L) {
      .stop_invalid_input("irreducible_pair requires m >= 2.")
    }
    mod_sets <- list(trait_names[1:2])
  }
  if (architecture == "multiple_modules" && is.null(mod_sets)) {
    if (m < 4L) {
      .stop_invalid_input("multiple_modules requires m >= 4.")
    }
    mod_sets <- list(trait_names[1L], trait_names[3:4])
  }
  if (!is.null(rep_set)) {
    rep_set <- .normalize_trait_set(rep_set, trait_names,
                                    "target_representative_set")
    if (length(rep_set) == 0L) {
      .stop_invalid_input("target_representative_set must be non-empty.")
    }
  }
  if (!is.null(mod_sets)) {
    if (!is.list(mod_sets)) mod_sets <- list(mod_sets)
    mod_sets <- lapply(mod_sets, .normalize_trait_set, trait_names,
                       "target_irreducible_modules")
    if (any(lengths(mod_sets) == 0L)) {
      .stop_invalid_input("target modules must be non-empty.")
    }
  }
  list(representative_set = rep_set, irreducible_modules = mod_sets)
}

#' Propose a population effect vector for one attempt
#'
#' Representative scenarios use the closed-form construction of
#' Methods eqs. (39)-(41); other structural scenarios draw from
#' proposal families that the acceptance rule filters.
#'
#' `effect_direction` shapes only the candidate_pair/candidate_dense
#' templates: `concordant` keeps every effect in the same direction,
#' `discordant` flips all non-leading effects, `mixed` alternates the
#' signs (identical to `discordant` for a pair).
#'
#' @param architecture Architecture label.
#' @param m Number of traits.
#' @param Sigma_P Population phenotypic covariance (m x m, dimnamed).
#' @param targets Resolved target sets.
#' @param target_loss Closed-form target loss r.
#' @param effect_direction One of `"default"`, `"concordant"`,
#'   `"discordant"`, `"mixed"`.
#' @return List with `beta` (length m) and `delta`.
#' @keywords internal
.propose_effects <- function(architecture, m, Sigma_P, targets,
                             target_loss,
                             effect_direction = "default") {
  trait_names <- colnames(Sigma_P)
  beta <- numeric(m)
  delta <- NA_real_
  switch(
    architecture,
    null = {},
    candidate_single = {
      beta[1L] <- 1
    },
    candidate_pair = {
      if (m < 2L) .stop_invalid_input("candidate_pair requires m >= 2.")
      beta[1:2] <- switch(effect_direction,
        default = c(1, 1),
        concordant = c(1, 0.8),
        discordant = c(1, -0.8),
        mixed = c(1, -0.8)     # alternating signs coincide for a pair
      )
    },
    candidate_dense = {
      base <- 1 - 0.2 * seq(0, length.out = m)
      beta <- switch(effect_direction,
        default = base,
        concordant = base,
        discordant = base * c(1, rep(-1, m - 1L)),
        mixed = base * rep(c(1, -1), length.out = m)
      )
    },
    representative_singleton = ,
    representative_pair = {
      S <- match(targets$representative_set, trait_names)
      R <- setdiff(seq_len(m), S)
      if (length(R) == 0L) {
        .stop_invalid_input(
          "the closed-form representative construction requires a non-empty complement."
        )
      }
      beta_S <- stats::rnorm(length(S))
      Sigma_SS <- Sigma_P[S, S, drop = FALSE]
      inv_SS <- .safe_inverse(Sigma_SS)
      Gamma <- inv_SS$inverse %*% Sigma_P[S, R, drop = FALSE]
      Omega <- Sigma_P[R, R, drop = FALSE] -
        crossprod(Sigma_P[S, R, drop = FALSE],
                  inv_SS$inverse %*% Sigma_P[S, R, drop = FALSE])
      Omega <- (Omega + t(Omega)) / 2
      w <- stats::rnorm(length(R))
      d <- drop(Omega %*% w)
      d <- d / sqrt(sum(w * d))          # d' Omega^{-1} d = 1
      q_S <- sum(beta_S * (inv_SS$inverse %*% beta_S))
      delta <- sqrt(target_loss * q_S / (1 - target_loss))
      beta[S] <- beta_S
      beta[R] <- drop(crossprod(Gamma, beta_S)) + delta * d
    },
    representative_full = {
      beta <- sign(stats::rnorm(m)) * (1 + 0.3 * abs(stats::rnorm(m)))
    },
    irreducible_singleton = {
      # target module trait carries the bulk; the rest stay small so
      # every set containing the target is feasible
      beta[1L] <- 1
      if (m > 1L) {
        beta[-1L] <- stats::runif(m - 1L, 0.05, 0.20) *
          sample(c(-1, 1), m - 1L, replace = TRUE)
      }
    },
    irreducible_pair = {
      # Traits 1-2 share an effect chunk. Each is pinned to the
      # best-linear prediction from the remaining traits plus a small
      # idiosyncratic residual (Gauss-Seidel passes), so deleting ONE
      # of them stays feasible; their joint residual against traits
      # 3..m (between-block correlation only 0.2) stays large, so
      # deleting BOTH is infeasible.
      beta <- numeric(m)
      beta[3:m] <- stats::rnorm(m - 2L, sd = 0.30)
      e <- stats::rnorm(2L, sd = 0.15)
      for (iter in 1:3) {
        for (j in 1:2) {
          others <- setdiff(seq_len(m), j)
          inv <- .safe_inverse(Sigma_P[others, others, drop = FALSE])
          g <- drop(inv$inverse %*% Sigma_P[others, j])
          beta[j] <- sum(g * beta[others]) + e[j]
        }
      }
    },
    multiple_modules = {
      # Module {1}: stand-alone effect. Module {3,4}: joint effect with
      # nearly EQUAL coefficients, so each of traits 3/4 is explainable
      # by the other through the within-block correlation (small
      # residual), while deleting both loses the whole chunk. Trait 2
      # (and any traits beyond 4) carry no effect and stay out of the
      # candidate set.
      beta <- numeric(m)
      beta[1L] <- 1
      beta[3L] <- stats::runif(1L, 0.58, 0.80) +
        stats::rnorm(1L, sd = 0.015)
      if (m >= 4L) {
        beta[4L] <- beta[3L] + stats::rnorm(1L, sd = 0.015)
      }
    }
  )
  names(beta) <- trait_names
  list(beta = beta, delta = delta)
}

#' Population structural truth for one causal variant
#'
#' Evaluates every subset of the candidate set with the unique Stage-1
#' loss path on population parameters and extracts the minimum
#' representative sets and irreducible modules per tolerance with the
#' Stage-2 extractors. This is the same deterministic mathematics the
#' estimation path uses.
#'
#' @param beta Named effect vector (length m).
#' @param Sigma_P Population phenotypic covariance.
#' @param tolerance Primary tolerance.
#' @param tolerances All tolerances to extract at.
#' @param marker_id Causal marker id.
#' @param trait_names Trait names.
#' @return List with `candidate_traits`, `subset_table`,
#'   `minimum_representative_sets`, `irreducible_modules`,
#'   `tolerance_path`, `margin_achieved`.
#' @keywords internal
.population_truth <- function(beta, Sigma_P, tolerance, tolerances,
                              marker_id, trait_names) {
  A <- trait_names[beta != 0]
  empty <- list(
    candidate_traits = character(),
    subset_table = data.frame(
      marker_id = character(), set_id = integer(),
      representing_set = I(list()), representing_key = character(),
      complement_set = I(list()), complement_key = character(),
      set_size = integer(), conditional_effect = I(list()),
      full_qform = numeric(),
      subset_qform = numeric(), residual_qform = numeric(),
      representation_loss = numeric(), feasible_primary = logical(),
      status = character(), stringsAsFactors = FALSE
    ),
    minimum_representative_sets = data.frame(
      marker_id = character(), tolerance = numeric(),
      trait_set = I(list()), trait_key = character(),
      set_size = integer(), representation_loss = numeric(),
      n_tied_solutions = integer(), stringsAsFactors = FALSE
    ),
    irreducible_modules = data.frame(
      marker_id = character(), tolerance = numeric(),
      trait_set = I(list()), trait_key = character(),
      module_size = integer(), complement_set = I(list()),
      complement_key = character(), complement_loss = numeric(),
      stringsAsFactors = FALSE
    ),
    tolerance_path = data.frame(
      marker_id = character(), tolerance = numeric(),
      minimum_set_size = integer(), n_minimum_sets = integer(),
      n_irreducible_modules = integer(), stringsAsFactors = FALSE
    ),
    margin_achieved = NA_real_
  )
  k <- length(A)
  if (k == 0L) return(empty)

  Sigma_AA <- Sigma_P[A, A, drop = FALSE]
  beta_A <- beta[A]
  enum <- .enumerate_subsets(A, "all")
  rows <- vector("list", length(enum$sets))
  for (si in seq_along(enum$sets)) {
    out <- .compute_subset_loss(beta_A, Sigma_AA, enum$sets[[si]],
                                trait_names = A)
    rows[[si]] <- data.frame(
      marker_id = marker_id, set_id = si,
      representing_set = I(list(out$representing_set)),
      representing_key = out$representing_key,
      complement_set = I(list(out$complement_set)),
      complement_key = out$complement_key,
      set_size = length(out$representing_set),
      conditional_effect = I(list(out$eta)),
      full_qform = out$full_qform,
      subset_qform = out$subset_qform,
      residual_qform = out$residual_qform,
      representation_loss = out$representation_loss,
      feasible_primary = !is.na(out$representation_loss) &&
        out$representation_loss <= tolerance + 1e-10,
      status = out$status$code,
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, rows)

  reps_list <- mods_list <- path_list <- list()
  for (tol in tolerances) {
    r <- .extract_minimum_representative_sets(tab, tol)
    mo <- .extract_irreducible_modules(tab, tol)
    if (nrow(r) > 0L) {
      reps_list[[length(reps_list) + 1L]] <-
        cbind(marker_id = marker_id, r[, setdiff(names(r), "status")],
              stringsAsFactors = FALSE)
    }
    if (nrow(mo) > 0L) {
      mods_list[[length(mods_list) + 1L]] <-
        cbind(marker_id = marker_id, mo[, setdiff(names(mo), "status")],
              stringsAsFactors = FALSE)
    }
    path_list[[length(path_list) + 1L]] <- data.frame(
      marker_id = marker_id, tolerance = tol,
      minimum_set_size = if (nrow(r) > 0L) r$set_size[1L] else NA_integer_,
      n_minimum_sets = nrow(r), n_irreducible_modules = nrow(mo),
      stringsAsFactors = FALSE
    )
  }
  bind <- function(lst, proto) {
    if (length(lst) == 0L) return(proto)
    out <- do.call(rbind, lst)
    rownames(out) <- NULL
    out
  }
  losses <- tab$representation_loss
  list(
    candidate_traits = A,
    subset_table = tab,
    minimum_representative_sets = bind(reps_list,
                                       empty$minimum_representative_sets),
    irreducible_modules = bind(mods_list, empty$irreducible_modules),
    tolerance_path = bind(path_list, empty$tolerance_path),
    margin_achieved = if (all(is.finite(losses))) {
      min(abs(losses - tolerance))
    } else {
      NA_real_
    }
  )
}

#' Population acceptance rule for the structural scenarios
#'
#' Reads only the population truth tables; never any sample quantity.
#'
#' @param architecture Architecture label.
#' @param tabs Output of [`.population_truth()`].
#' @param targets Resolved target sets.
#' @param tolerance Primary tolerance.
#' @param margin Required distance from the tolerance boundary.
#' @return `TRUE` when the generation satisfies its frozen acceptance
#'   conditions.
#' @keywords internal
.accept_scene <- function(architecture, tabs, targets, tolerance,
                          margin) {
  tab <- tabs$subset_table
  if (nrow(tab) == 0L) return(FALSE)
  losses <- tab$representation_loss
  if (any(!is.finite(losses))) return(FALSE)
  if (margin > 0 && any(abs(losses - tolerance) < margin)) {
    return(FALSE)
  }
  loss_of <- function(key) {
    i <- match(key, tab$representing_key)
    if (is.na(i)) return(NA_real_)
    tab$representation_loss[i]
  }
  reps <- tabs$minimum_representative_sets
  reps <- reps[reps$tolerance == tolerance, ]
  mods <- tabs$irreducible_modules
  mods <- mods[mods$tolerance == tolerance, ]
  switch(
    architecture,
    representative_singleton = ,
    representative_pair = {
      S <- targets$representative_set
      S_key <- .trait_set_key(S)
      # target feasible and of minimum size; no smaller feasible set
      is_min <- nrow(reps) > 0L &&
        reps$set_size[1L] == length(S) && S_key %in% reps$trait_key
      if (!is_min) return(FALSE)
      if (architecture == "representative_pair") {
        singles <- tab[tab$set_size == 1L, ]
        if (any(singles$representation_loss <= tolerance)) return(FALSE)
      }
      TRUE
    },
    representative_full = {
      proper <- tab[tab$set_size < max(tab$set_size), ]
      all(proper$representation_loss > tolerance)
    },
    irreducible_singleton = ,
    irreducible_pair = {
      T <- targets$irreducible_modules[[1L]]
      comp_key <- .trait_set_key(setdiff(tabs$candidate_traits, T))
      # deleting the whole module is infeasible
      if (!is.finite(loss_of(comp_key)) ||
          loss_of(comp_key) <= tolerance) {
        return(FALSE)
      }
      # deleting any proper subset stays feasible (minimality)
      for (sz in seq_len(length(T) - 1L)) {
        for (sub in utils::combn(T, sz, simplify = FALSE)) {
          ck <- .trait_set_key(setdiff(tabs$candidate_traits, sub))
          if (!is.finite(loss_of(ck)) || loss_of(ck) > tolerance) {
            return(FALSE)
          }
        }
      }
      # and the module is actually extracted at the primary tolerance
      .trait_set_key(T) %in% mods$trait_key
    },
    multiple_modules = {
      want <- sort(vapply(targets$irreducible_modules,
                          .trait_set_key, character(1)))
      identical(sort(mods$trait_key), want)
    },
    TRUE
  )
}
