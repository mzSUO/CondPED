#' Fit the multi-trait REML null model
#'
#' Fits the multi-trait null model \eqn{y = X_0 b + u + \varepsilon} with
#' covariance \eqn{V = \Sigma_G \otimes K + \Sigma_E \otimes I_n} by
#' restricted maximum likelihood, following the interface contract
#' (section 4.2) and Methods equations (2)-(6). Returns the rotation
#' objects and the global conditional contrasts reused by the whole
#' pipeline.
#'
#' @details
#' The genomic relationship matrix is eigen-decomposed exactly once,
#' \eqn{K = U \Lambda U^\top}, and the data are rotated:
#' \eqn{\widetilde Y = U^\top Y}, \eqn{\widetilde W = U^\top W}. In the
#' rotated space every coordinate \eqn{j} has the \eqn{m \times m}
#' covariance \eqn{V_j = \lambda_j \Sigma_G + \Sigma_E}, so the full
#' \eqn{nm \times nm} matrices \eqn{V} and \eqn{P} are never constructed.
#' Both covariance matrices use the log-Cholesky parameterisation
#' \eqn{\Sigma = L L^\top} with \eqn{L_{ii} = \exp(\theta_{ii})} (floored
#' at `chol_floor`), which keeps them positive definite throughout the
#' optimisation. The restricted log-likelihood includes the fixed-effect
#' correction term \eqn{\log|X_0^\top V^{-1} X_0|} (contract section 2.3),
#' and the optimiser receives an analytic gradient of the rotated
#' objective.
#'
#' Starting values split the sample phenotypic covariance by `h2_start`
#' (projected to the positive-definite space); additional starts randomly
#' re-split it. The solution with the largest restricted log-likelihood
#' among the converged, legally estimated runs is retained.
#'
#' When `K` is (numerically) the identity or its eigenvalue dispersion is
#' very low, the two variance components are only weakly identifiable;
#' this is reported through `diagnostics$identifiable_warning` and a status
#' warning instead of being silently ignored (contract section 2.1).
#'
#' @param Y Numeric `n x m` phenotype matrix (rows are individuals).
#' @param W Optional numeric `n x q` fixed-effect design matrix; when
#'   `NULL` an intercept-only design is used.
#' @param K Numeric `n x n` genomic relationship matrix, scaled to average
#'   diagonal 1.
#' @param add_intercept Logical; add an intercept column to `W` when it is
#'   not already present.
#' @param optimizer Optimizer used for the restricted log-likelihood;
#'   `"BFGS"` (via [stats::optim()]) or `"nlminb"` (via
#'   [stats::nlminb()]).
#' @param n_starts Integer; number of starting points (the first is the
#'   deterministic `h2_start` split, the rest are random re-splits; the
#'   random draws follow the current RNG stream, so `set.seed()` before the
#'   call gives full reproducibility).
#' @param h2_start Numeric in (0, 1); heritability used to split the sample
#'   phenotypic covariance into starting values.
#' @param start_values Optional explicit first start: either a numeric
#'   log-Cholesky parameter vector of length `m * (m + 1)` (first
#'   `m * (m + 1) / 2` for `Sigma_G`, then the same for `Sigma_E`), or a
#'   list with components `Sigma_G` and `Sigma_E`.
#' @param eig_tol Tolerance used when flooring the eigenvalues of `K` and
#'   computing its numerical rank.
#' @param inverse_tol Relative tolerance forwarded to `.safe_inverse()`.
#' @param chol_floor Lower bound for the log-Cholesky diagonal (standard
#'   deviation scale).
#' @param control List of optimizer control settings; `maxit` (default 500)
#'   and `reltol` (default 1e-8) are honoured.
#' @param compute_gamma Logical; compute the global conditional projection
#'   coefficients from \eqn{\widehat\Sigma_P = \widehat\Sigma_G +
#'   \widehat\Sigma_E} (computed once, fixed for the whole pipeline).
#' @param return_rotation Logical; return the rotation object (eigenvectors,
#'   eigenvalues, rotated data, per-coordinate inverse covariances and the
#'   GLS components).
#' @param verbose Logical; print per-start optimisation progress.
#'
#' @return An object of class `"condped_mt_null"`: a list with components
#'   \describe{
#'     \item{Sigma_G, Sigma_E, Sigma_P_ref}{Estimated `m x m` covariance
#'       matrices. `Sigma_P_ref = Sigma_G + Sigma_E` is the fixed
#'       reference phenotypic covariance shared by all loci and subset
#'       analyses.}
#'     \item{fixed_effects}{`q x m` matrix of GLS fixed-effect estimates.}
#'     \item{gamma}{`m x (m - 1)` matrix of projection coefficients
#'       (row `i` is \eqn{\gamma_{i,-i}}); `NULL` when
#'       `compute_gamma = FALSE`.}
#'     \item{contrasts}{`m x m` matrix whose column `i` is the contrast
#'       vector \eqn{c_i}.}
#'     \item{conditional_variance}{Length `m` vector of Schur complements.}
#'     \item{logLik}{Maximised restricted log-likelihood.}
#'     \item{npar}{Number of variance parameters, `m * (m + 1)`.}
#'     \item{optimizer}{Optimizer actually used.}
#'     \item{trait_names}{Trait names (from `colnames(Y)` or generated).}
#'     \item{individual_ids}{Individual ids (from `rownames(Y)` or
#'       generated).}
#'     \item{convergence}{List with `code`, `message`, `gradient_norm`,
#'       `iterations`, `selected_start` and `starts` (a data frame recording
#'       the objective value, convergence code, legality and selection flag
#'       of every starting point).}
#'     \item{rotation}{List with `U`, `lambda`, `Y_tilde`, `W_tilde`,
#'       `Vinv` (an `m x m x n` array of per-coordinate inverse
#'       covariances), `logdet_Vj`, `XtVinvX`, `XtVinvX_inv`,
#'       `beta_fixed` (vectorised fixed effects) and `residual_tilde`;
#'       `NULL` when `return_rotation = FALSE`.}
#'     \item{diagnostics}{List with `K_rank`, `K_condition`, `K_eigen_sd`,
#'       `min_eigen_Sigma_G`, `min_eigen_Sigma_E`, `min_eigen_Sigma_P_ref`,
#'       `condition_Sigma_P_ref` and `identifiable_warning`.}
#'     \item{status}{Standard CondPED status list.}
#'   }
#' @export
fit_mt_null <- function(
  Y,
  W = NULL,
  K,
  add_intercept = TRUE,
  optimizer = c("BFGS", "nlminb"),
  n_starts = 5L,
  h2_start = 0.50,
  start_values = NULL,
  eig_tol = 1e-8,
  inverse_tol = sqrt(.Machine$double.eps),
  chol_floor = 1e-8,
  control = list(maxit = 500L, reltol = 1e-8),
  compute_gamma = TRUE,
  return_rotation = TRUE,
  verbose = FALSE
) {
  optimizer <- match.arg(optimizer)
  .validate_dimensions(Y = Y, K = K)
  Y <- as.matrix(Y)
  n <- nrow(Y)
  m <- ncol(Y)
  trait_names <- colnames(Y)
  if (is.null(trait_names)) trait_names <- paste0("Trait", seq_len(m))

  if (!is.numeric(n_starts) || length(n_starts) != 1L ||
      !is.finite(n_starts) || n_starts != round(n_starts) || n_starts < 1) {
    .stop_invalid_input("n_starts must be a single positive integer.")
  }
  n_starts <- as.integer(n_starts)
  if (!is.numeric(h2_start) || length(h2_start) != 1L ||
      !is.finite(h2_start) || h2_start <= 0 || h2_start >= 1) {
    .stop_invalid_input("h2_start must be a single number in (0, 1).")
  }
  maxit <- if (!is.null(control$maxit)) control$maxit else 500L
  reltol <- if (!is.null(control$reltol)) control$reltol else 1e-8

  # ---- fixed-effect design ---------------------------------------------------
  W <- .prepare_design(W, n, add_intercept)
  q <- ncol(W)
  .validate_dimensions(Y = Y, W = W, K = K)

  # ---- K: exactly one eigen decomposition -------------------------------------
  if (!isSymmetric(unname(K), tol = 1e-8)) {
    .stop_invalid_input("K must be symmetric (max |K - t(K)| exceeds 1e-8).")
  }
  K <- (K + t(K)) / 2
  eig_K <- eigen(K, symmetric = TRUE)
  if (min(eig_K$values) < -eig_tol) {
    .stop_invalid_input(
      "K is not positive semi-definite (smallest eigenvalue < -eig_tol)."
    )
  }
  lambda <- pmax(eig_K$values, 0)
  U <- eig_K$vectors
  K_rank <- sum(lambda > eig_tol * max(lambda))
  lambda_pos <- lambda[lambda > eig_tol * max(lambda)]
  K_condition <- max(lambda) / min(lambda_pos)

  Y_tilde <- crossprod(U, Y)
  W_tilde <- crossprod(U, W)

  # ---- starting values ---------------------------------------------------------
  starts <- .reml_starts(
    Y = Y, m = m, n_starts = n_starts, h2_start = h2_start,
    start_values = start_values, chol_floor = chol_floor
  )

  objective <- .make_reml_objective(
    m = m, lambda = lambda, Y_tilde = Y_tilde, W_tilde = W_tilde,
    chol_floor = chol_floor
  )

  # ---- multi-start optimisation -------------------------------------------------
  runs <- vector("list", length(starts))
  for (s in seq_along(starts)) {
    theta0 <- starts[[s]]
    if (optimizer == "BFGS") {
      fit <- stats::optim(
        theta0, fn = objective$fn, gr = objective$gr, method = "BFGS",
        control = list(fnscale = -1, maxit = maxit, reltol = reltol,
                       trace = if (isTRUE(verbose)) 1L else 0L))
      runs[[s]] <- list(
        par = fit$par, logLik = fit$value, code = fit$convergence,
        message = if (is.null(fit$message)) "" else fit$message,
        iterations = unname(fit$counts[["function"]])
      )
    } else {
      fit <- stats::nlminb(
        theta0, objective = objective$fn, gradient = objective$gr,
        control = list(iter.max = maxit, rel.tol = reltol,
                       trace = if (isTRUE(verbose)) 1L else 0L))
      runs[[s]] <- list(
        par = fit$par, logLik = -fit$objective, code = fit$convergence,
        message = fit$message, iterations = fit$iterations
      )
    }
    if (isTRUE(verbose)) {
      cat(sprintf("start %d: logLik = %.6f (code %d)\n",
                  s, runs[[s]]$logLik, runs[[s]]$code))
    }
  }

  # ---- select the best legal run --------------------------------------------------
  # Legality: a run is legal only if it produced a genuine objective value.
  # The sentinel -1e10 marks numerically illegal parameters; an optimiser
  # may still report convergence inside that flat illegal region, so the
  # sentinel check takes precedence over the convergence code.
  bad_ll <- -1e10
  legal <- vapply(runs, function(r) {
    is.finite(r$logLik) && r$logLik > bad_ll / 2
  }, logical(1))
  converged <- legal & vapply(runs, function(r) r$code == 0, logical(1))
  pool <- if (any(converged)) which(converged) else which(legal)
  all_failed <- length(pool) == 0L
  if (all_failed) {
    # No legal run at all: fall back to the first start purely so that a
    # complete, honestly flagged object can be returned (status$ok = FALSE).
    pool <- 1L
  }
  selected <- pool[which.max(vapply(runs[pool], function(r) r$logLik,
                                    numeric(1)))]
  best <- runs[[selected]]
  theta_hat <- best$par
  starts_table <- data.frame(
    start = seq_along(runs),
    logLik = vapply(runs, function(r) r$logLik, numeric(1)),
    code = vapply(runs, function(r) r$code, numeric(1)),
    converged = converged,
    legal = legal,
    selected = seq_along(runs) == selected
  )

  # ---- rebuild all components at the optimum --------------------------------------
  ncp <- m * (m + 1L) / 2L
  Sigma_G <- .logchol_unpack(theta_hat[seq_len(ncp)], m, chol_floor)
  Sigma_E <- .logchol_unpack(theta_hat[ncp + seq_len(ncp)], m, chol_floor)
  Sigma_P <- Sigma_G + Sigma_E
  dimnames(Sigma_G) <- dimnames(Sigma_E) <- dimnames(Sigma_P) <-
    list(trait_names, trait_names)

  gls <- .gls_components_rotated(
    Sigma_G = Sigma_G, Sigma_E = Sigma_E, lambda = lambda,
    Y_tilde = Y_tilde, W_tilde = W_tilde, inverse_tol = inverse_tol,
    keep_vinv = return_rotation
  )
  beta_vec <- gls$beta_fixed
  # beta_vec is indexed as (covariate a, trait i) -> (a - 1) * m + i,
  # matching X0 = I_m %x% W; reshape to the q x m fixed-effect matrix.
  fixed_effects <- t(matrix(beta_vec, nrow = m, ncol = q,
                            dimnames = list(trait_names, colnames(W))))
  residual_tilde <- Y_tilde - W_tilde %*% fixed_effects

  # Global conditional contrasts from the fitted Sigma_P (computed once,
  # fixed for the whole pipeline).
  if (isTRUE(compute_gamma)) {
    if (m > 1L) {
      ctr <- .contrasts_from_cov(Sigma_P, trait_names)
    } else {
      ctr <- list(
        gamma = matrix(numeric(0), nrow = 1L, ncol = 0L),
        C = matrix(1, 1, 1, dimnames = list(trait_names, trait_names)),
        conditional_variance = Sigma_P[1, 1]
      )
    }
    gamma <- ctr$gamma
    contrasts <- ctr$C
    conditional_variance <- ctr$conditional_variance
  } else {
    gamma <- NULL
    contrasts <- NULL
    conditional_variance <- NULL
  }

  # ---- diagnostics and status ------------------------------------------------------
  K_eigen_sd <- stats::sd(lambda)
  # Weak identifiability: K numerically the identity or near-zero
  # eigenvalue dispersion (contract section 2.1). The threshold is on the
  # sd of the eigenvalues of K (mean diagonal 1 scale); 0.05 is far below
  # any mildly structured or background-marker GRM.
  identifiable_warning <- K_eigen_sd < 0.05
  min_eigen_G <- .min_eigen_sym(Sigma_G)
  min_eigen_E <- .min_eigen_sym(Sigma_E)
  condition_P <- {
    ev <- eigen(Sigma_P, symmetric = TRUE, only.values = TRUE)$values
    max(ev) / min(ev)
  }

  gradient_norm <- sqrt(sum(objective$gr(theta_hat)^2))

  ok <- !all_failed && isTRUE(converged[selected])
  warnings <- character()
  if (all_failed) {
    warnings <- c(warnings, paste(
      "All optimisation starts failed to produce a legal log-likelihood;",
      "the returned estimates come from the first start and must not be",
      "interpreted as a fitted model."
    ))
  } else if (!ok) {
    warnings <- c(warnings, sprintf(
      "No start fully converged; returning the best legal run (start %d, code %d).",
      selected, best$code))
  }
  if (identifiable_warning) {
    warnings <- c(warnings, paste(
      "K has near-zero eigenvalue dispersion (sd(eigen(K)) < 0.05):",
      "Sigma_G and Sigma_E are only weakly identifiable. Interpret",
      "variance-component estimates with caution."
    ))
  }
  status <- .new_status(
    ok = ok,
    code = if (ok) "ok" else "non_convergence",
    message = if (ok) {
      ""
    } else if (all_failed) {
      "All optimisation starts failed to produce a legal log-likelihood."
    } else {
      "Optimiser did not report full convergence."
    },
    warnings = warnings
  )

  structure(
    list(
      Sigma_G = Sigma_G,
      Sigma_E = Sigma_E,
      Sigma_P_ref = Sigma_P,
      fixed_effects = fixed_effects,
      gamma = gamma,
      contrasts = contrasts,
      conditional_variance = conditional_variance,
      logLik = best$logLik,
      npar = m * (m + 1L),
      optimizer = optimizer,
      trait_names = trait_names,
      individual_ids = if (!is.null(rownames(Y))) {
        rownames(Y)
      } else {
        paste0("Ind", seq_len(n))
      },
      convergence = list(
        code = best$code,
        message = best$message,
        gradient_norm = gradient_norm,
        iterations = best$iterations,
        selected_start = selected,
        starts = starts_table
      ),
      rotation = if (isTRUE(return_rotation)) {
        list(
          U = U,
          lambda = lambda,
          Y_tilde = Y_tilde,
          W_tilde = W_tilde,
          Vinv = gls$Vinv,
          logdet_Vj = gls$logdet_Vj,
          XtVinvX = gls$XtVinvX,
          XtVinvX_inv = gls$XtVinvX_inv,
          beta_fixed = beta_vec,
          residual_tilde = residual_tilde
        )
      },
      diagnostics = list(
        K_rank = K_rank,
        K_condition = K_condition,
        K_eigen_sd = K_eigen_sd,
        min_eigen_Sigma_G = min_eigen_G,
        min_eigen_Sigma_E = min_eigen_E,
        min_eigen_Sigma_P_ref = min(ev),
        condition_Sigma_P_ref = condition_P,
        identifiable_warning = identifiable_warning
      ),
      status = status
    ),
    class = "condped_mt_null"
  )
}

# ---- internal helpers -----------------------------------------------------------

#' Prepare the fixed-effect design matrix
#' @param W User-supplied design or `NULL`.
#' @param n Number of individuals.
#' @param add_intercept Logical; add an intercept when absent.
#' @return An `n x q` numeric matrix with column names.
#' @keywords internal
.prepare_design <- function(W, n, add_intercept) {
  if (is.null(W)) {
    return(matrix(1, nrow = n, ncol = 1,
                  dimnames = list(NULL, "Intercept")))
  }
  W <- as.matrix(W)
  if (!is.numeric(W) || nrow(W) != n) {
    .stop_invalid_input("W must be a numeric matrix with n rows.")
  }
  if (is.null(colnames(W))) colnames(W) <- paste0("W", seq_len(ncol(W)))
  has_intercept <- any(apply(W, 2L, function(x) all(x == x[1L])))
  if (add_intercept && !has_intercept) {
    W <- cbind(Intercept = 1, W)
  }
  W
}

#' Pack a covariance matrix into log-Cholesky parameters
#'
#' @param Sigma Positive definite `m x m` matrix.
#' @return Numeric vector of length `m * (m + 1) / 2`: the log of the
#'   Cholesky diagonal, followed by the sub-diagonal entries of the lower
#'   Cholesky factor in column-major order.
#' @keywords internal
.logchol_pack <- function(Sigma) {
  L <- t(chol(Sigma))
  c(log(diag(L)), L[lower.tri(L)])
}

#' Unpack log-Cholesky parameters into a covariance matrix
#'
#' @param theta Numeric vector of length `m * (m + 1) / 2` (see
#'   [`.logchol_pack()`]).
#' @param m Matrix dimension.
#' @param chol_floor Lower bound for the diagonal of the Cholesky factor
#'   (standard-deviation scale), keeping the result positive definite.
#' @return The `m x m` covariance matrix \eqn{\Sigma = L L^\top}.
#' @keywords internal
.logchol_unpack <- function(theta, m, chol_floor = 1e-8) {
  .logchol_unpack_factor(theta, m, chol_floor)$Sigma
}

#' Unpack log-Cholesky parameters, returning the factor as well
#' @inheritParams .logchol_unpack
#' @return List with `L` (lower Cholesky factor, diagonal floored at
#'   `chol_floor`), `Sigma` and `diag_scale` (the chain-rule factor
#'   \eqn{d L_{aa} / d \theta_a}).
#' @keywords internal
.logchol_unpack_factor <- function(theta, m, chol_floor = 1e-8) {
  ncp <- m * (m + 1L) / 2L
  if (length(theta) != ncp) {
    stop(sprintf("theta must have length %d for m = %d.", ncp, m),
         call. = FALSE)
  }
  L <- matrix(0, m, m)
  d_raw <- exp(theta[seq_len(m)])
  diag(L) <- pmax(d_raw, chol_floor)
  if (m > 1L) L[lower.tri(L)] <- theta[(m + 1L):ncp]
  list(
    L = L,
    Sigma = tcrossprod(L),
    diag_scale = ifelse(d_raw > chol_floor, d_raw, 0)
  )
}

#' Rotated restricted (or full) log-likelihood of the multi-trait null model
#'
#' Evaluates the (restricted) log-likelihood of
#' \eqn{y = X_0 b + u + \varepsilon},
#' \eqn{V = \Sigma_G \otimes K + \Sigma_E \otimes I_n}, in the eigen-rotated
#' space where coordinate \eqn{j} has covariance
#' \eqn{V_j = \lambda_j \Sigma_G + \Sigma_E}. The REML version includes the
#' fixed-effect correction \eqn{\log|X_0^\top V^{-1} X_0|}. The full
#' \eqn{nm \times nm} matrices are never constructed.
#'
#' @param theta Log-Cholesky parameters: first `m(m+1)/2` for `Sigma_G`,
#'   then the same number for `Sigma_E`.
#' @param m Number of traits.
#' @param lambda Eigenvalues of `K` (length n, non-negative).
#' @param Y_tilde Rotated phenotypes `U'Y` (`n x m`).
#' @param W_tilde Rotated design `U'W` (`n x q`).
#' @param chol_floor Log-Cholesky diagonal floor.
#' @param reml Logical; `TRUE` (default) computes the restricted
#'   log-likelihood, `FALSE` the full (profile) ML log-likelihood.
#' @return Scalar log-likelihood (higher is better); `-1e10` for
#'   numerically illegal parameters.
#' @keywords internal
.reml_objective_rotated <- function(theta, m, lambda, Y_tilde, W_tilde,
                                    chol_floor = 1e-8, reml = TRUE) {
  .reml_value_grad(theta, m, lambda, Y_tilde, W_tilde, chol_floor, reml,
                   need_grad = FALSE)$value
}

#' Rotated REML/ML objective with analytic gradient
#'
#' Shared workhorse behind [`.reml_objective_rotated()`] and the optimiser
#' interface. The gradient uses the standard REML derivative
#' \eqn{\partial \ell / \partial \theta = -\tfrac12 [\operatorname{tr}(P
#' \, \partial V) - y^\top P \, \partial V \, P y]}, evaluated block-wise in
#' the rotated space: with \eqn{D_j = M_j (X^\top V^{-1} X)^{-1} M_j^\top},
#' \eqn{B_j = V_j^{-1} - V_j^{-1} D_j V_j^{-1}} and \eqn{s_j = V_j^{-1}
#' \widetilde r_j}, every parameter of the log-Cholesky factor
#' \eqn{\partial \Sigma = E L^\top + L E^\top} gives
#' \eqn{\operatorname{tr}(B_j \partial \Sigma) = 2 (L^\top B_j)_{b,a}} and
#' \eqn{s^\top \partial \Sigma s = 2 s_a (L^\top s)_b}.
#'
#' @inheritParams .reml_objective_rotated
#' @param need_grad Logical; skip gradient assembly when `FALSE`.
#' @return List with `value` (scalar log-likelihood) and `gradient`
#'   (numeric vector, `NULL` when `need_grad = FALSE`; both are
#'   `-1e10` / zeros for numerically illegal parameters).
#' @keywords internal
.reml_value_grad <- function(theta, m, lambda, Y_tilde, W_tilde,
                             chol_floor = 1e-8, reml = TRUE,
                             need_grad = TRUE) {
  ncp <- m * (m + 1L) / 2L
  bad <- list(value = -1e10,
              gradient = if (need_grad) numeric(2L * ncp) else NULL)
  if (length(theta) != 2L * ncp || any(!is.finite(theta))) return(bad)
  fac_G <- .logchol_unpack_factor(theta[seq_len(ncp)], m, chol_floor)
  fac_E <- .logchol_unpack_factor(theta[ncp + seq_len(ncp)], m, chol_floor)
  Sigma_G <- fac_G$Sigma
  Sigma_E <- fac_E$Sigma
  # extreme trial points can overflow the exp() inside the log-Cholesky
  # unpack, yielding non-finite covariance entries; treat as illegal
  if (any(!is.finite(Sigma_G)) || any(!is.finite(Sigma_E))) return(bad)

  n <- nrow(Y_tilde)
  q <- ncol(W_tilde)
  logdet_V <- 0
  y_Vinv_y <- 0
  XtVinvX <- matrix(0, nrow = q * m, ncol = q * m)
  XtVinvy <- numeric(q * m)
  A_arr <- array(0, dim = c(m, m, n))

  for (j in seq_len(n)) {
    V_j <- lambda[j] * Sigma_G + Sigma_E
    ch <- tryCatch(chol(V_j), error = function(e) NULL)
    if (is.null(ch)) return(bad)
    logdet_V <- logdet_V + 2 * sum(log(diag(ch)))
    A_j <- chol2inv(ch)
    A_arr[, , j] <- A_j
    y_j <- Y_tilde[j, ]
    w_j <- W_tilde[j, ]
    Ay <- A_j %*% y_j
    y_Vinv_y <- y_Vinv_y + sum(y_j * Ay)
    XtVinvX <- XtVinvX + kronecker(outer(w_j, w_j), A_j)
    XtVinvy <- XtVinvy + as.vector(outer(Ay, w_j))
  }

  if (any(!is.finite(XtVinvX)) || any(!is.finite(XtVinvy))) return(bad)
  inv <- .safe_inverse(XtVinvX)
  if (inv$status == "failed" || inv$rank < q * m) return(bad)
  G_x <- inv$inverse
  beta_hat <- drop(G_x %*% XtVinvy)
  quad <- y_Vinv_y - sum(XtVinvy * beta_hat)
  logdet_XtVinvX <- sum(log(inv$eigenvalues[inv$eigenvalues > 0]))

  value <- if (reml) {
    -0.5 * ((n * m - q * m) * log(2 * pi) + logdet_V + logdet_XtVinvX + quad)
  } else {
    -0.5 * (n * m * log(2 * pi) + logdet_V + quad)
  }
  if (!need_grad) return(list(value = value, gradient = NULL))

  # ---- analytic gradient ----------------------------------------------------
  # Residual r_j = y_j - B_hat' w_j and per-coordinate score s_j = V_j^{-1} r_j.
  B_hat <- matrix(beta_hat, nrow = m, ncol = q)  # B_hat[i, a]
  resid <- Y_tilde - W_tilde %*% t(B_hat)        # n x m, rows r_j
  # S_mat[j, i] = sum_k A_arr[i, k, j] * resid[j, k], i.e. s_j = A_j r_j.
  # (matrix(..., nrow = m) guards the m = 1 case against dimension dropping.)
  S_mat <- matrix(0, nrow = n, ncol = m)
  for (k in seq_len(m)) {
    S_mat <- S_mat + t(matrix(A_arr[, k, ], nrow = m)) * resid[, k]
  }

  # D_j = M_j G_x M_j' with M_j = I_m %x% w_j':
  # D_arr[i, i', j] = sum_{a,b} w_j[a] w_j[b] G_x[(a, i), (b, i')].
  G4 <- array(G_x, dim = c(m, q, m, q))          # [i, a, i', b]
  Gm <- matrix(aperm(G4, c(1L, 3L, 2L, 4L)), m * m, q * q)  # [(i,i'), (a,b)]
  S_w <- matrix(0, nrow = n, ncol = q * q)       # rows: vec(w_j w_j')
  for (b in seq_len(q)) {
    for (a in seq_len(q)) {
      S_w[, a + (b - 1L) * q] <- W_tilde[, a] * W_tilde[, b]
    }
  }
  D_arr <- array(Gm %*% t(S_w), dim = c(m, m, n))

  # B_j = A_j - A_j D_j A_j (per coordinate), vectorised over j.
  DA_arr <- array(0, dim = c(m, m, n))
  for (l in seq_len(m)) {
    for (i in seq_len(m)) {
      DA_arr[i, l, ] <- colSums(matrix(D_arr[i, , ], nrow = m) *
                                  matrix(A_arr[, l, ], nrow = m))
    }
  }
  B_arr <- A_arr
  for (l in seq_len(m)) {
    for (i in seq_len(m)) {
      B_arr[i, l, ] <- A_arr[i, l, ] -
        colSums(matrix(A_arr[i, , ], nrow = m) *
                  matrix(DA_arr[, l, ], nrow = m))
    }
  }

  Bmat <- matrix(B_arr, nrow = m * m, ncol = n)
  Bsum_G <- matrix(Bmat %*% lambda, nrow = m, ncol = m)
  Bsum_E <- matrix(rowSums(Bmat), nrow = m, ncol = m)

  grad_one <- function(fac, Bsum, weight) {
    T_mat <- S_mat %*% fac$L                    # rows: (L' s_j)'
    Q <- crossprod(weight * S_mat, T_mat)       # Q[a, b] = sum_j w s_a t_b
    G1 <- t(fac$L) %*% Bsum                     # G1[b, a]
    dL <- Q - t(G1)                             # dL[a, b]
    c(diag(dL) * fac$diag_scale, dL[lower.tri(dL)])
  }
  gradient <- c(
    grad_one(fac_G, Bsum_G, lambda),
    grad_one(fac_E, Bsum_E, rep(1, n))
  )
  if (any(!is.finite(gradient))) return(bad)
  list(value = value, gradient = gradient)
}

#' Build a cached fn/gr pair for the optimiser
#'
#' `optim`/`nlminb` evaluate the objective and the gradient at the same
#' trial point; caching halves the work.
#'
#' @inheritParams .reml_objective_rotated
#' @return List with functions `fn(theta)` (log-likelihood, to maximise)
#'   and `gr(theta)` (its gradient).
#' @keywords internal
.make_reml_objective <- function(m, lambda, Y_tilde, W_tilde, chol_floor) {
  cache <- new.env(parent = emptyenv())
  cache$theta <- NULL
  cache$vg <- NULL
  vg <- function(theta) {
    if (!is.null(cache$theta) && identical(theta, cache$theta)) {
      return(cache$vg)
    }
    value_grad <- .reml_value_grad(
      theta, m = m, lambda = lambda, Y_tilde = Y_tilde, W_tilde = W_tilde,
      chol_floor = chol_floor, reml = TRUE, need_grad = TRUE
    )
    cache$theta <- theta
    cache$vg <- value_grad
    value_grad
  }
  list(
    fn = function(theta) vg(theta)$value,
    gr = function(theta) vg(theta)$gradient
  )
}

#' GLS components of the rotated null model at fixed covariance matrices
#'
#' Recomputes, at the fitted (Sigma_G, Sigma_E), the per-coordinate inverse
#' covariances, their log-determinants and the normal-equation components
#' \eqn{X_0^\top V^{-1} X_0}, \eqn{X_0^\top V^{-1} y} and the GLS fixed
#' effects, without constructing the full \eqn{nm \times nm} matrices. All
#' inversions go through [`.safe_inverse()`] with explicit tolerances; the
#' log-determinants come from the same eigen decompositions.
#'
#' @param Sigma_G, Sigma_E Fitted `m x m` covariance matrices.
#' @param lambda Eigenvalues of `K` (length n).
#' @param Y_tilde, W_tilde Rotated data matrices.
#' @param inverse_tol Tolerance forwarded to `.safe_inverse()`.
#' @param keep_vinv Logical; keep the `m x m x n` array of inverses.
#' @return List with `Vinv`, `logdet_Vj`, `XtVinvX`, `XtVinvX_inv` and
#'   `beta_fixed`.
#' @keywords internal
.gls_components_rotated <- function(Sigma_G, Sigma_E, lambda, Y_tilde,
                                    W_tilde, inverse_tol, keep_vinv = TRUE) {
  n <- nrow(Y_tilde)
  m <- ncol(Y_tilde)
  q <- ncol(W_tilde)
  logdet_Vj <- numeric(n)
  Vinv <- if (keep_vinv) array(0, dim = c(m, m, n)) else NULL
  XtVinvX <- matrix(0, nrow = q * m, ncol = q * m)
  XtVinvy <- numeric(q * m)
  for (j in seq_len(n)) {
    V_j <- lambda[j] * Sigma_G + Sigma_E
    inv_j <- .safe_inverse(V_j, tol = inverse_tol)
    if (inv_j$status == "failed") {
      stop(sprintf("V_j inversion failed at rotated coordinate %d.", j),
           call. = FALSE)
    }
    logdet_Vj[j] <- sum(log(inv_j$eigenvalues[inv_j$eigenvalues > 0]))
    A_j <- inv_j$inverse
    if (keep_vinv) Vinv[, , j] <- A_j
    w_j <- W_tilde[j, ]
    XtVinvX <- XtVinvX + kronecker(outer(w_j, w_j), A_j)
    XtVinvy <- XtVinvy + as.vector(outer(A_j %*% Y_tilde[j, ], w_j))
  }
  inv <- .safe_inverse(XtVinvX, tol = inverse_tol)
  list(
    Vinv = Vinv,
    logdet_Vj = logdet_Vj,
    XtVinvX = XtVinvX,
    XtVinvX_inv = inv$inverse,
    beta_fixed = drop(inv$inverse %*% XtVinvy)
  )
}

#' Build the list of log-Cholesky starting points
#' @param Y Phenotype matrix.
#' @param m Number of traits.
#' @param n_starts Number of starts (first is deterministic).
#' @param h2_start Heritability split for the first start.
#' @param start_values Optional explicit first start (theta vector or
#'   `list(Sigma_G, Sigma_E)`).
#' @param chol_floor Diagonal floor used when projecting to PD.
#' @return List of theta vectors.
#' @keywords internal
.reml_starts <- function(Y, m, n_starts, h2_start, start_values, chol_floor) {
  S_P <- stats::cov(Y)
  S_P <- (S_P + t(S_P)) / 2
  project_pd <- function(S) {
    eig <- eigen(S, symmetric = TRUE)
    sweep(eig$vectors, 2L, pmax(eig$values, chol_floor^2), `*`) %*%
      t(eig$vectors)
  }
  theta_from_split <- function(h2) {
    SG <- project_pd(h2 * S_P)
    SE <- project_pd((1 - h2) * S_P)
    c(.logchol_pack(SG), .logchol_pack(SE))
  }
  starts <- list()
  if (!is.null(start_values)) {
    if (is.list(start_values)) {
      starts[[1L]] <- c(.logchol_pack(project_pd(start_values$Sigma_G)),
                        .logchol_pack(project_pd(start_values$Sigma_E)))
    } else {
      if (!is.numeric(start_values) ||
          length(start_values) != m * (m + 1L)) {
        stop("start_values must be a theta vector of length m * (m + 1) ",
             "or a list with Sigma_G and Sigma_E.", call. = FALSE)
      }
      starts[[1L]] <- as.numeric(start_values)
    }
    starts[[2L]] <- theta_from_split(h2_start)
  } else {
    starts[[1L]] <- theta_from_split(h2_start)
  }
  while (length(starts) < n_starts) {
    starts[[length(starts) + 1L]] <-
      theta_from_split(stats::runif(1L, min = 0.1, max = 0.9))
  }
  starts[seq_len(n_starts)]
}
