# Conditional projection and conditional multi-trait scan (Stage 6B-2).
#
# Implements Methods section 3.1: conditioning on a set of SNPs by
# augmenting the fixed-effect design, X*_C = [X0, Z_C] with
# Z_C = I_m x X_C, while the global-null covariance V stays fixed.
# P_C = V^-1 - V^-1 X*_C (X*_C' V^-1 X*_C)^+ X*_C' V^-1 is rebuilt for
# every conditioning set; Sigma_G, Sigma_E and V never change. No
# locus-wise REML anywhere. All linear algebra reuses the Stage-2/3
# rotated-block machinery (`.safe_inverse()` tolerances and
# `.gls_block_components()`); the score formulas are not
# re-implemented.

#' Build the conditional projection for a conditioning genotype set
#'
#' Augments the null-model fixed-effect design with
#' \eqn{Z_C = I_m \otimes X_C} and constructs the residual-projection
#' components under the FIXED global-null \eqn{\widehat V}. An empty
#' `X_C` reduces exactly to the original scan projection.
#'
#' In the eigen-rotated space every per-coordinate block uses the same
#' \eqn{A_j = V_j^{-1}} blocks from `null_fit$rotation$Vinv`; the
#' augmented design values per coordinate are
#' \eqn{d_j = (\widetilde w_j, \widetilde x_{C,j})}, giving
#' \eqn{X^{*\top}V^{-1}X^* = \sum_j d_j d_j^\top \otimes A_j} (in the
#' trait-fastest \eqn{(a-1)m + t} ordering of the stacked design).
#'
#' @param null_fit A `"condped_mt_null"` fit with `rotation`.
#' @param X_C `NULL` or an `n x c` conditioning genotype matrix (or a
#'   length-`n` vector for one SNP). Missing dosages are mean-imputed
#'   as in [scan_mt_omnibus()].
#' @param rank_tol Tolerance forwarded to [`.safe_inverse()`].
#'
#' @return A list with components `AM_arr`, `Ar`, `A_arr` (the fixed
#'   \eqn{V^{-1}} blocks, identical to `null_fit$rotation$Vinv`),
#'   `G_inv` (pseudo-inverse of the augmented
#'   \eqn{X^{*\top}V^{-1}X^*}), `S_XX`, `beta_star`, `residual_tilde`,
#'   `U`, `lambda`, `n_conditioning`, `Sigma_G`, `Sigma_E`
#'   (references, unchanged), `rank_S_XX`, `condition_S_XX`,
#'   `used_pseudoinverse` and `status`.
#' @keywords internal
.build_conditional_projection <- function(null_fit, X_C = NULL,
                                          rank_tol = sqrt(.Machine$double.eps)) {
  .validate_null_fit(null_fit)
  if (is.null(null_fit$rotation)) {
    .stop_invalid_input(
      "null_fit does not contain the rotation object; re-fit with return_rotation = TRUE."
    )
  }
  rot <- null_fit$rotation
  m <- ncol(rot$Y_tilde)
  n <- nrow(rot$Y_tilde)
  q <- ncol(rot$W_tilde)
  A_arr <- rot$Vinv   # fixed V^{-1} blocks; never modified

  if (is.null(X_C)) X_C <- matrix(0, nrow = n, ncol = 0L)
  if (is.vector(X_C)) X_C <- matrix(X_C, ncol = 1L)
  if (!is.numeric(X_C) || !is.matrix(X_C) || nrow(X_C) != n) {
    .stop_invalid_input("X_C must be an n x c numeric matrix (or NULL).")
  }
  c_cond <- ncol(X_C)
  if (c_cond > 0L && anyNA(X_C)) {
    col_mean <- colMeans(X_C, na.rm = TRUE)
    for (s in seq_len(c_cond)) {
      if (anyNA(X_C[, s])) X_C[is.na(X_C[, s]), s] <- col_mean[s]
    }
  }
  XC_tilde <- if (c_cond > 0L) {
    crossprod(rot$U, X_C)
  } else {
    matrix(0, nrow = n, ncol = 0L)
  }

  # Combined per-coordinate design values and the A_j y_j products.
  D <- cbind(rot$W_tilde, XC_tilde)              # n x (q + c)
  k <- m * (q + c_cond)
  Ay <- matrix(0, nrow = m, ncol = n)
  for (j in seq_len(n)) {
    Ay[, j] <- A_arr[, , j] %*% rot$Y_tilde[j, ]
  }

  S_XX <- matrix(0, k, k)
  S_Xy <- numeric(k)
  AM_arr <- array(0, dim = c(m, k, n))
  for (j in seq_len(n)) {
    dj <- D[j, ]
    Aj <- A_arr[, , j]
    S_XX <- S_XX + kronecker(tcrossprod(dj), Aj)
    S_Xy <- S_Xy + kronecker(dj, Ay[, j])
    AM_arr[, , j] <- kronecker(matrix(dj, nrow = 1L), Aj)
  }
  S_XX <- (S_XX + t(S_XX)) / 2

  inv <- .safe_inverse(S_XX, tol = rank_tol, symmetric = TRUE)
  if (inv$status == "failed") {
    return(list(
      AM_arr = NULL, Ar = NULL, A_arr = A_arr, G_inv = NULL,
      S_XX = S_XX, beta_star = NULL, residual_tilde = NULL,
      U = rot$U, lambda = rot$lambda, n_conditioning = c_cond,
      Sigma_G = null_fit$Sigma_G, Sigma_E = null_fit$Sigma_E,
      rank_S_XX = NA_integer_, condition_S_XX = NA_real_,
      used_pseudoinverse = NA,
      status = .new_status(ok = FALSE, code = "unstable",
                           message = "Failed to invert the augmented fixed-effect information matrix.")
    ))
  }
  beta_star <- drop(inv$inverse %*% S_Xy)
  B_mat <- matrix(beta_star, nrow = m, ncol = q + c_cond)
  residual_tilde <- rot$Y_tilde - t(B_mat %*% t(D))
  Ar <- matrix(0, nrow = m, ncol = n)
  for (j in seq_len(n)) {
    Ar[, j] <- A_arr[, , j] %*% residual_tilde[j, ]
  }

  list(
    AM_arr = AM_arr,
    Ar = Ar,
    A_arr = A_arr,
    G_inv = inv$inverse,
    S_XX = S_XX,
    beta_star = beta_star,
    residual_tilde = residual_tilde,
    U = rot$U,
    lambda = rot$lambda,
    n_conditioning = c_cond,
    Sigma_G = null_fit$Sigma_G,
    Sigma_E = null_fit$Sigma_E,
    rank_S_XX = inv$rank,
    condition_S_XX = inv$condition_number,
    used_pseudoinverse = inv$used_pseudoinverse,
    status = .new_status(
      ok = TRUE,
      code = if (inv$rank < k) "rank_deficient" else "ok",
      message = if (inv$rank < k) {
        "The augmented fixed-effect information matrix is rank deficient (collinear conditioning input); pseudo-inverse used."
      } else {
        ""
      }
    )
  )
}

#' Conditional multi-trait scan under a fixed conditional projection
#'
#' For each candidate SNP, computes the conditional score
#' \eqn{U_{j|C} = Z_j^\top P_C y}, information
#' \eqn{J_{j|C} = Z_j^\top P_C Z_j}, statistic
#' \eqn{Q_{j|C} = U^\top J^+ U} with `df = rank(J)` and the chi-square
#' p-value, reusing [`.gls_block_components()`]. The returned
#' `beta`/`covariance` are conditional-scan temporary estimates only
#' and must not be interpreted as final signal-specific effects.
#'
#' @param proj A projection from [`.build_conditional_projection()`].
#' @param G `n x p` genotype matrix of candidate SNPs.
#' @param marker_ids Marker ids; defaults to `colnames(G)`.
#' @param rank_tol Tolerance forwarded to [`.safe_inverse()`].
#' @param return_effects Logical; also return per-SNP `beta`/`covariance`.
#'
#' @return A list with `conditional` (data.frame: `marker_id`, `Q`,
#'   `df`, `p_value`, `rank_J`, `condition_J`, `status`), `effects`
#'   (list with `beta` matrix and `covariance` array, or `NULL`),
#'   `n_conditioning`, `status` and `diagnostics`.
#' @keywords internal
.conditional_mt_scan <- function(proj, G, marker_ids = colnames(G),
                                 rank_tol = sqrt(.Machine$double.eps),
                                 return_effects = TRUE) {
  if (is.null(proj$AM_arr) || is.null(proj$G_inv)) {
    .stop_invalid_input(
      "proj does not contain usable projection components (build failed?)."
    )
  }
  n <- nrow(G)
  p <- ncol(G)
  if (is.null(marker_ids)) marker_ids <- paste0("marker", seq_len(p))
  if (length(marker_ids) != p || anyDuplicated(marker_ids)) {
    .stop_invalid_input("marker_ids must be unique and match ncol(G).")
  }
  m <- nrow(proj$Ar)
  trait_names <- rownames(proj$A_arr)
  if (is.null(trait_names)) trait_names <- paste0("Trait", seq_len(m))

  rows <- vector("list", p)
  effects_list <- vector("list", p)
  n_rank_deficient <- 0L
  for (i in seq_len(p)) {
    x <- G[, i]
    if (anyNA(x)) x[is.na(x)] <- mean(x, na.rm = TRUE)
    x_tilde <- crossprod(proj$U, x)
    block <- .gls_block_components(x_tilde, proj$AM_arr, proj$Ar,
                                   proj$A_arr, proj$G_inv, rank_tol)
    qs <- .q_from_block(block, m)
    if (block$rank < m) n_rank_deficient <- n_rank_deficient + 1L
    rows[[i]] <- data.frame(
      marker_id = marker_ids[i],
      Q = qs$Q,
      df = as.integer(block$rank),
      p_value = qs$p_value,
      rank_J = as.integer(block$rank),
      condition_J = block$condition,
      status = qs$status,
      stringsAsFactors = FALSE
    )
    if (isTRUE(return_effects)) {
      effects_list[[i]] <- .format_one_effect(
        marker_ids[i], block$beta, block$J_inv, block$rank, trait_names
      )
    }
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL

  list(
    conditional = tab,
    effects = if (isTRUE(return_effects)) {
      .combine_effects(effects_list)
    } else {
      NULL
    },
    n_conditioning = proj$n_conditioning,
    status = .new_status(ok = TRUE, code = "ok"),
    diagnostics = list(
      n_tested = p,
      n_rank_deficient = n_rank_deficient,
      proj_status = proj$status$code
    )
  )
}

#' Score statistic and status from one GLS block (shared convention)
#'
#' Same truncation/status convention as [scan_mt_omnibus()]:
#' numerically tiny negative Q values are truncated to zero, clearly
#' negative ones are reported as `numerical_error` with NA p-value.
#'
#' @param block Output of [`.gls_block_components()`].
#' @param m Number of traits.
#' @return List with `Q`, `p_value`, `status`.
#' @keywords internal
.q_from_block <- function(block, m) {
  if (block$rank == 0L) {
    return(list(Q = NA_real_, p_value = NA_real_,
                status = "rank_deficient"))
  }
  Q_raw <- sum(block$U * (block$J_inv %*% block$U))
  neg_tol <- sqrt(.Machine$double.eps) * max(1, abs(Q_raw))
  if (Q_raw < 0 && Q_raw > -neg_tol) {
    Q <- 0
  } else if (Q_raw < 0) {
    return(list(Q = NA_real_, p_value = NA_real_,
                status = "numerical_error"))
  } else {
    Q <- Q_raw
  }
  status <- if (block$rank < m) "rank_deficient" else "ok"
  list(
    Q = Q,
    p_value = stats::pchisq(Q, df = block$rank, lower.tail = FALSE),
    status = status
  )
}
