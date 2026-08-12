# Final joint signal-specific effect estimation (Stage 6B-4).
#
# THE single joint effect-estimation path. Both
# `estimate_mt_effects(conditioning_sets = ...)` and
# `resolve_locus_signals()` call this function. For a locus's final
# representatives S = {s1..sK}, with X_S the n x K genotype matrix and
# Z_S = I_m x X_S, and the FIXED global-null V and original residual
# maker P:
#
#   beta_S = (Z_S' P Z_S)^+ Z_S' P y,   Var(beta_S) = (Z_S' P Z_S)^+.
#
# All quantities are assembled with the Stage-3 rotated blocks; the
# null covariance is never re-estimated.

#' Joint GLS effects for a set of signal markers
#'
#' Computes \eqn{U_S = Z_S^\top P y}, \eqn{J_S = Z_S^\top P Z_S} for
#' the joint marker set with the fixed null projection (empty
#' conditioning set, i.e. the original P), then
#' \eqn{\widehat\beta_S = J_S^+ U_S} and \eqn{J_S^+}.
#'
#' @param null_fit A `"condped_mt_null"` fit with `rotation`.
#' @param G `n x p` genotype matrix.
#' @param signals Character vector (or column indices) of signal
#'   markers, in the desired block order.
#' @param marker_ids Marker order of `G`.
#' @param rank_tol Tolerance forwarded to [`.safe_inverse()`].
#'
#' @return A list with `signals`, `U`, `J`, `J_inv`, `beta` (K x m
#'   matrix, dimnamed by signal and trait), `se` (K x m),
#'   `covariance` (m x m x K array of per-signal blocks),
#'   `joint_covariance` (the full Km x Km matrix), `rank`,
#'   `condition`, `used_pseudoinverse` and `status`.
#' @keywords internal
.fit_joint_signal_effects <- function(null_fit, G, signals,
                                      marker_ids = colnames(G),
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
  trait_names <- null_fit$trait_names
  if (is.null(trait_names)) trait_names <- paste0("Trait", seq_len(m))

  if (is.null(marker_ids)) marker_ids <- paste0("M", seq_len(ncol(G)))
  if (length(signals) == 0L) {
    .stop_invalid_input("signals must contain at least one marker.")
  }
  if (is.character(signals)) {
    idx <- match(signals, marker_ids)
    if (anyNA(idx)) {
      .stop_invalid_input("signal markers not found in G: %s.",
                          signals[which(is.na(idx))[1L]])
    }
  } else {
    idx <- as.integer(signals)
    if (any(idx < 1L | idx > ncol(G))) {
      .stop_invalid_input("signals indices out of range.")
    }
  }
  signals <- marker_ids[idx]
  K <- length(idx)

  AM_arr <- .precompute_AM(rot)
  Ar <- .precompute_Ar(rot)
  A_arr <- rot$Vinv
  G_inv <- rot$XtVinvX_inv
  qm <- ncol(AM_arr)

  X_tilde <- matrix(0, nrow = n, ncol = K)
  for (i in seq_len(K)) {
    x <- G[, idx[i]]
    if (anyNA(x)) x[is.na(x)] <- mean(x, na.rm = TRUE)
    X_tilde[, i] <- crossprod(rot$U, x)
  }

  # U blocks: U_i = sum_j x_tilde[j, i] * (A_j r_j)
  U <- drop(Ar %*% X_tilde)          # m x K (column i = U_i)
  U_vec <- as.vector(U)              # (i-1)*m + t ordering

  # C_i = sum_j x_tilde[j, i] * (A_j M_j), m x qm per marker
  C_arr <- array(0, dim = c(m, qm, K))
  for (j in seq_len(n)) {
    AM <- AM_arr[, , j]
    for (i in seq_len(K)) {
      C_arr[, , i] <- C_arr[, , i] + X_tilde[j, i] * AM
    }
  }
  C_mat <- matrix(NA_real_, nrow = m * K, ncol = qm)
  for (i in seq_len(K)) {
    C_mat[(i - 1L) * m + seq_len(m), ] <- C_arr[, , i]
  }

  # J1 block (i, j) = sum_c x_i,c * x_j,c * A_c
  J1 <- matrix(0, m * K, m * K)
  for (i in seq_len(K)) {
    for (jj in i:K) {
      w <- X_tilde[, i] * X_tilde[, jj]
      blk <- matrix(0, m, m)
      for (c in seq_len(n)) {
        blk <- blk + w[c] * A_arr[, , c]
      }
      ri <- (i - 1L) * m + seq_len(m)
      rj <- (jj - 1L) * m + seq_len(m)
      J1[ri, rj] <- blk
      if (jj != i) J1[rj, ri] <- t(blk)
    }
  }
  J <- J1 - C_mat %*% G_inv %*% t(C_mat)
  J <- (J + t(J)) / 2

  inv <- .safe_inverse(J, tol = rank_tol, symmetric = TRUE)
  if (inv$status == "failed") {
    return(list(
      signals = signals, U = U_vec, J = J,
      J_inv = matrix(NA_real_, m * K, m * K),
      beta = matrix(NA_real_, K, m, dimnames = list(signals, trait_names)),
      se = matrix(NA_real_, K, m, dimnames = list(signals, trait_names)),
      covariance = array(NA_real_, dim = c(m, m, K)),
      joint_covariance = matrix(NA_real_, m * K, m * K),
      rank = 0L, condition = Inf, used_pseudoinverse = NA,
      status = .new_status(ok = FALSE, code = "unstable",
                           message = "Failed to invert the joint information matrix.")
    ))
  }
  J_inv <- inv$inverse
  beta_vec <- if (inv$rank > 0L) drop(J_inv %*% U_vec) else {
    rep(NA_real_, m * K)
  }
  beta <- matrix(NA_real_, K, m, dimnames = list(signals, trait_names))
  se <- beta
  covariance <- array(NA_real_, dim = c(m, m, K),
                      dimnames = list(trait_names, trait_names, signals))
  for (i in seq_len(K)) {
    rows <- (i - 1L) * m + seq_len(m)
    beta[i, ] <- beta_vec[rows]
    blk <- J_inv[rows, rows, drop = FALSE]
    covariance[, , i] <- blk
    se[i, ] <- if (inv$rank > 0L && all(is.finite(diag(blk)))) {
      sqrt(pmax(diag(blk), 0))
    } else {
      rep(NA_real_, m)
    }
  }

  list(
    signals = signals,
    U = U_vec,
    J = J,
    J_inv = J_inv,
    beta = beta,
    se = se,
    covariance = covariance,
    joint_covariance = J_inv,
    rank = as.integer(inv$rank),
    condition = inv$condition_number,
    used_pseudoinverse = inv$used_pseudoinverse,
    status = .new_status(
      ok = TRUE,
      code = if (inv$rank < m * K) "rank_deficient" else "ok",
      message = if (inv$rank < m * K) {
        "Joint signal information matrix is rank deficient; pseudo-inverse used."
      } else {
        ""
      }
    )
  )
}
