#' Multi-trait omnibus score scan
#'
#' Runs the joint multi-trait score test for every marker, strictly following
#' equations (9)-(12) of the Methods: \eqn{U_l = Z_l^\top P y},
#' \eqn{J_l = Z_l^\top P Z_l}, \eqn{Q_l = U_l^\top J_l^+ U_l}. The full
#' \eqn{Z_l} and \eqn{P} are never constructed explicitly; block computations
#' use the rotation object of the null fit.
#'
#' @param null_fit An object of class `"condped_mt_null"` returned by
#'   [fit_mt_null()].
#' @param G Numeric `n x p` genotype dosage matrix coded 0/1/2.
#' @param marker_ids Character vector of marker identifiers; defaults to
#'   `colnames(G)`.
#' @param maf_min Minimum minor allele frequency; markers below the threshold
#'   are filtered.
#' @param chunk_size Integer; number of markers processed per block.
#' @param rank_tol Relative tolerance used for the rank and generalised
#'   inverse of \eqn{J_l}.
#' @param return_score Logical; also return the score vectors.
#' @param return_effects Logical; also return marker effect estimates.
#' @param bootstrap_p Logical; compute bootstrap p-values.
#' @param n_boot Integer; number of bootstrap replicates when
#'   `bootstrap_p = TRUE`.
#' @param seed Optional random seed.
#'
#' @return A list with components `omnibus` (a data.frame with columns
#'   `marker_id`, `maf`, `n_eff`, `Q`, `df`, `p_value`, `rank_J`,
#'   `condition_J`, `status`, plus `filter_reason` for filtered markers),
#'   optional `score` and `effects`, plus `status` and `diagnostics`, as
#'   specified in the interface contract.
#' @export
scan_mt_omnibus <- function(
    null_fit,
    G,
    marker_ids = colnames(G),
    maf_min = 0.05,
    chunk_size = 2000L,
    rank_tol = sqrt(.Machine$double.eps),
    return_score = FALSE,
    return_effects = FALSE,
    bootstrap_p = FALSE,
    n_boot = 0L,
    seed = NULL
) {
  t0 <- proc.time()[["elapsed"]]
  .validate_null_fit(null_fit)
  if (is.null(null_fit$rotation)) {
    .stop_invalid_input(
      "null_fit does not contain the rotation object; re-fit with return_rotation = TRUE."
    )
  }
  rot <- null_fit$rotation
  n <- nrow(rot$Y_tilde)
  m <- ncol(rot$Y_tilde)

  G <- as.matrix(G)
  if (!is.numeric(G) || !is.matrix(G)) {
    .stop_invalid_input("G must be a numeric matrix.")
  }
  if (nrow(G) != n) {
    .stop_invalid_input("G must have the same number of rows (individuals) as the null-fit phenotypes.")
  }
  if (any(!is.na(G) & !is.finite(G))) {
    .stop_invalid_input("G must not contain non-finite values other than NA.")
  }
  p <- ncol(G)

  if (is.null(marker_ids)) {
    marker_ids <- paste0("M", seq_len(p))
  }
  if (length(marker_ids) != p) {
    .stop_invalid_input("length(marker_ids) must equal ncol(G).")
  }
  if (!is.numeric(maf_min) || length(maf_min) != 1L ||
      !is.finite(maf_min) || maf_min < 0 || maf_min > 0.5) {
    .stop_invalid_input("maf_min must be a single number in [0, 0.5].")
  }
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || !is.finite(chunk_size) || chunk_size < 1L) {
    .stop_invalid_input("chunk_size must be a positive integer.")
  }
  if (length(n_boot) != 1L || !is.finite(n_boot) || n_boot < 0L) {
    .stop_invalid_input("n_boot must be a non-negative integer.")
  }

  warnings <- character()
  if (isTRUE(bootstrap_p) && n_boot > 0L) {
    warnings <- c(warnings,
                  "Bootstrap p-values are not implemented in S3; ignored.")
  }

  # Precompute the rotated score information components once.
  A_arr <- rot$Vinv                  # m x m x n
  AM_arr <- .precompute_AM(rot)      # m x qm x n, coordinate j = A_j M_j
  Ar <- .precompute_Ar(rot)          # m x n,   column j = A_j r_tilde_j

  # Result containers
  omnibus <- vector("list", p)
  if (isTRUE(return_score)) score_list <- vector("list", p)
  if (isTRUE(return_effects)) effects_list <- vector("list", p)

  n_filtered <- 0L
  n_rank_deficient <- 0L
  n_tested <- 0L

  G_inv <- rot$XtVinvX_inv  # qm x qm

  # Chunked marker processing
  chunk_starts <- seq(1L, p, by = chunk_size)
  for (chunk_start in chunk_starts) {
    chunk_end <- min(p, chunk_start + chunk_size - 1L)
    idx <- chunk_start:chunk_end
    G_chunk <- G[, idx, drop = FALSE]
    ids_chunk <- marker_ids[idx]
    # Mean-dosage imputation of missing genotypes for the score computation;
    # MAF and n_eff below are computed from the original (unimputed) dosages.
    G_imp <- G_chunk
    if (anyNA(G_imp)) {
      col_mean <- colMeans(G_imp, na.rm = TRUE)
      for (s in seq_len(ncol(G_imp))) {
        if (anyNA(G_imp[, s]) && is.finite(col_mean[s])) {
          G_imp[is.na(G_imp[, s]), s] <- col_mean[s]
        }
      }
    }
    X_tilde <- crossprod(rot$U, G_imp)  # n x k

    for (s in seq_along(idx)) {
      marker_id <- ids_chunk[s]
      x <- G_chunk[, s]
      x_tilde <- X_tilde[, s]

      # Filtering based on the original dosage vector
      n_eff <- sum(!is.na(x))
      filter_reason <- ""
      passed <- TRUE
      maf <- NA_real_
      genotype_variance <- if (n_eff > 1L) {
        stats::var(x, na.rm = TRUE)
      } else {
        NA_real_
      }
      if (n_eff == 0L) {
        passed <- FALSE
        filter_reason <- "all_missing"
      } else {
        allele_freq <- mean(x, na.rm = TRUE) / 2
        maf <- min(allele_freq, 1 - allele_freq)
        if (!is.finite(maf) || maf <= 0) {
          passed <- FALSE
          filter_reason <- "monomorphic"
        } else if (maf < maf_min) {
          passed <- FALSE
          filter_reason <- "low_maf"
        }
      }

      if (!passed) {
        n_filtered <- n_filtered + 1L
        omnibus[[idx[s]]] <- data.frame(
          marker_id = marker_id,
          maf = maf,
          genotype_variance = genotype_variance,
          n_eff = n_eff,
          Q = NA_real_,
          df = 0L,
          p_value = NA_real_,
          rank_J = 0L,
          condition_J = Inf,
          status = "filtered",
          filter_reason = filter_reason,
          stringsAsFactors = FALSE
        )
        if (isTRUE(return_score)) score_list[[idx[s]]] <- NA_real_
        if (isTRUE(return_effects)) {
          effects_list[[idx[s]]] <- .format_one_effect(
            marker_id, beta = rep(NA_real_, m), J_inv = diag(NA_real_, m),
            rank = 0L, trait_names = colnames(null_fit$Sigma_P_ref)
          )
        }
        next
      }

      n_tested <- n_tested + 1L
      block <- .gls_block_components(x_tilde, AM_arr, Ar, A_arr, G_inv, rank_tol)

      if (block$rank == 0L) {
        n_rank_deficient <- n_rank_deficient + 1L
        status <- "rank_deficient"
        Q <- NA_real_
        df <- 0L
        p_value <- NA_real_
      } else {
        Q_raw <- sum(block$U * (block$J_inv %*% block$U))
        # Numerical negative tiny values are truncated to zero.
        neg_tol <- sqrt(.Machine$double.eps) * max(1, abs(Q_raw))
        if (Q_raw < 0 && Q_raw > -neg_tol) {
          Q <- 0
        } else if (Q_raw < 0) {
          Q <- NA_real_
          status <- "numerical_error"
          p_value <- NA_real_
        } else {
          Q <- Q_raw
          status <- "ok"
        }
        if (status == "ok" || !is.na(Q)) {
          df <- block$rank
          p_value <- stats::pchisq(Q, df = df, lower.tail = FALSE)
        } else {
          df <- block$rank
        }
        if (block$rank < m) {
          n_rank_deficient <- n_rank_deficient + 1L
          if (status == "ok") status <- "rank_deficient"
        }
      }

      omnibus[[idx[s]]] <- data.frame(
        marker_id = marker_id,
        maf = maf,
        genotype_variance = genotype_variance,
        n_eff = n_eff,
        Q = Q,
        df = as.integer(df),
        p_value = p_value,
        rank_J = as.integer(block$rank),
        condition_J = block$condition,
        status = status,
        filter_reason = filter_reason,
        stringsAsFactors = FALSE
      )
      if (isTRUE(return_score)) score_list[[idx[s]]] <- block$U
      if (isTRUE(return_effects)) {
        effects_list[[idx[s]]] <- .format_one_effect(
          marker_id, beta = block$beta, J_inv = block$J_inv,
          rank = block$rank, trait_names = colnames(null_fit$Sigma_P_ref)
        )
      }
    }
  }

  omnibus_df <- do.call(rbind, omnibus)
  rownames(omnibus_df) <- NULL

  out <- list(
    omnibus = omnibus_df,
    score = if (isTRUE(return_score)) score_list else NULL,
    effects = if (isTRUE(return_effects)) {
      .combine_effects(effects_list)
    } else {
      NULL
    },
    status = .new_status(
      ok = TRUE,
      code = "ok",
      message = "",
      warnings = warnings
    ),
    diagnostics = list(
      n_tested = n_tested,
      n_filtered = n_filtered,
      n_rank_deficient = n_rank_deficient,
      elapsed = proc.time()[["elapsed"]] - t0
    )
  )
  structure(out, class = "condped_omnibus")
}

#' Precompute the rotated components A_j M_j for all coordinates
#'
#' In the eigen-rotated space, \eqn{V^{-1}} is block-diagonal with blocks
#' \eqn{A_j = V_j^{-1}} and \eqn{X_0 = I_m \otimes \widetilde W}. For each
#' coordinate \eqn{j}, \eqn{M_j = I_m \otimes \widetilde w_j^\top} is the
#' `m x qm` local design block. The product \eqn{A_j M_j} is needed for the
#' score information correction term.
#'
#' @param rot The rotation object returned by [fit_mt_null()].
#' @return An `m x qm x n` array whose `j`-th slice is `A_j %*% M_j`.
#' @keywords internal
.precompute_AM <- function(rot) {
  m <- ncol(rot$Y_tilde)
  q <- ncol(rot$W_tilde)
  n <- nrow(rot$Y_tilde)
  A_arr <- rot$Vinv  # m x m x n

  # A_j M_j: for covariate a, columns (a-1)*m + (1:m) equal A_j * w_j[a]
  AM_arr <- array(0, dim = c(m, m * q, n))
  for (a in seq_len(q)) {
    cols <- (a - 1L) * m + seq_len(m)
    for (j in seq_len(n)) {
      AM_arr[, cols, j] <- A_arr[, , j] * rot$W_tilde[j, a]
    }
  }
  AM_arr
}

#' Precompute A_j %*% r_tilde_j for all rotated coordinates
#'
#' @param rot The rotation object.
#' @return An `m x n` matrix whose column `j` is `A_j %*% r_tilde_j`.
#' @keywords internal
.precompute_Ar <- function(rot) {
  m <- ncol(rot$Y_tilde)
  n <- nrow(rot$Y_tilde)
  A_arr <- rot$Vinv
  Ar <- matrix(0, nrow = m, ncol = n)
  for (j in seq_len(n)) {
    Ar[, j] <- A_arr[, , j] %*% rot$residual_tilde[j, ]
  }
  Ar
}

#' Compute the GLS score/information block for one rotated marker
#'
#' Implements \eqn{U_l = \sum_j \widetilde x_{lj} A_j \widetilde r_j} and
#' \eqn{J_l = \sum_j \widetilde x_{lj}^2 A_j -
#'   (\sum_j \widetilde x_{lj} A_j M_j) G
#'   (\sum_j \widetilde x_{lj} M_j^\top A_j)}.
#' The generalised inverse of \eqn{J_l} uses explicit rank tolerance.
#'
#' @param x_tilde Length-`n` rotated genotype vector for one marker.
#' @param AM_arr `m x qm x n` array from [`.precompute_AM()`].
#' @param Ar `m x n` matrix from [`.precompute_Ar()`].
#' @param G_inv The `qm x qm` matrix \eqn{(X_0^\top V^{-1} X_0)^{-1}}.
#' @param rank_tol Tolerance forwarded to [`.safe_inverse()`].
#' @return A list with components `U`, `J`, `J_inv`, `rank`, `condition`,
#'   `status`, `beta`.
#' @keywords internal
.gls_block_components <- function(x_tilde, AM_arr, Ar, A_arr, G_inv, rank_tol) {
  m <- nrow(Ar)
  qm <- ncol(AM_arr)
  n <- length(x_tilde)

  U <- drop(Ar %*% x_tilde)

  # J1 = sum_j x_tilde_j^2 A_j
  w <- x_tilde^2
  J1 <- matrix(0, m, m)
  for (j in seq_len(n)) {
    J1 <- J1 + w[j] * A_arr[, , j]
  }

  # C = sum_j x_tilde_j (A_j M_j), size m x qm
  C <- matrix(0, m, qm)
  for (j in seq_len(n)) {
    C <- C + x_tilde[j] * AM_arr[, , j]
  }

  # J = J1 - C G C'
  J <- J1 - C %*% G_inv %*% t(C)

  inv <- .safe_inverse(J, tol = rank_tol, symmetric = TRUE)
  rank <- inv$rank
  condition <- inv$condition_number

  if (inv$status == "failed") {
    return(list(
      U = U, J = J, J_inv = matrix(NA_real_, m, m),
      rank = 0L, condition = Inf, status = "failed",
      beta = rep(NA_real_, m)
    ))
  }

  J_inv <- inv$inverse
  beta <- if (rank > 0L) drop(J_inv %*% U) else rep(NA_real_, m)

  list(
    U = U,
    J = J,
    J_inv = J_inv,
    rank = as.integer(rank),
    condition = condition,
    status = inv$status,
    beta = beta
  )
}

#' Validate a null fit object for S3
#'
#' @param null_fit Candidate object.
#' @keywords internal
.validate_null_fit <- function(null_fit) {
  if (!inherits(null_fit, "condped_mt_null")) {
    .stop_invalid_input("null_fit must inherit class 'condped_mt_null'.")
  }
  if (!isTRUE(null_fit$status$ok)) {
    .stop_invalid_input(
      "null_fit did not converge (status$ok = FALSE); cannot proceed."
    )
  }
  req <- c("Sigma_G", "Sigma_E", "Sigma_P_ref", "fixed_effects", "rotation",
           "diagnostics", "status")
  missing <- setdiff(req, names(null_fit))
  if (length(missing) > 0L) {
    .stop_invalid_input(
      "null_fit is missing required components: %s.",
      paste(missing, collapse = ", ")
    )
  }
}

#' Format effect estimates for one marker
#'
#' @param marker_id Marker identifier.
#' @param beta Length-`m` effect vector.
#' @param J_inv `m x m` effect covariance matrix.
#' @param rank Numerical rank of `J`.
#' @param trait_names Character vector of trait names.
#' @return A list with `effects_long`, `beta`, `covariance`.
#' @keywords internal
.format_one_effect <- function(marker_id, beta, J_inv, rank, trait_names) {
  m <- length(beta)
  if (is.null(trait_names)) trait_names <- paste0("Trait", seq_len(m))
  if (rank == 0L || any(!is.finite(diag(J_inv)))) {
    se <- rep(NA_real_, m)
    z <- rep(NA_real_, m)
    p_value <- rep(NA_real_, m)
  } else {
    se <- sqrt(pmax(diag(J_inv), 0))
    z <- beta / se
    p_value <- stats::pchisq(z^2, df = 1, lower.tail = FALSE)
  }
  list(
    effects_long = data.frame(
      marker_id = marker_id,
      trait = trait_names,
      beta = beta,
      se = se,
      z = z,
      p_value = p_value,
      stringsAsFactors = FALSE
    ),
    beta = beta,
    covariance = J_inv
  )
}

#' Combine per-marker effect lists into scan/estimate return format
#'
#' @param effects_list List of outputs from [`.format_one_effect()`].
#' @return List with `effects_long`, `beta`, `covariance`.
#' @keywords internal
.combine_effects <- function(effects_list) {
  effects_long <- do.call(rbind, lapply(effects_list, `[[`, "effects_long"))
  rownames(effects_long) <- NULL
  beta <- do.call(rbind, lapply(effects_list, `[[`, "beta"))
  if (!is.matrix(beta)) beta <- matrix(beta, nrow = 1L)
  # Stable ID contract: rows are marker ids, columns trait names.
  dimnames(beta) <- list(
    vapply(effects_list, function(e) e$effects_long$marker_id[1L],
           character(1)),
    as.character(effects_list[[1L]]$effects_long$trait)
  )
  covariance <- simplify2array(lapply(effects_list, `[[`, "covariance"))
  list(
    effects_long = effects_long,
    beta = beta,
    covariance = covariance
  )
}
