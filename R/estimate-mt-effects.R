#' Multi-trait marker effect estimation by GLS
#'
#' Estimates the multi-trait effects of selected loci by generalised least
#' squares, \eqn{\widehat\beta = J^+ U}, reusing the rotation object of the
#' null fit. The returned covariance equals \eqn{J^+}.
#'
#' @param null_fit An object of class `"condped_mt_null"` returned by
#'   [fit_mt_null()].
#' @param G Numeric `n x p` genotype dosage matrix coded 0/1/2.
#' @param targets Vector of target marker indices (integer positions in
#'   `G`) or identifiers (matching `marker_ids`) to estimate.
#' @param marker_ids Character vector of marker identifiers; defaults to
#'   `colnames(G)`.
#' @param conditioning_sets Optional list with one entry per target:
#'   a (possibly empty) character vector of conditioning markers. When
#'   supplied, each target's estimate is its coefficient/covariance
#'   block in the joint GLS of (conditioning set + target), computed
#'   through the single shared path
#'   [`.fit_joint_signal_effects()`]; `NULL` gives marginal effects.
#' @param rank_tol Relative tolerance used for the generalised inverse of
#'   \eqn{J}.
#' @param return_covariance Logical; return the per-locus `m x m` effect
#'   covariance matrices.
#'
#' @return A list with components `effects_long` (a data.frame with columns
#'   `marker_id`, `trait`, `beta`, `se`, `z`, `p_value`), `beta` (an
#'   n_loci x m matrix with marker and trait dimnames), `covariance` (an
#'   m x m x n_loci array), `se` (an n_loci x m matrix),
#'   `genotype_variance` (named numeric vector of per-locus variances of
#'   the non-missing centred genotypes), `locus_table` (a data.frame with
#'   `marker_id`, `maf`, `genotype_variance`, `rank_J`, `condition_J`,
#'   `status`), `trait_names`, `status` and `diagnostics`, as specified
#'   in the interface contract.
#' @export
estimate_mt_effects <- function(
    null_fit,
    G,
    targets,
    marker_ids = colnames(G),
    conditioning_sets = NULL,
    rank_tol = sqrt(.Machine$double.eps),
    return_covariance = TRUE
) {
  t0 <- proc.time()[["elapsed"]]
  .validate_null_fit(null_fit)
  if (is.null(null_fit$rotation)) {
    .stop_invalid_input(
      "null_fit does not contain the rotation object; re-fit with return_rotation = TRUE."
    )
  }
  rot <- null_fit$rotation
  m <- ncol(rot$Y_tilde)

  G <- as.matrix(G)
  if (!is.numeric(G) || !is.matrix(G)) {
    .stop_invalid_input("G must be a numeric matrix.")
  }
  if (nrow(G) != nrow(rot$Y_tilde)) {
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
  if (length(targets) == 0L) {
    .stop_invalid_input("targets must contain at least one marker.")
  }

  # Resolve targets to integer column indices
  if (is.character(targets)) {
    loci_idx <- match(targets, marker_ids)
    if (anyNA(loci_idx)) {
      .stop_invalid_input("Some target identifiers were not found in marker_ids.")
    }
  } else {
    loci_idx <- as.integer(targets)
    if (any(loci_idx < 1L | loci_idx > p)) {
      .stop_invalid_input("targets indices are out of range.")
    }
  }
  loci_idx <- unique(loci_idx)
  n_loci <- length(loci_idx)

  trait_names <- colnames(null_fit$Sigma_P_ref)
  A_arr <- rot$Vinv
  AM_arr <- .precompute_AM(rot)
  Ar <- .precompute_Ar(rot)
  G_inv <- rot$XtVinvX_inv

  # conditioning_sets: optional list, one entry per target; each entry is
  # a (possibly empty) character vector of conditioning markers. The
  # conditioned estimate of a target is its coefficient/covariance block
  # in the joint GLS of (conditioning set + target), computed through the
  # single shared path .fit_joint_signal_effects().
  conditional <- !is.null(conditioning_sets)
  if (conditional) {
    if (!is.list(conditioning_sets) ||
        length(conditioning_sets) != length(loci_idx)) {
      .stop_invalid_input(
        "conditioning_sets must be a list with one entry per target."
      )
    }
    conditioning_sets <- lapply(seq_along(conditioning_sets), function(i) {
      cs <- conditioning_sets[[i]]
      if (is.null(cs)) return(character())
      if (!is.character(cs) || anyNA(cs)) {
        .stop_invalid_input(
          "conditioning_sets entries must be character vectors of marker ids."
        )
      }
      cs <- unique(cs)
      tgt <- marker_ids[loci_idx[i]]
      if (tgt %in% cs) {
        .stop_invalid_input(
          "a conditioning set must not contain its own target (\"%s\").", tgt
        )
      }
      missing_cs <- setdiff(cs, marker_ids)
      if (length(missing_cs) > 0L) {
        .stop_invalid_input(
          "conditioning markers not found in marker_ids: %s.", missing_cs[1L]
        )
      }
      cs
    })
  }

  effects_list <- vector("list", n_loci)
  n_rank_deficient <- 0L
  locus_rows <- vector("list", n_loci)
  gv <- numeric(n_loci)

  for (i in seq_len(n_loci)) {
    idx <- loci_idx[i]
    x <- G[, idx]
    # Marker summaries use the original (unimputed) dosages; the score
    # uses mean-dosage imputation.
    n_eff <- sum(!is.na(x))
    allele_freq <- if (n_eff > 0L) mean(x, na.rm = TRUE) / 2 else NA_real_
    maf_l <- if (is.finite(allele_freq)) min(allele_freq, 1 - allele_freq) else NA_real_
    gv[i] <- if (n_eff > 1L) stats::var(x, na.rm = TRUE) else NA_real_

    if (!conditional) {
      if (anyNA(x)) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
      }
      x_tilde <- crossprod(rot$U, x)
      block <- .gls_block_components(x_tilde, AM_arr, Ar, A_arr, G_inv, rank_tol)
      beta_i <- block$beta
      cov_i <- block$J_inv
      rank_i <- block$rank
      status_i <- if (block$rank < m) "rank_deficient" else "ok"
      if (block$rank < m) n_rank_deficient <- n_rank_deficient + 1L
    } else {
      tgt_id <- marker_ids[idx]
      joint <- .fit_joint_signal_effects(
        null_fit, G, c(conditioning_sets[[i]], tgt_id),
        marker_ids = marker_ids, rank_tol = rank_tol
      )
      pos <- length(conditioning_sets[[i]]) + 1L
      beta_i <- joint$beta[pos, ]
      cov_i <- joint$covariance[, , pos]
      rank_i <- if (joint$rank > 0L) min(joint$rank, m) else 0L
      status_i <- joint$status$code
      if (status_i != "ok") n_rank_deficient <- n_rank_deficient + 1L
    }

    locus_rows[[i]] <- data.frame(
      marker_id = marker_ids[idx],
      maf = maf_l,
      genotype_variance = gv[i],
      rank_J = as.integer(rank_i),
      condition_J = if (!conditional) block$condition else joint$condition,
      status = status_i,
      stringsAsFactors = FALSE
    )

    effects_list[[i]] <- .format_one_effect(
      marker_id = marker_ids[idx],
      beta = beta_i,
      J_inv = if (isTRUE(return_covariance)) cov_i else matrix(NA_real_, m, m),
      rank = rank_i,
      trait_names = trait_names
    )
  }

  combined <- .combine_effects(effects_list)
  locus_ids <- marker_ids[loci_idx]
  dimnames(combined$beta) <- list(locus_ids, trait_names)
  se_mat <- matrix(combined$effects_long$se, nrow = n_loci, byrow = TRUE,
                   dimnames = list(locus_ids, trait_names))
  locus_table <- do.call(rbind, locus_rows)
  rownames(locus_table) <- NULL

  out <- list(
    effects_long = combined$effects_long,
    beta = combined$beta,
    covariance = if (isTRUE(return_covariance)) combined$covariance else NULL,
    se = se_mat,
    genotype_variance = stats::setNames(gv, locus_ids),
    locus_table = locus_table,
    trait_names = trait_names,
    status = .new_status(ok = TRUE, code = "ok", message = "",
                         warnings = character()),
    diagnostics = list(
      n_loci = n_loci,
      n_rank_deficient = n_rank_deficient,
      n_conditioned = if (conditional) {
        sum(lengths(conditioning_sets) > 0L)
      } else {
        0L
      },
      elapsed = proc.time()[["elapsed"]] - t0
    )
  )
  structure(out, class = "condped_effects")
}
