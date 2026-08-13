# Stage 6C: comparison pipelines (lead-ASSET / resolved-ASSET / CondPED).
#
# One orchestrator `.run_comparison_pipelines()` plus three thin
# per-pipeline workers:
#   * lead_asset:     per locus, the MARGINAL trait-wise effects of the
#                     locus lead SNP (estimate_mt_effects, no conditional
#                     resolution) are fed to the ASSET adapter
#                     (.run_asset_comparison).
#   * resolved_asset: per RESOLVED signal, the FINAL JOINT signal effects
#                     from resolve_locus_signals()$beta are fed to the
#                     adapter. Stepwise temporary conditional betas are
#                     never used.
#   * condped_full:   the full CondPED chain from the SAME resolution:
#                     attribute_traits() (within-signal Holm) ->
#                     derive_conditional_contrasts(fit$Sigma_P_ref) ->
#                     decompose_conditional_effects().
# The resolution is computed ONCE and shared by resolved_asset and
# condped_full. The ASSET backend hook is threaded through so tests can
# inject a fake ASSET.

#' Aggregate per-unit adapter statuses into one pipeline status
#'
#' @param codes Character vector of per-locus/per-signal statuses.
#' @return `"ok"` when every unit is ok (or there are none), otherwise
#'   the sorted unique non-ok codes joined by `";"`.
#' @keywords internal
.pipeline_status <- function(codes) {
  codes <- as.character(codes)
  bad <- codes[codes != "ok"]
  if (length(bad) == 0L) return("ok")
  paste(sort(unique(bad)), collapse = ";")
}

#' Lead-ASSET pipeline: marginal lead-SNP effects through ASSET
#'
#' For every locus, takes the locus lead SNP's MARGINAL trait-wise
#' beta/se from [estimate_mt_effects()] (no conditional resolution) and
#' runs the ASSET adapter once per locus.
#'
#' @param fit The global null fit (with rotation).
#' @param G `n x p` genotype matrix (must have marker-id colnames).
#' @param locus_object List with `loci` and `membership` data.frames.
#' @param omnibus Omnibus scan data.frame (from [scan_mt_omnibus()]),
#'   used only to annotate the lead's marginal omnibus p-value.
#' @param trait_names Trait names (order of `fit$Sigma_P_ref`).
#' @param asset_backend Backend hook for [`.run_asset_comparison()`].
#'
#' The ASSET correlation matrix is derived per locus from the SAME
#' marginal effect fit as beta/se: `cov2cor()` of the marginal
#' covariance block — never from `Sigma_P_ref`.
#'
#' @return List with `loci` (named list per `locus_id`; each entry has
#'   `locus_id`, `lead_snp`, `lead_p`, `omnibus_p`, named `beta`/`se`
#'   vectors, `Sigma_Z` and the frozen adapter fields), `n_loci` and
#'   `status`.
#' @keywords internal
.pipeline_lead_asset <- function(fit, G, locus_object, omnibus,
                                 trait_names, asset_backend) {
  loci <- locus_object$loci
  per_locus <- vector("list", nrow(loci))
  ids <- character(nrow(loci))
  for (i in seq_len(nrow(loci))) {
    lid <- loci$locus_id[i]
    lead <- loci$lead_snp[i]
    eff <- estimate_mt_effects(fit, G, targets = lead)
    el <- eff$effects_long
    beta <- stats::setNames(el$beta, el$trait)[trait_names]
    se <- stats::setNames(el$se, el$trait)[trait_names]
    Sigma_Z <- .sigma_z_from_cov(eff$covariance[, , 1L])
    omp <- omnibus$p_value[match(lead, omnibus$marker_id)]
    adapter <- if (is.null(Sigma_Z)) {
      c(.asset_empty_result("rank_deficient"),
        list(settings = list(scr_pthr = 0.05, two_sided = TRUE,
                             Neff = nrow(fit$rotation$Y_tilde))))
    } else {
      .run_asset_comparison(
        beta, se = se, Sigma_Z = Sigma_Z, trait_names = trait_names,
        sample_size = nrow(fit$rotation$Y_tilde),
        backend = asset_backend
      )
    }
    per_locus[[i]] <- c(
      list(
        locus_id = lid,
        lead_snp = lead,
        lead_p = loci$lead_p[i],
        omnibus_p = if (length(omp) == 0L) NA_real_ else unname(omp),
        beta = beta,
        se = se,
        Sigma_Z = Sigma_Z
      ),
      adapter
    )
    ids[i] <- lid
  }
  names(per_locus) <- ids
  list(
    loci = per_locus,
    n_loci = nrow(loci),
    status = .pipeline_status(vapply(per_locus, `[[`, character(1),
                                     "status"))
  )
}

#' Resolved-ASSET pipeline: final joint signal effects through ASSET
#'
#' For every RESOLVED signal of a [resolve_locus_signals()] result,
#' feeds the signal's FINAL JOINT trait-wise beta/se (from the
#' resolution's `beta` long table — never stepwise temporary betas) to
#' the ASSET adapter. The ASSET correlation matrix is `cov2cor()` of
#' that signal's final joint covariance block.
#'
#' @param resolution A [resolve_locus_signals()] result.
#' @param trait_names Trait names (order of `fit$Sigma_P_ref`).
#' @param asset_backend Backend hook for [`.run_asset_comparison()`].
#' @param sample_size Total analysed sample size per trait.
#'
#' @return List with `signals` (named list per `signal_id`; each entry
#'   has `locus_id`, `signal_id`, `representative_snp`, named
#'   `beta`/`se` vectors, `Sigma_Z` and the frozen adapter fields),
#'   `n_signals` and `status`.
#' @keywords internal
.pipeline_resolved_asset <- function(resolution, trait_names,
                                     asset_backend, sample_size) {
  b <- resolution$beta
  if (is.null(b) || nrow(b) == 0L) {
    return(list(signals = stats::setNames(list(), character()),
                n_signals = 0L, status = "ok"))
  }
  m <- length(trait_names)
  per_signal <- lapply(seq_len(nrow(resolution$signals)), function(k) {
    sig_row <- resolution$signals[k, ]
    sid <- sig_row$signal_id
    rows <- b[b$signal_id == sid, , drop = FALSE]
    beta <- stats::setNames(rows$beta, rows$trait)[trait_names]
    se <- stats::setNames(rows$se, rows$trait)[trait_names]
    # this signal's block of its locus's final joint covariance
    joint_cov <- resolution$covariance[[sig_row$locus_id]]
    ord <- sig_row$signal_order
    blk_idx <- (ord - 1L) * m + seq_len(m)
    cov_block <- if (!is.null(joint_cov) &&
                     nrow(joint_cov) >= max(blk_idx)) {
      joint_cov[blk_idx, blk_idx, drop = FALSE]
    } else {
      matrix(NA_real_, m, m)
    }
    Sigma_Z <- .sigma_z_from_cov(cov_block)
    adapter <- if (is.null(Sigma_Z)) {
      c(.asset_empty_result("rank_deficient"),
        list(settings = list(scr_pthr = 0.05, two_sided = TRUE,
                             Neff = sample_size)))
    } else {
      .run_asset_comparison(
        beta, se = se, Sigma_Z = Sigma_Z, trait_names = trait_names,
        sample_size = sample_size,
        backend = asset_backend
      )
    }
    c(
      list(
        locus_id = rows$locus_id[1L],
        signal_id = sid,
        representative_snp = rows$representative_snp[1L],
        beta = beta,
        se = se,
        Sigma_Z = Sigma_Z
      ),
      adapter
    )
  })
  sig_ids <- vapply(per_signal, `[[`, character(1), "signal_id")
  names(per_signal) <- sig_ids
  list(
    signals = per_signal,
    n_signals = length(sig_ids),
    status = .pipeline_status(vapply(per_signal, `[[`, character(1),
                                     "status"))
  )
}

#' ASSET correlation matrix from an effect covariance block
#'
#' `cov2cor()` of the block, or `NULL` when the block is unusable
#' (non-finite or non-positive diagonal).
#'
#' @param cov_block `m x m` covariance matrix.
#' @return Correlation matrix or `NULL`.
#' @keywords internal
.sigma_z_from_cov <- function(cov_block) {
  if (is.null(cov_block) || any(!is.finite(cov_block))) return(NULL)
  d <- diag(cov_block)
  if (any(d <= 0)) return(NULL)
  out <- stats::cov2cor(cov_block)
  if (any(!is.finite(out))) return(NULL)
  out
}

#' CondPED full pipeline: attribution + conditional subset decomposition
#'
#' Continues from the SAME signal resolution used by the resolved-ASSET
#' pipeline: [attribute_traits()] (within-signal Holm) defines candidate
#' trait sets, [derive_conditional_contrasts()] freezes the reference
#' covariance from `fit$Sigma_P_ref`, and
#' [decompose_conditional_effects()] evaluates the conditional
#' representation loss over candidate subsets.
#'
#' @param fit The global null fit (uses `fit$Sigma_P_ref`).
#' @param resolution A [resolve_locus_signals()] result.
#' @param alpha_trait Within-signal Holm level for the attribution.
#' @param tolerance Primary representation-loss tolerance tau.
#'
#' @return List with `candidate_sets`, `signal_table`, `trait_table`,
#'   `subset_analysis` (list with `subset_table`,
#'   `minimum_representative_sets`, `irreducible_modules`,
#'   `tolerance_path`) and `status` (aggregation of the attribution and
#'   decomposition status codes).
#' @keywords internal
.pipeline_condped_full <- function(fit, resolution, alpha_trait,
                                   tolerance) {
  attribution <- attribute_traits(
    resolution, candidate_mode = "holm_fwer", alpha_trait = alpha_trait
  )
  contrasts <- derive_conditional_contrasts(fit$Sigma_P_ref)
  decomp <- decompose_conditional_effects(
    resolution, attribution, contrasts, tolerance = tolerance
  )
  list(
    candidate_sets = attribution$candidate_sets,
    signal_table = attribution$signal_table,
    trait_table = attribution$trait_table,
    subset_analysis = list(
      subset_table = decomp$subset_table,
      minimum_representative_sets = decomp$minimum_representative_sets,
      irreducible_modules = decomp$irreducible_modules,
      tolerance_path = decomp$tolerance_path
    ),
    status = .pipeline_status(c(attribution$status$code,
                                decomp$status$code))
  )
}

#' Run the three Stage-6C comparison pipelines
#'
#' Orchestrates the lead-ASSET, resolved-ASSET and CondPED-full
#' comparison pipelines over a common locus set. The signal resolution
#' ([resolve_locus_signals()]) is computed once and shared by the
#' resolved-ASSET and CondPED-full pipelines. The trait z-score
#' correlation matrix handed to ASSET is `cov2cor(fit$Sigma_P_ref)`
#' (correlated traits supported; never diagonalised). ASSET itself is
#' optional: with `asset_backend = NULL` and ASSET not installed, the
#' ASSET pipelines report status `"asset_not_available"` without
#' erroring.
#'
#' @param fit The single global null fit (with rotation).
#' @param G `n x p` genotype dosage matrix with marker-id colnames.
#' @param locus_object A [define_associated_loci()] result (or a list
#'   with `loci` and `membership` data.frames; `loci` must carry
#'   `locus_id`, `lead_snp`, `lead_p`).
#' @param scan A [scan_mt_omnibus()] result or its omnibus data.frame;
#'   used to annotate lead omnibus p-values.
#' @param trait_names Trait names; must match
#'   `colnames(fit$Sigma_P_ref)` (canonical order is taken from the
#'   fit).
#' @param alpha_signal Signal-level alpha: forwarded to the resolution
#'   screening rule and used as the within-signal Holm level in the
#'   CondPED-full attribution.
#' @param signal_adjust Local screening rule for the resolution.
#' @param max_signals Maximum signals per locus in the resolution.
#' @param asset_backend Backend hook forwarded to
#'   [`.run_asset_comparison()`] (test injection point for a fake
#'   ASSET).
#' @param tolerance Primary representation-loss tolerance for the
#'   CondPED-full decomposition.
#' @param ... Forwarded to [resolve_locus_signals()] (e.g. `rank_tol`).
#'
#' @return A list with components
#' \describe{
#'   \item{lead_asset}{Result of [`.pipeline_lead_asset()`].}
#'   \item{resolved_asset}{Result of [`.pipeline_resolved_asset()`].}
#'   \item{condped_full}{Result of [`.pipeline_condped_full()`].}
#'   \item{status}{`"ok"` when all three pipeline statuses are `"ok"`,
#'     otherwise the sorted unique non-ok codes joined by `";"`.}
#' }
#' @keywords internal
.run_comparison_pipelines <- function(fit, G, locus_object, scan,
                                      trait_names, alpha_signal = 0.05,
                                      signal_adjust = "within_locus_bonferroni",
                                      max_signals = 10L,
                                      asset_backend = NULL,
                                      tolerance = 0.10, ...) {
  .validate_null_fit(fit)
  if (is.null(fit$rotation)) {
    .stop_invalid_input(
      "fit does not contain the rotation object; re-fit with return_rotation = TRUE."
    )
  }
  if (is.null(fit$Sigma_P_ref)) {
    .stop_invalid_input("fit does not contain Sigma_P_ref.")
  }
  fit_traits <- colnames(fit$Sigma_P_ref)
  if (!is.character(trait_names) || anyNA(trait_names) ||
      anyDuplicated(trait_names) ||
      !setequal(trait_names, fit_traits)) {
    .stop_invalid_input(
      "trait_names must match the trait names of fit$Sigma_P_ref."
    )
  }
  trait_names <- fit_traits   # canonical order from the fit

  G <- as.matrix(G)
  if (!is.numeric(G) || is.null(colnames(G))) {
    .stop_invalid_input("G must be a numeric matrix with marker-id colnames.")
  }
  if (!is.list(locus_object) || is.null(locus_object$loci) ||
      is.null(locus_object$membership)) {
    .stop_invalid_input(
      "locus_object must provide loci and membership data.frames."
    )
  }
  if (!all(c("locus_id", "lead_snp", "lead_p") %in%
           names(locus_object$loci))) {
    .stop_invalid_input(
      "locus_object$loci must contain locus_id, lead_snp, lead_p."
    )
  }
  omnibus <- .as_omnibus_table(scan)
  if (!is.null(asset_backend) && !is.function(asset_backend)) {
    .stop_invalid_input("asset_backend must be NULL or a function.")
  }
  .check_prob(alpha_signal, "alpha_signal")
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0 || tolerance >= 1) {
    .stop_invalid_input("tolerance must be a single number in [0, 1).")
  }

  # Pipeline 1: marginal lead-SNP effects through ASSET (no resolution).
  # Each pipeline derives its ASSET correlation matrix from the SAME
  # effect fit that produced the beta/se it submits.
  lead_asset <- .pipeline_lead_asset(
    fit, G, locus_object, omnibus, trait_names, asset_backend
  )

  # Shared signal resolution (computed once).
  resolution <- resolve_locus_signals(
    fit, G, locus_object,
    signal_adjust = signal_adjust,
    alpha_signal = alpha_signal,
    max_signals = max_signals,
    ...
  )

  # Pipeline 2: final joint signal effects through ASSET.
  resolved_asset <- .pipeline_resolved_asset(
    resolution, trait_names, asset_backend,
    sample_size = nrow(fit$rotation$Y_tilde)
  )

  # Pipeline 3: full CondPED chain from the same resolution.
  condped_full <- .pipeline_condped_full(
    fit, resolution, alpha_trait = alpha_signal, tolerance = tolerance
  )

  list(
    lead_asset = lead_asset,
    resolved_asset = resolved_asset,
    condped_full = condped_full,
    status = .pipeline_status(c(lead_asset$status,
                                resolved_asset$status,
                                condped_full$status))
  )
}
