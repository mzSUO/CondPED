#' End-to-end CondPED analysis (signal-level, revised v1.0)
#'
#' The thinnest possible orchestrator: it only calls the existing
#' public functions in the frozen order and never re-implements any
#' formula. `fit_mt_null()` is called exactly once; there is no LOCO
#' branch; `chromosome`/`position` are marker annotations.
#'
#' @details
#' ## `signal_mode = "resolve"`
#'
#' \preformatted{
#' validate inputs
#' -> construct/read one fixed K
#' -> fit_mt_null() exactly once
#' -> scan_mt_omnibus()
#' -> genome-wide p adjustment -> selected_markers
#' -> define_associated_loci()
#' -> resolve_locus_signals() (stepwise + final joint refit)
#' -> attribute_traits()
#' -> signal_summary
#' -> derive_conditional_contrasts()
#' -> decompose_conditional_effects()
#' -> optional P2 modules
#' }
#'
#' Every downstream step uses the FINAL joint signal-specific effects;
#' stepwise temporary conditional betas never propagate.
#'
#' ## `signal_mode = "predefined"`
#'
#' Uses `control$predefined_signals` (a data.frame with `locus_id`,
#' `representative_snp` and `conditioning_snps`); the omnibus scan and
#' locus construction are skipped. The signal+trait oracle combines
#' `signal_mode = "predefined"`, `candidate_mode = "predefined"` and
#' `control$predefined_candidate_sets`.
#'
#' @param Y Numeric `n x m` phenotype matrix.
#' @param W Optional numeric `n x q` fixed-effect design matrix.
#' @param G Numeric `n x p` genotype dosage matrix coded 0/1/2.
#' @param K Optional numeric `n x n` genomic relationship matrix.
#' @param chromosome Optional length-`p` marker annotation.
#' @param position Optional length-`p` physical positions; required for
#'   `signal_mode = "resolve"`.
#' @param signal_mode `"resolve"` or `"predefined"`.
#' @param alpha_omnibus Omnibus significance level.
#' @param omnibus_adjust Adjustment method for the omnibus p-values.
#' @param candidate_mode Candidate-set mode: `"holm_fwer"`,
#'   `"all_traits"` or `"predefined"`.
#' @param alpha_trait Within-signal Holm level.
#' @param subset_mode Subset evaluation mode passed to
#'   [decompose_conditional_effects()].
#' @param custom_sets Optional list of trait sets for
#'   `subset_mode = "custom"`.
#' @param tolerance Primary representation-loss tolerance.
#' @param sensitivity_tolerance Additional tolerances.
#' @param max_traits_exact Maximum candidate-set size for exact
#'   enumeration.
#' @param pve,crossfit,bootstrap P2 switches (all `FALSE` by default;
#'   `TRUE` warns and returns `NULL` until Stage 10).
#' @param control List with optional fields `locus_method`,
#'   `window_bp`, `r2_threshold`, `merge_overlaps`, `signal_adjust`,
#'   `alpha_signal`, `fixed_p_threshold`, `max_signals`,
#'   `predefined_signals`, `predefined_candidate_sets`,
#'   `null_control`, `inverse_tol`, `qform_tol`,
#'   `decomposition_rel_tol`, `monotonicity_tol`,
#'   `condition_number_warning`. Unrecognised fields warn.
#' @param seed Random seed (recorded).
#'
#' @return An object of class `"condped_fit"` with components
#'   `null_fit`, `omnibus`, `loci`, `signals`, `effects`,
#'   `candidate_traits`, `signal_summary`, `subset_analysis`,
#'   `pve`, `crossfit`, `bootstrap`, `settings`, `status`,
#'   `diagnostics`.
#' @export
condped <- function(
  Y,
  W = NULL,
  G,
  K = NULL,
  chromosome = NULL,
  position = NULL,
  signal_mode = c("resolve", "predefined"),
  alpha_omnibus = 0.05,
  omnibus_adjust = c("BH", "bonferroni", "none"),
  candidate_mode = c("holm_fwer", "all_traits", "predefined"),
  alpha_trait = 0.05,
  subset_mode = c("all", "singleton", "custom"),
  custom_sets = NULL,
  tolerance = 0.10,
  sensitivity_tolerance = c(0.05, 0.10, 0.20),
  max_traits_exact = 10L,
  pve = FALSE,
  crossfit = FALSE,
  bootstrap = FALSE,
  control = list(),
  seed = 1L
) {
  t0 <- proc.time()[["elapsed"]]
  signal_mode <- match.arg(signal_mode)
  omnibus_adjust <- match.arg(omnibus_adjust)
  candidate_mode <- match.arg(candidate_mode)
  subset_mode <- match.arg(subset_mode)
  .check_prob(alpha_omnibus, "alpha_omnibus")
  .check_prob(alpha_trait, "alpha_trait")

  control_allowed <- c("locus_method", "window_bp", "r2_threshold",
                       "merge_overlaps", "signal_adjust", "alpha_signal",
                       "fixed_p_threshold", "max_signals",
                       "predefined_signals", "predefined_candidate_sets",
                       "null_control", "inverse_tol", "qform_tol",
                       "decomposition_rel_tol", "monotonicity_tol",
                       "condition_number_warning")
  if (!is.list(control)) {
    .stop_invalid_input("control must be a list.")
  }
  unknown_control <- setdiff(names(control), control_allowed)
  warnings <- character()
  if (length(unknown_control) > 0L) {
    msg <- sprintf("Unrecognised control field(s) ignored: %s.",
                   paste(unknown_control, collapse = ", "))
    warnings <- c(warnings, msg)
    warning(msg, call. = FALSE)
  }
  ctrl <- function(name, default) {
    if (is.null(control[[name]])) default else control[[name]]
  }
  null_control <- if (is.list(control$null_control)) {
    control$null_control
  } else {
    list()
  }

  # ---- input validation -------------------------------------------------------
  .validate_dimensions(Y = Y, W = W, G = G, K = K)
  if (anyNA(Y)) .stop_invalid_input("Y must not contain missing values.")
  marker_ids <- colnames(G)
  if (is.null(marker_ids)) marker_ids <- paste0("marker", seq_len(ncol(G)))
  if (anyDuplicated(marker_ids)) {
    .stop_invalid_input("marker ids (column names of G) must be unique.")
  }
  colnames(G) <- marker_ids
  if (is.null(chromosome)) chromosome <- rep("chr1", ncol(G))
  if (length(chromosome) != ncol(G)) {
    .stop_invalid_input("chromosome must have length ncol(G).")
  }
  if (signal_mode == "resolve") {
    if (is.null(position)) {
      .stop_invalid_input(
        "position is required for signal_mode = \"resolve\"."
      )
    }
    if (length(position) != ncol(G) || anyNA(position)) {
      .stop_invalid_input("position must have one value per marker.")
    }
  }
  if (candidate_mode == "predefined" &&
      is.null(control$predefined_candidate_sets)) {
    .stop_invalid_input(
      "control$predefined_candidate_sets is required for candidate_mode = \"predefined\"."
    )
  }
  if (!is.null(seed)) set.seed(seed)

  # ---- one fixed K, one null model ---------------------------------------------
  if (is.null(K)) {
    maf_hat <- colMeans(G, na.rm = TRUE) / 2
    sd_hat <- sqrt(2 * maf_hat * (1 - maf_hat))
    keep <- sd_hat > 0
    Z <- sweep(G[, keep, drop = FALSE], 2L, 2 * maf_hat[keep], `-`)
    Z <- sweep(Z, 2L, sd_hat[keep], `/`)
    K <- .make_grm(Z)
  }
  fit <- fit_mt_null(Y = Y, W = W, K = K, control = null_control)
  if (!isTRUE(fit$status$ok)) {
    return(.condped_failed_fit(fit, warnings, t0, signal_mode,
                               alpha_omnibus, omnibus_adjust,
                               candidate_mode, alpha_trait, subset_mode,
                               custom_sets, tolerance,
                               sensitivity_tolerance, max_traits_exact,
                               pve, crossfit, bootstrap, control, seed))
  }

  basis <- derive_conditional_contrasts(fit$Sigma_P_ref)

  empty_long <- data.frame(
    locus_id = character(), signal_id = character(),
    representative_snp = character(), trait = character(),
    beta = numeric(), se = numeric(), p_value = numeric(),
    stringsAsFactors = FALSE
  )

  if (signal_mode == "resolve") {
    # ---- omnibus scan + explicit genome-wide selection -------------------------
    scan <- scan_mt_omnibus(fit, G, marker_ids = marker_ids)
    om <- scan$omnibus
    valid <- !is.na(om$p_value)
    om$p_adjusted <- NA_real_
    om$p_adjusted[valid] <- switch(omnibus_adjust,
      BH = stats::p.adjust(om$p_value[valid], method = "BH"),
      bonferroni = stats::p.adjust(om$p_value[valid], method = "bonferroni"),
      none = om$p_value[valid]
    )
    om$selected <- !is.na(om$p_adjusted) & om$p_adjusted <= alpha_omnibus
    selected_markers <- om$marker_id[om$selected]
    scan$omnibus <- om

    if (length(selected_markers) == 0L) {
      loci_obj <- NULL
      resolution <- NULL
      effects_obj <- NULL
      attr <- attribute_traits(empty_long, candidate_mode = "all_traits")
      sig_summary <- .build_signal_summary(attr)
      dec <- decompose_conditional_effects(
        empty_long, attr, basis,
        restrict_to_candidates = TRUE,
        subset_mode = subset_mode, custom_sets = custom_sets,
        tolerance = tolerance,
        sensitivity_tolerance = sensitivity_tolerance,
        max_traits_exact = max_traits_exact
      )
      status_code <- "empty_selection"
    } else {
      loci_obj <- define_associated_loci(
        scan, G, chromosome, position,
        selected_markers = selected_markers,
        marker_ids = marker_ids,
        method = ctrl("locus_method", "physical"),
        window_bp = ctrl("window_bp", 1e6),
        r2_threshold = ctrl("r2_threshold", 0.2),
        locus_map = control$locus_map,
        merge_overlaps = ctrl("merge_overlaps", TRUE)
      )
      resolution <- resolve_locus_signals(
        fit, G, loci_obj,
        marker_ids = marker_ids,
        signal_adjust = ctrl("signal_adjust", "within_locus_bonferroni"),
        alpha_signal = ctrl("alpha_signal", 0.05),
        fixed_p_threshold = control$fixed_p_threshold,
        max_signals = ctrl("max_signals", 10L),
        rank_tol = ctrl("inverse_tol", sqrt(.Machine$double.eps))
      )
      effects_obj <- resolution
      attr <- attribute_traits(
        resolution, candidate_mode = candidate_mode,
        alpha_trait = alpha_trait,
        predefined_sets = control$predefined_candidate_sets
      )
      sig_summary <- .build_signal_summary(attr, resolution)
      dec <- decompose_conditional_effects(
        resolution, attr, basis,
        restrict_to_candidates = TRUE,   # candidate sets drive the scope
        subset_mode = subset_mode, custom_sets = custom_sets,
        tolerance = tolerance,
        sensitivity_tolerance = sensitivity_tolerance,
        max_traits_exact = max_traits_exact
      )
      status_code <- "ok"
    }
    loci_out <- if (is.null(loci_obj)) NULL else loci_obj$loci
    signals_out <- if (is.null(resolution)) NULL else resolution$signals
    omnibus_out <- scan
  } else {
    # ---- predefined signals ---------------------------------------------------
    ps <- control$predefined_signals
    if (is.null(ps) || !is.data.frame(ps) ||
        !all(c("locus_id", "representative_snp", "conditioning_snps") %in%
               names(ps))) {
      .stop_invalid_input(
        paste0("control$predefined_signals must be a data.frame with ",
               "locus_id, representative_snp and conditioning_snps.")
      )
    }
    missing_snps <- setdiff(ps$representative_snp, marker_ids)
    if (length(missing_snps) > 0L) {
      .stop_invalid_input(
        "predefined signal markers not found in G: %s.", missing_snps[1L]
      )
    }
    pre <- .effects_from_predefined(fit, G, ps, marker_ids,
                                    ctrl("inverse_tol",
                                         sqrt(.Machine$double.eps)))
    loci_out <- NULL
    signals_out <- pre$signals
    effects_obj <- pre
    omnibus_out <- NULL
    attr <- attribute_traits(
      pre, candidate_mode = candidate_mode,
      alpha_trait = alpha_trait,
      predefined_sets = control$predefined_candidate_sets
    )
    sig_summary <- .build_signal_summary(attr, resolution = NULL,
                                         locus_info = pre$locus_info)
    dec <- decompose_conditional_effects(
      pre, attr, basis,
      restrict_to_candidates = TRUE,
      subset_mode = subset_mode, custom_sets = custom_sets,
      tolerance = tolerance,
      sensitivity_tolerance = sensitivity_tolerance,
      max_traits_exact = max_traits_exact
    )
    status_code <- "ok"
  }

  # ---- P2 modules (never called when FALSE) --------------------------------------
  for (p2 in c("pve", "crossfit", "bootstrap")) {
    if (isTRUE(get(p2))) {
      msg <- sprintf(
        paste0("%s = TRUE requested, but the corresponding P2 module is ",
               "not implemented yet (Stage 10); returning NULL."), p2
      )
      warnings <- c(warnings, msg)
      warning(msg, call. = FALSE)
    }
  }

  structure(
    list(
      null_fit = fit,
      omnibus = omnibus_out,
      loci = loci_out,
      signals = signals_out,
      effects = effects_obj,
      candidate_traits = attr,
      signal_summary = sig_summary,
      subset_analysis = list(
        subset_table = dec$subset_table,
        minimum_representative_sets = dec$minimum_representative_sets,
        irreducible_modules = dec$irreducible_modules,
        tolerance_path = dec$tolerance_path
      ),
      pve = NULL,
      crossfit = NULL,
      bootstrap = NULL,
      settings = list(
        signal_mode = signal_mode,
        alpha_omnibus = alpha_omnibus,
        omnibus_adjust = omnibus_adjust,
        candidate_mode = candidate_mode,
        alpha_trait = alpha_trait,
        subset_mode = subset_mode,
        custom_sets = custom_sets,
        tolerance = tolerance,
        sensitivity_tolerance = sensitivity_tolerance,
        max_traits_exact = max_traits_exact,
        pve = pve, crossfit = crossfit, bootstrap = bootstrap,
        control = control, seed = seed
      ),
      status = .new_status(
        ok = TRUE, code = status_code,
        message = if (status_code == "empty_selection") {
          "No marker passed the omnibus screen; returning empty structures."
        } else {
          ""
        },
        warnings = warnings
      ),
      diagnostics = list(
        n_markers = ncol(G),
        n_selected_markers = if (signal_mode == "resolve") {
          length(selected_markers)
        } else {
          0L
        },
        n_loci = if (!is.null(loci_out)) nrow(loci_out) else 0L,
        n_signals = if (!is.null(signals_out)) nrow(signals_out) else 0L,
        decompose = dec$diagnostics,
        elapsed = proc.time()[["elapsed"]] - t0
      )
    ),
    class = "condped_fit"
  )
}

#' Effects from predefined signals (oracle path)
#'
#' Per signal s with conditioning set C_s, the joint model of
#' (C_s + s) is fitted through the shared path
#' [`.fit_joint_signal_effects()`] and the s block is returned. No
#' selection and no screening happen here.
#'
#' @param fit The single null fit.
#' @param G Genotype matrix.
#' @param ps data.frame with `locus_id`, `representative_snp`,
#'   `conditioning_snps` (character vector per row, or list-column).
#' @param marker_ids Marker order of `G`.
#' @param rank_tol Tolerance for [`.safe_inverse()`].
#' @return A resolve-like list with `signals`, `beta`, `covariance`,
#'   `locus_info`, `status`.
#' @keywords internal
.effects_from_predefined <- function(fit, G, ps, marker_ids, rank_tol) {
  signals_rows <- list()
  beta_rows <- list()
  cov_list <- list()
  locus_info_rows <- list()
  locus_ids <- unique(as.character(ps$locus_id))
  for (lid in locus_ids) {
    sub <- ps[as.character(ps$locus_id) == lid, , drop = FALSE]
    num_status <- "ok"
    for (k in seq_len(nrow(sub))) {
      snp <- as.character(sub$representative_snp[k])
      cond <- sub$conditioning_snps[[k]]
      if (is.null(cond)) cond <- character()
      cond <- setdiff(as.character(cond), snp)
      joint <- .fit_joint_signal_effects(
        fit, G, c(cond, snp), marker_ids = marker_ids, rank_tol = rank_tol
      )
      pos <- length(cond) + 1L
      sid <- if ("signal_id" %in% names(sub) &&
                 !is.na(sub$signal_id[k])) {
        as.character(sub$signal_id[k])
      } else {
        paste0(lid, "::S", k)
      }
      signals_rows[[length(signals_rows) + 1L]] <- data.frame(
        locus_id = lid, signal_order = k, representative_snp = snp,
        signal_id = sid, stringsAsFactors = FALSE
      )
      beta_rows[[length(beta_rows) + 1L]] <- data.frame(
        locus_id = lid, signal_order = k, representative_snp = snp,
        signal_id = sid, trait = colnames(joint$beta),
        beta = joint$beta[pos, ], se = joint$se[pos, ],
        stringsAsFactors = FALSE
      )
      cov_list[[sid]] <- joint$joint_covariance
      if (joint$status$code != "ok") num_status <- joint$status$code
    }
    locus_info_rows[[length(locus_info_rows) + 1L]] <- data.frame(
      locus_id = lid, n_signals = nrow(sub),
      numerical_status = num_status, stringsAsFactors = FALSE
    )
  }
  signals_df <- do.call(rbind, signals_rows)
  beta_df <- do.call(rbind, beta_rows)
  # .as_signal_effects() prefers a signal_id column when present
  list(
    signals = signals_df,
    beta = beta_df,
    covariance = cov_list,
    locus_info = do.call(rbind, locus_info_rows),
    status = .new_status(ok = TRUE, code = "ok"),
    diagnostics = list(n_loci = length(locus_ids),
                       n_signals = nrow(signals_df))
  )
}

#' Failed-fit early return for condped()
#' @keywords internal
.condped_failed_fit <- function(fit, warnings, t0, signal_mode,
                                alpha_omnibus, omnibus_adjust,
                                candidate_mode, alpha_trait, subset_mode,
                                custom_sets, tolerance,
                                sensitivity_tolerance, max_traits_exact,
                                pve, crossfit, bootstrap, control, seed) {
  structure(
    list(
      null_fit = fit, omnibus = NULL, loci = NULL, signals = NULL,
      effects = NULL, candidate_traits = NULL, signal_summary = NULL,
      subset_analysis = list(
        subset_table = NULL, minimum_representative_sets = NULL,
        irreducible_modules = NULL, tolerance_path = NULL
      ),
      pve = NULL, crossfit = NULL, bootstrap = NULL,
      settings = list(
        signal_mode = signal_mode, alpha_omnibus = alpha_omnibus,
        omnibus_adjust = omnibus_adjust, candidate_mode = candidate_mode,
        alpha_trait = alpha_trait, subset_mode = subset_mode,
        custom_sets = custom_sets, tolerance = tolerance,
        sensitivity_tolerance = sensitivity_tolerance,
        max_traits_exact = max_traits_exact,
        pve = pve, crossfit = crossfit, bootstrap = bootstrap,
        control = control, seed = seed
      ),
      status = .new_status(
        ok = FALSE, code = "non_convergence",
        message = "fit_mt_null() did not converge; downstream steps skipped.",
        warnings = warnings
      ),
      diagnostics = list(elapsed = proc.time()[["elapsed"]] - t0)
    ),
    class = "condped_fit"
  )
}
