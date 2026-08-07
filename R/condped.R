#' End-to-end CondPED analysis (v1.0)
#'
#' The thinnest possible orchestrator: it only calls the existing
#' public functions in the frozen order and never re-implements any
#' formula. `fit_mt_null()` is called exactly once; there is no LOCO
#' branch; `chromosome` is marker annotation only. When the P2 flags
#' are `FALSE` (default), the P2 functions are never called and the
#' corresponding fields are `NULL`.
#'
#' @details
#' Frozen pipeline:
#' \preformatted{
#' validate inputs
#' -> construct/read one fixed K
#' -> fit_mt_null() exactly once
#' -> scan_mt_omnibus()
#' -> determine follow-up loci
#' -> estimate_mt_effects()
#' -> attribute_traits()
#' -> derive_conditional_contrasts()
#' -> decompose_conditional_effects()
#' -> optional P2 modules
#' }
#'
#' Follow-up loci: `control$followup_loci` when given (validated
#' against the marker ids; the family for layer-1 adjustment then
#' consists of these loci), otherwise the loci passing the omnibus
#' screen (`omnibus_adjust` at `alpha_omnibus`). With no follow-up
#' loci the function returns a structurally complete empty result with
#' status code `empty_selection`; `omnibus_only` and
#' `trait_restricted` loci are legal results, not failures.
#'
#' @param Y Numeric `n x m` phenotype matrix (rows are individuals).
#' @param W Optional numeric `n x q` fixed-effect design matrix; when
#'   `NULL` an intercept-only design is used.
#' @param G Numeric `n x p` genotype dosage matrix coded 0/1/2.
#' @param K Optional numeric `n x n` genomic relationship matrix; when
#'   `NULL` it is built once from `G`.
#' @param chromosome Optional length-`p` marker annotation vector; only
#'   carried through to `settings`, never used to switch null models.
#' @param alpha_omnibus Omnibus significance level.
#' @param omnibus_adjust Adjustment method for the omnibus p-values.
#' @param candidate_mode Candidate-set mode passed to
#'   [attribute_traits()]: `"holm_fwer"` (formal) or `"all_traits"`.
#' @param alpha_trait Within-locus Holm level.
#' @param subset_mode Subset evaluation mode passed to
#'   [decompose_conditional_effects()].
#' @param custom_sets Optional list of trait sets for
#'   `subset_mode = "custom"`.
#' @param tolerance Primary representation-loss tolerance.
#' @param sensitivity_tolerance Additional tolerances for the
#'   sensitivity path.
#' @param max_traits_exact Maximum candidate-set size for exact
#'   enumeration.
#' @param pve,crossfit,bootstrap Logical P2 switches; all `FALSE` by
#'   default. The P2 modules are scheduled for a later stage: when set
#'   to `TRUE` a warning is emitted and the field stays `NULL`.
#' @param control List with optional fields `followup_loci`,
#'   `locus_unit`, `inverse_tol`, `qform_tol`,
#'   `decomposition_rel_tol`, `monotonicity_tol`,
#'   `condition_number_warning`, `null_control`. Unrecognised fields
#'   produce a warning.
#' @param seed Random seed (recorded; the pipeline itself is
#'   deterministic given the data).
#'
#' @return An object of class `"condped_fit"`: a list with components
#'   `null_fit`, `omnibus`, `effects`, `candidate_traits`,
#'   `subset_analysis` (list with `subset_table`,
#'   `minimum_representative_sets`, `irreducible_modules`,
#'   `tolerance_path`), `pve`, `crossfit`, `bootstrap`, `settings`,
#'   `status` and `diagnostics`.
#' @export
condped <- function(
  Y,
  W = NULL,
  G,
  K = NULL,
  chromosome = NULL,
  alpha_omnibus = 0.05,
  omnibus_adjust = c("BH", "bonferroni", "none"),
  candidate_mode = c("holm_fwer", "all_traits"),
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
  omnibus_adjust <- match.arg(omnibus_adjust)
  candidate_mode <- match.arg(candidate_mode)
  subset_mode <- match.arg(subset_mode)
  .check_prob(alpha_omnibus, "alpha_omnibus")
  .check_prob(alpha_trait, "alpha_trait")

  # ---- control ---------------------------------------------------------------
  control_allowed <- c("followup_loci", "locus_unit", "inverse_tol",
                       "qform_tol", "decomposition_rel_tol",
                       "monotonicity_tol", "condition_number_warning",
                       "null_control")
  if (!is.list(control)) {
    .stop_invalid_input("control must be a list.")
  }
  unknown_control <- setdiff(names(control), control_allowed)
  warnings <- character()
  if (length(unknown_control) > 0L) {
    msg <- sprintf(
      "Unrecognised control field(s) ignored: %s.",
      paste(unknown_control, collapse = ", ")
    )
    warnings <- c(warnings, msg)
    warning(msg, call. = FALSE)
  }
  followup_loci <- control$followup_loci
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
  if (!is.null(chromosome) && length(chromosome) != ncol(G)) {
    .stop_invalid_input("chromosome must have length ncol(G).")
  }
  if (!is.null(followup_loci)) {
    followup_loci <- as.character(followup_loci)
    missing_fu <- setdiff(followup_loci, marker_ids)
    if (length(missing_fu) > 0L) {
      .stop_invalid_input(
        "control$followup_loci not found in G: %s.",
        paste(utils::head(missing_fu, 3L), collapse = ", ")
      )
    }
    followup_loci <- unique(followup_loci)
  }
  if (!is.null(seed)) set.seed(seed)

  # ---- one fixed K -------------------------------------------------------------
  if (is.null(K)) {
    maf_hat <- colMeans(G, na.rm = TRUE) / 2
    sd_hat <- sqrt(2 * maf_hat * (1 - maf_hat))
    keep <- sd_hat > 0
    Z <- sweep(G[, keep, drop = FALSE], 2L, 2 * maf_hat[keep], `-`)
    Z <- sweep(Z, 2L, sd_hat[keep], `/`)
    K <- .make_grm(Z)
  }

  # ---- one null model ------------------------------------------------------------
  fit <- fit_mt_null(Y = Y, W = W, K = K, control = null_control)
  if (!isTRUE(fit$status$ok)) {
    return(structure(
      list(
        null_fit = fit, omnibus = NULL, effects = NULL,
        candidate_traits = NULL,
        subset_analysis = list(
          subset_table = NULL, minimum_representative_sets = NULL,
          irreducible_modules = NULL, tolerance_path = NULL
        ),
        pve = NULL, crossfit = NULL, bootstrap = NULL,
        settings = .condped_settings(
          alpha_omnibus, omnibus_adjust, candidate_mode, alpha_trait,
          subset_mode, custom_sets, tolerance, sensitivity_tolerance,
          max_traits_exact, pve, crossfit, bootstrap,
          followup_loci, chromosome, seed
        ),
        status = .new_status(
          ok = FALSE, code = "non_convergence",
          message = "fit_mt_null() did not converge; downstream steps skipped.",
          warnings = warnings
        ),
        diagnostics = list(elapsed = proc.time()[["elapsed"]] - t0)
      ),
      class = "condped_fit"
    ))
  }

  # ---- omnibus scan ---------------------------------------------------------------
  scan <- scan_mt_omnibus(fit, G, marker_ids = marker_ids)
  om <- scan$omnibus

  # ---- follow-up loci --------------------------------------------------------------
  valid <- !is.na(om$p_value)
  p_adj <- switch(omnibus_adjust,
    BH = stats::p.adjust(om$p_value[valid], method = "BH"),
    bonferroni = stats::p.adjust(om$p_value[valid], method = "bonferroni"),
    none = om$p_value[valid]
  )
  omnibus_hits <- om$marker_id[valid][p_adj <= alpha_omnibus]
  if (is.null(followup_loci)) {
    followup <- omnibus_hits
    om_for_attr <- om
  } else {
    followup <- followup_loci
    om_for_attr <- om[om$marker_id %in% followup, , drop = FALSE]
  }
  n_followup_significant <- sum(followup %in% omnibus_hits)

  # ---- effects, candidates, subset analysis ------------------------------------------
  if (length(followup) > 0L) {
    est <- estimate_mt_effects(fit, G, loci = followup,
                               marker_ids = marker_ids)
    attr <- attribute_traits(om_for_attr, est,
                             omnibus_method = omnibus_adjust,
                             alpha_omnibus = alpha_omnibus,
                             candidate_mode = candidate_mode,
                             alpha_trait = alpha_trait)
    basis <- derive_conditional_contrasts(fit$Sigma_P_ref)
    dec <- decompose_conditional_effects(
      est, attr, basis,
      restrict_to_candidates = TRUE,
      subset_mode = subset_mode, custom_sets = custom_sets,
      tolerance = tolerance,
      sensitivity_tolerance = sensitivity_tolerance,
      max_traits_exact = max_traits_exact
    )
    effects_out <- est
    candidate_out <- attr
  } else {
    # structurally complete empty result; estimate_mt_effects() is not
    # called on an empty locus set
    empty_long <- data.frame(
      marker_id = character(), trait = character(),
      beta = numeric(), se = numeric(), z = numeric(),
      p_value = numeric(), stringsAsFactors = FALSE
    )
    attr <- attribute_traits(om_for_attr, empty_long,
                             omnibus_method = omnibus_adjust,
                             alpha_omnibus = alpha_omnibus,
                             candidate_mode = candidate_mode,
                             alpha_trait = alpha_trait)
    basis <- derive_conditional_contrasts(fit$Sigma_P_ref)
    dec <- decompose_conditional_effects(
      empty_long, attr, basis,
      restrict_to_candidates = TRUE,
      subset_mode = subset_mode, custom_sets = custom_sets,
      tolerance = tolerance,
      sensitivity_tolerance = sensitivity_tolerance,
      max_traits_exact = max_traits_exact
    )
    effects_out <- NULL
    candidate_out <- attr
  }

  # ---- P2 modules (off by default; never called when FALSE) ---------------------------
  pve_res <- crossfit_res <- bootstrap_res <- NULL
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

  status_code <- if (length(followup) == 0L) "empty_selection" else "ok"
  structure(
    list(
      null_fit = fit,
      omnibus = scan,
      effects = effects_out,
      candidate_traits = candidate_out,
      subset_analysis = list(
        subset_table = dec$subset_table,
        minimum_representative_sets = dec$minimum_representative_sets,
        irreducible_modules = dec$irreducible_modules,
        tolerance_path = dec$tolerance_path
      ),
      pve = pve_res,
      crossfit = crossfit_res,
      bootstrap = bootstrap_res,
      settings = .condped_settings(
        alpha_omnibus, omnibus_adjust, candidate_mode, alpha_trait,
        subset_mode, custom_sets, tolerance, sensitivity_tolerance,
        max_traits_exact, pve, crossfit, bootstrap,
        followup_loci, chromosome, seed
      ),
      status = .new_status(
        ok = TRUE, code = status_code,
        message = if (status_code == "empty_selection") {
          "No locus passed the omnibus screen; returning empty structures."
        } else {
          ""
        },
        warnings = warnings
      ),
      diagnostics = list(
        n_markers = ncol(G),
        n_omnibus_significant = length(omnibus_hits),
        n_followup = length(followup),
        n_followup_significant = n_followup_significant,
        n_selected_candidates = attr$diagnostics$n_selected,
        locus_status_counts = attr$diagnostics[c("n_omnibus_only",
                                                 "n_trait_restricted",
                                                 "n_multi_trait")],
        decompose = dec$diagnostics,
        elapsed = proc.time()[["elapsed"]] - t0
      )
    ),
    class = "condped_fit"
  )
}

#' Collect condped() settings for the return object
#' @keywords internal
.condped_settings <- function(alpha_omnibus, omnibus_adjust,
                              candidate_mode, alpha_trait,
                              subset_mode, custom_sets, tolerance,
                              sensitivity_tolerance, max_traits_exact,
                              pve, crossfit, bootstrap,
                              followup_loci, chromosome, seed) {
  list(
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
    followup_loci = followup_loci,
    chromosome = chromosome,
    seed = seed
  )
}
