#' Evaluate CondPED simulation replicates (v1.0, Stage 6C)
#'
#' Computes per-replicate recovery metrics and grouped summaries with
#' Monte Carlo standard errors. Truth signals and estimated signals are
#' first matched by genomic position (never by signal-id string): loci
#' are matched by chromosome plus interval overlap, and signals within
#' overlapping loci are matched one-to-one by representative-to-causal
#' distance (or by genotype r^2 when a genotype matrix is supplied).
#' All per-signal metrics are averages over the MATCHED signal pairs;
#' unrecovered signals are absent from beta/candidate/eta/rho/Rep/Irr
#' metrics (never coded as beta = 0). All set comparisons use
#' normalised trait keys; families are compared as whole sets of keys
#' (never just the first solution). The empty-set rules of the frozen
#' contract are applied verbatim.
#'
#' @details
#' ## Replicate file format
#'
#' Each entry of `files` is an RDS containing one replicate result:
#' \describe{
#'   \item{truth}{The `truth` component of a
#'     [simulate_condped_data()] call: `loci`, `signals`, `beta` (long
#'     format), `candidate_traits`, `effect_breadth`,
#'     `effect_direction`, `minimum_representative_sets`,
#'     `irreducible_modules`, `subset_table` and friends.}
#'   \item{estimates}{A list with `signals` (data.frame with
#'     `locus_id`, `signal_id`, `signal_order`, `representative_snp`
#'     and optionally `chromosome`/`position`), `beta` (long
#'     data.frame with `locus_id`, `signal_order`, `trait`, `beta`,
#'     `se`), `candidate_sets` (named list signal_id -> traits),
#'     `minimum_representative_sets`, `irreducible_modules`,
#'     `subset_table`, `omnibus_summary` (data.frame with `marker_id`
#'     and `p_value`) and optionally `positions` (named marker ->
#'     position vector), `G` (genotype matrix) and `asset` (per-pipeline
#'     ASSET results, see below).}
#'   \item{settings}{List with the scenario descriptors named in
#'     `group_by` plus `tolerance` and `alpha_omnibus`; may also carry
#'     `positions` and `G`.}
#'   \item{status}{Standard status list.}
#'   \item{runtime}{Numeric seconds.}
#' }
#'
#' The optional `estimates$asset` component is a list with entries
#' `lead`, `resolved` and `condped`; each entry is a list with `sets`
#' (named list signal_id -> selected trait subset) and `p` (named
#' p-values), and optionally `direction` (named direction patterns).
#' When absent, the nine ASSET/ablation metrics are `NA`.
#'
#' Rep metrics are the PRIMARY set-level summary; Irr metrics are the
#' secondary set-level summary (recorded in `settings$metric_levels`).
#'
#' ## Metric definitions (per replicate, primary tolerance)
#'
#' Signal-resolution metrics:
#' \describe{
#'   \item{locus_detection_rate}{Fraction of truth loci overlapped by
#'     at least one estimated locus.}
#'   \item{signal_count_exact_recovery}{Over detected truth loci,
#'     fraction whose estimated signal count equals the truth signal
#'     count.}
#'   \item{secondary_signal_power}{Fraction of truth signals with
#'     `signal_order >= 2` that are matched (`NA` when the truth has a
#'     single signal).}
#'   \item{missed_signal_rate}{Missed truth signals / total truth
#'     signals.}
#'   \item{extra_signal_rate}{Extra estimated signals /
#'     max(1, total estimated signals).}
#'   \item{signal_coverage}{Matched truth signals / total truth
#'     signals.}
#'   \item{representative_to_causal_r2}{Mean squared genotype
#'     correlation between estimated representative and causal marker
#'     over matched pairs (`NA` without a genotype matrix).}
#'   \item{representative_to_causal_distance}{Mean distance in bp over
#'     matched pairs.}
#'   \item{residual_local_false_positive_rate}{Fraction of detected
#'     truth loci containing at least one unmatched (extra) estimated
#'     signal.}
#'   \item{signal_resolution_attempt_rate}{Indicator that at least one
#'     truth locus was detected and therefore entered conditional signal
#'     resolution (`NA` when the truth has no loci).}
#'   \item{conditional_signal_count_exact_recovery}{Alias of
#'     `signal_count_exact_recovery`, making the conditioning on locus
#'     detection explicit (`NA` when no truth locus is detected).}
#'   \item{conditional_secondary_signal_power}{Fraction of truth
#'     secondary signals (`signal_order >= 2`) in detected truth loci
#'     that are matched (`NA` when no detected locus carries a
#'     secondary signal).}
#'   \item{conditional_missed_signal_rate}{Missed truth signals in
#'     detected truth loci / truth signals in detected loci (`NA` when
#'     no truth locus is detected).}
#'   \item{truth_locus_split_rate}{Fraction of detected truth loci
#'     overlapped by more than one estimated locus (`NA` when no truth
#'     locus is detected).}
#' }
#'
#' Quantitative metrics (averaged over matched signal pairs):
#' \describe{
#'   \item{type1_omnibus}{Fraction of non-causal markers with
#'     `p_value <= alpha_omnibus`.}
#'   \item{power_omnibus}{Indicator that any causal marker passes;
#'     `NA` for the null architecture.}
#'   \item{beta_bias, beta_rmse, beta_coverage,
#'     beta_sign_accuracy}{Signed mean error, RMSE, 95% interval
#'     coverage (only when `se` is supplied) and sign accuracy (over
#'     traits with nonzero true effect), computed per matched pair
#'     over traits and then averaged.}
#'   \item{candidate_tpr, candidate_fdp, candidate_precision,
#'     candidate_exact_recovery, candidate_jaccard}{TPR / FDP /
#'     precision (1 - FDP) / exact recovery / Jaccard of the estimated
#'     candidate set against `truth$candidate_traits`, per matched
#'     pair.}
#'   \item{effect_breadth_bias}{Estimated candidate-set size minus
#'     truth breadth, averaged over matched pairs.}
#'   \item{direction_recovery}{Fraction of matched pairs whose
#'     estimated direction pattern (frozen rule: 0 candidates ->
#'     not_applicable, 1 -> single_trait, >= 2 same sign ->
#'     concordant, 2 opposite -> antagonistic, >= 3 both signs ->
#'     mixed) equals `truth$effect_direction`.}
#'   \item{conditional_effect_bias, conditional_effect_rmse,
#'     conditional_effect_sign_accuracy}{Elementwise bias, RMSE and
#'     sign accuracy (over elements with nonzero true eta) of the
#'     conditional effects eta over subset-table rows matched by
#'     `representing_key` within the matched pair.}
#'   \item{representation_loss_bias, representation_loss_rmse,
#'     representation_loss_mae, representation_threshold_accuracy,
#'     representation_map_mae}{Bias, RMSE and MAE of the
#'     representation loss over matched subset rows; the fraction of
#'     rows whose feasibility side (loss <= tolerance) matches the
#'     truth; and the row-wise MAE of the subset map.}
#'   \item{rep_min_cardinality_recovery, rep_exact_family_recovery,
#'     rep_tie_coverage}{Indicator that the estimated minimum set size
#'     equals the truth minimum set size; representative-family exact
#'     recovery; and the fraction of truth tied minimum sets whose key
#'     appears in the estimated family.}
#'   \item{irr_family_recovery, irr_missed_rate,
#'     irr_extra_split_rate}{Irreducible-module family exact recovery;
#'     fraction of truth joint modules (size > 1) whose key is not
#'     exactly recovered; fraction of truth joint modules split into
#'     two or more strict sub-modules (both `NA` when the truth has no
#'     joint module).}
#'   \item{asset_power, asset_trait_set_jaccard,
#'     asset_trait_set_precision, asset_trait_set_recall,
#'     asset_direction_recovery, lead_trait_breadth_inflation,
#'     resolved_trait_breadth_inflation,
#'     pseudo_multitrait_interpretation_rate}{ASSET/ablation metrics
#'     read from `estimates$asset`; all `NA` when absent.}
#'   \item{overall_recovery}{Candidate exact AND representative family
#'     exact AND module family exact, averaged over matched pairs.}
#'   \item{conditional_recovery}{Representative family exact AND
#'     module family exact, evaluated only for pairs whose candidate
#'     set was exactly recovered (`NA` otherwise).}
#' }
#'
#' Empty-set rules (frozen): truth empty + estimate empty -> exact = 1,
#' Jaccard = 1, TPR = NA, FDP = 0; truth non-empty + estimate empty ->
#' TPR = 0, FDP = 0, Jaccard = 0, exact = 0; truth empty + estimate
#' non-empty -> TPR = NA, FDP = 1, Jaccard = 0, exact = 0.
#'
#' @param files Character vector of RDS file paths, one per replicate.
#' @param metrics Metrics to report; defaults to the full frozen list.
#' @param include_unstable Logical; keep replicates whose status code
#'   is `unstable` in the metric summaries (they are always listed in
#'   `failures`).
#' @param group_by Scenario columns to group summaries by; columns
#'   absent from the replicate settings are dropped with a warning.
#' @param conf_level Confidence level for the summary intervals.
#'
#' @return A list with components `replicate_metrics` (one row per
#'   evaluated replicate), `summary` (grouped means), `mcse` (grouped
#'   Monte Carlo standard errors), `failures` (read failures and
#'   non-ok replicates), `settings` and `status`.
#' @export
evaluate_condped_simulation <- function(
  files,
  metrics = c(
    "locus_detection_rate",
    "signal_count_exact_recovery",
    "secondary_signal_power",
    "missed_signal_rate",
    "extra_signal_rate",
    "signal_coverage",
    "representative_to_causal_r2",
    "representative_to_causal_distance",
    "residual_local_false_positive_rate",
    "signal_resolution_attempt_rate",
    "conditional_signal_count_exact_recovery",
    "conditional_secondary_signal_power",
    "conditional_missed_signal_rate",
    "truth_locus_split_rate",
    "type1_omnibus",
    "power_omnibus",
    "beta_bias",
    "beta_rmse",
    "beta_coverage",
    "beta_sign_accuracy",
    "candidate_tpr",
    "candidate_fdp",
    "candidate_precision",
    "candidate_exact_recovery",
    "candidate_jaccard",
    "effect_breadth_bias",
    "direction_recovery",
    "conditional_effect_bias",
    "conditional_effect_rmse",
    "conditional_effect_sign_accuracy",
    "representation_loss_bias",
    "representation_loss_rmse",
    "representation_loss_mae",
    "representation_threshold_accuracy",
    "representation_map_mae",
    "rep_min_cardinality_recovery",
    "rep_exact_family_recovery",
    "rep_tie_coverage",
    "irr_family_recovery",
    "irr_missed_rate",
    "irr_extra_split_rate",
    "asset_power",
    "asset_trait_set_jaccard",
    "asset_trait_set_precision",
    "asset_trait_set_recall",
    "asset_direction_recovery",
    "lead_trait_breadth_inflation",
    "resolved_trait_breadth_inflation",
    "pseudo_multitrait_interpretation_rate",
    "overall_recovery",
    "conditional_recovery",
    "runtime"
  ),
  include_unstable = TRUE,
  group_by = c(
    "experiment",
    "n",
    "architecture",
    "locus_pve",
    "correlation",
    "target_loss",
    "tolerance"
  ),
  conf_level = 0.95
) {
  metric_all <- c(
    "locus_detection_rate", "signal_count_exact_recovery",
    "secondary_signal_power", "missed_signal_rate", "extra_signal_rate",
    "signal_coverage", "representative_to_causal_r2",
    "representative_to_causal_distance",
    "residual_local_false_positive_rate",
    "signal_resolution_attempt_rate",
    "conditional_signal_count_exact_recovery",
    "conditional_secondary_signal_power",
    "conditional_missed_signal_rate",
    "truth_locus_split_rate",
    "type1_omnibus",
    "power_omnibus", "beta_bias", "beta_rmse", "beta_coverage",
    "beta_sign_accuracy", "candidate_tpr", "candidate_fdp",
    "candidate_precision", "candidate_exact_recovery",
    "candidate_jaccard", "effect_breadth_bias", "direction_recovery",
    "conditional_effect_bias", "conditional_effect_rmse",
    "conditional_effect_sign_accuracy", "representation_loss_bias",
    "representation_loss_rmse", "representation_loss_mae",
    "representation_threshold_accuracy", "representation_map_mae",
    "rep_min_cardinality_recovery", "rep_exact_family_recovery",
    "rep_tie_coverage", "irr_family_recovery", "irr_missed_rate",
    "irr_extra_split_rate", "asset_power", "asset_trait_set_jaccard",
    "asset_trait_set_precision", "asset_trait_set_recall",
    "asset_direction_recovery", "lead_trait_breadth_inflation",
    "resolved_trait_breadth_inflation",
    "pseudo_multitrait_interpretation_rate", "overall_recovery",
    "conditional_recovery", "runtime"
  )
  metrics <- match.arg(metrics, several.ok = TRUE)
  if (!is.logical(include_unstable) || length(include_unstable) != 1L ||
      is.na(include_unstable)) {
    .stop_invalid_input("include_unstable must be TRUE or FALSE.")
  }
  .check_prob(conf_level, "conf_level")
  if (!is.character(files) || length(files) == 0L) {
    .stop_invalid_input("files must be a non-empty character vector of paths.")
  }

  warnings <- character()
  rep_rows <- list()
  fail_rows <- list()

  for (i in seq_along(files)) {
    f <- files[[i]]
    obj <- tryCatch(readRDS(f), error = function(e) e)
    if (inherits(obj, "error")) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        file = f, reason = paste("read error:", conditionMessage(obj)),
        included = FALSE, stringsAsFactors = FALSE
      )
      next
    }
    if (is.null(obj$truth) || is.null(obj$estimates)) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        file = f, reason = "missing truth or estimates component",
        included = FALSE, stringsAsFactors = FALSE
      )
      next
    }
    ok <- isTRUE(obj$status$ok)
    code <- if (!is.null(obj$status$code)) obj$status$code else "ok"
    unstable <- !ok || code != "ok"
    included <- !unstable || isTRUE(include_unstable)
    if (unstable) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        file = f, reason = sprintf("replicate status: %s", code),
        included = included, stringsAsFactors = FALSE
      )
    }
    if (!included) next
    G <- obj$G
    if (is.null(G) && !is.null(obj$estimates$G)) G <- obj$estimates$G
    if (is.null(G) && !is.null(obj$settings$G)) G <- obj$settings$G
    row <- tryCatch(
      .evaluate_one_replicate(obj$truth, obj$estimates,
                              settings = obj$settings,
                              runtime = obj$runtime, G = G),
      error = function(e) e
    )
    if (inherits(row, "error")) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        file = f,
        reason = paste("evaluation error:", conditionMessage(row)),
        included = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }
    row$file <- f
    rep_rows[[length(rep_rows) + 1L]] <- row
  }

  replicate_metrics <- if (length(rep_rows) > 0L) {
    out <- do.call(rbind, rep_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame()
  }
  failures <- if (length(fail_rows) > 0L) {
    out <- do.call(rbind, fail_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame(file = character(), reason = character(),
               included = logical(), stringsAsFactors = FALSE)
  }

  # ---- grouped summary -------------------------------------------------------
  # A group_by column is honoured only when at least one replicate
  # actually carries it; all-missing columns are dropped with a warning.
  present <- intersect(group_by, names(replicate_metrics))
  informative <- present[vapply(present, function(g) {
    nrow(replicate_metrics) > 0L && !all(is.na(replicate_metrics[[g]]))
  }, logical(1))]
  dropped <- setdiff(group_by, informative)
  if (length(dropped) > 0L && nrow(replicate_metrics) > 0L) {
    msg <- sprintf(
      "group_by column(s) not present in replicate settings and dropped: %s.",
      paste(dropped, collapse = ", ")
    )
    warnings <- c(warnings, msg)
    warning(msg, call. = FALSE)
  }
  group_cols <- informative
  summary_df <- data.frame()
  mcse_df <- data.frame()
  if (nrow(replicate_metrics) > 0L) {
    summ <- .summarise_metrics(replicate_metrics, metrics, group_cols,
                               conf_level)
    summary_df <- summ$summary
    mcse_df <- summ$mcse
  }

  list(
    replicate_metrics = replicate_metrics,
    summary = summary_df,
    mcse = mcse_df,
    failures = failures,
    settings = list(
      metrics = metrics,
      metric_levels = list(
        primary = c("rep_min_cardinality_recovery",
                    "rep_exact_family_recovery",
                    "rep_tie_coverage"),
        secondary = c("irr_family_recovery",
                      "irr_missed_rate",
                      "irr_extra_split_rate")
      ),
      include_unstable = include_unstable,
      group_by = group_by,
      conf_level = conf_level,
      n_files = length(files)
    ),
    status = .new_status(
      ok = TRUE,
      code = if (length(rep_rows) == 0L) "empty_selection" else "ok",
      message = if (length(rep_rows) == 0L) {
        "No evaluable replicates."
      } else {
        ""
      },
      warnings = warnings
    )
  )
}

#' Per-replicate metrics against the population truth
#'
#' Matches truth and estimated signals with [.match_signals()] and
#' averages the per-signal metrics over the matched pairs.
#'
#' @param truth `simulate_condped_data()` truth component.
#' @param estimates Estimates list (see
#'   [evaluate_condped_simulation()]).
#' @param settings Replicate settings (scenario descriptors).
#' @param runtime Numeric seconds.
#' @param G Optional n x p genotype matrix with marker-id column
#'   names; enables the representative-to-causal r^2 metric and
#'   r^2-based signal matching.
#' @return A one-row data.frame of metrics plus scenario columns.
#' @keywords internal
.evaluate_one_replicate <- function(truth, estimates,
                                    settings = list(),
                                    runtime = NA_real_,
                                    G = NULL) {
  tol <- if (!is.null(settings$tolerance)) settings$tolerance else 0.10
  alpha <- if (!is.null(settings$alpha_omnibus)) {
    settings$alpha_omnibus
  } else {
    0.05
  }
  if (is.null(estimates$positions) && !is.null(settings$positions)) {
    estimates$positions <- settings$positions
  }

  # ---- key uniqueness (never silently aggregate) ------------------------------
  .check_key_uniqueness(truth$signals, c("locus_id", "signal_id"),
                        "truth$signals")
  .check_key_uniqueness(estimates$signals, c("locus_id", "signal_id"),
                        "estimates$signals")
  .check_key_uniqueness(truth$beta, c("locus_id", "signal_id", "trait"),
                        "truth$beta")
  eb_sig_key <- if (!is.null(estimates$beta) &&
                    !("signal_id" %in% names(estimates$beta)) &&
                    "signal_order" %in% names(estimates$beta)) {
    "signal_order"
  } else {
    "signal_id"
  }
  .check_key_uniqueness(estimates$beta,
                        c("locus_id", eb_sig_key, "trait"),
                        "estimates$beta")
  for (side in c("truth", "estimates")) {
    obj <- if (side == "truth") truth else estimates
    sub <- obj$subset_table
    .check_key_uniqueness(
      sub,
      c("locus_id", "signal_id", "representing_key",
        intersect(c("complement_key", "complement_trait"), names(sub))),
      paste0(side, "$subset_table")
    )
    for (nm in c("minimum_representative_sets",
                 "irreducible_modules", "tolerance_path")) {
      .check_key_uniqueness(obj[[nm]],
                            c("locus_id", "signal_id", "tolerance",
                              "trait_key"),
                            paste0(side, "$", nm))
    }
    # the representation map is keyed per subset row, not per tolerance
    .check_key_uniqueness(
      obj$representation_map,
      c("locus_id", "signal_id", "representing_key"),
      paste0(side, "$representation_map")
    )
  }

  # ---- signal matching ---------------------------------------------------------
  m <- .match_signals(truth, estimates, G = G)
  t_sig <- m$truth_signals
  n_t <- nrow(t_sig)
  n_e <- nrow(m$est_signals)
  n_loci_t <- nrow(m$truth_loci)
  matches <- m$matches
  np <- nrow(matches)

  # ---- signal-resolution metrics -----------------------------------------------
  detected_loci <- unique(m$locus_matches$truth_locus_id)
  locus_detection_rate <- if (n_loci_t == 0L) {
    NA_real_
  } else {
    m$matched_locus_count / n_loci_t
  }

  count_exact <- NA_real_
  local_fpr <- NA_real_
  split_rate <- NA_real_
  if (length(detected_loci) > 0L) {
    exact_v <- numeric(length(detected_loci))
    fpr_v <- numeric(length(detected_loci))
    n_est_per_locus <- integer(length(detected_loci))
    for (li in seq_along(detected_loci)) {
      lid <- detected_loci[[li]]
      e_loci <- unique(m$locus_matches$est_locus_id[
        m$locus_matches$truth_locus_id == lid
      ])
      n_est_per_locus[[li]] <- length(e_loci)
      e_ids <- m$est_signals$signal_id[m$est_signals$locus_id %in% e_loci]
      exact_v[[li]] <- as.numeric(
        length(e_ids) == sum(t_sig$locus_id == lid)
      )
      fpr_v[[li]] <- as.numeric(any(e_ids %in% m$extra))
    }
    count_exact <- mean(exact_v)
    local_fpr <- mean(fpr_v)
    split_rate <- mean(n_est_per_locus > 1L)
  }

  sec_power <- NA_real_
  if ("signal_order" %in% names(t_sig) && any(t_sig$signal_order >= 2)) {
    sec_ids <- t_sig$signal_id[t_sig$signal_order >= 2]
    sec_power <- mean(sec_ids %in% matches$truth_signal_id)
  }

  missed_rate <- if (n_t == 0L) NA_real_ else length(m$missed) / n_t
  extra_rate <- length(m$extra) / max(1, n_e)

  # ---- conditional (given locus detection) resolution metrics -------------------
  attempt_rate <- if (n_loci_t == 0L) {
    NA_real_
  } else {
    as.numeric(length(detected_loci) > 0L)
  }
  # signal_count_exact_recovery is already conditioned on detection
  cond_count_exact <- count_exact
  cond_sec_power <- NA_real_
  cond_missed <- NA_real_
  if (length(detected_loci) > 0L) {
    in_det <- t_sig$locus_id %in% detected_loci
    det_ids <- t_sig$signal_id[in_det]
    if (length(det_ids) > 0L) {
      cond_missed <- mean(!det_ids %in% matches$truth_signal_id)
    }
    if ("signal_order" %in% names(t_sig)) {
      sec_det <- t_sig$signal_id[in_det & t_sig$signal_order >= 2]
      if (length(sec_det) > 0L) {
        cond_sec_power <- mean(sec_det %in% matches$truth_signal_id)
      }
    }
  }
  coverage <- if (n_t == 0L) {
    NA_real_
  } else {
    length(unique(matches$truth_signal_id)) / n_t
  }
  mean_r2 <- if (is.null(G) || np == 0L) {
    NA_real_
  } else {
    .mean_or_na(matches$representative_to_causal_r2)
  }
  mean_dist <- if (np == 0L) {
    NA_real_
  } else {
    .mean_or_na(matches$representative_to_causal_distance)
  }

  # ---- omnibus -----------------------------------------------------------------
  causal <- truth$causal_markers
  if (is.null(causal)) causal <- character()
  om <- estimates$omnibus_summary
  type1 <- NA_real_
  power <- NA_real_
  if (!is.null(om) && nrow(om) > 0L) {
    tested <- !is.na(om$p_value)
    null_rows <- tested & !(om$marker_id %in% causal)
    if (any(null_rows)) {
      type1 <- mean(om$p_value[null_rows] <= alpha)
    }
    if (length(causal) > 0L) {
      p_causal <- om$p_value[om$marker_id %in% causal & tested]
      if (length(p_causal) > 0L) power <- as.numeric(any(p_causal <= alpha))
    }
  }

  # ---- per matched-pair metrics --------------------------------------------------
  v0 <- rep(NA_real_, np)
  c_tpr <- v0; c_fdp <- v0; c_prec <- v0; c_exact <- v0; c_jac <- v0
  b_bias <- v0; b_rmse <- v0; b_cov <- v0; b_sign <- v0
  e_bias <- v0; e_rmse <- v0; e_sign <- v0
  r_bias <- v0; r_rmse <- v0; r_mae <- v0; r_acc <- v0
  rep_card <- v0; rep_exact <- v0; rep_tie <- v0
  irr_exact <- v0; irr_miss <- v0; irr_split <- v0
  breadth <- v0; dir_ok <- v0
  overall_v <- v0; cond_v <- v0

  for (i in seq_len(np)) {
    ts <- matches$truth_signal_id[[i]]
    es <- matches$est_signal_id[[i]]
    t_loc <- matches$truth_locus_id[[i]]
    e_loc <- matches$est_locus_id[[i]]

    # candidate set
    t_cand <- truth$candidate_traits[[ts]]
    if (is.null(t_cand)) t_cand <- character()
    e_cand <- estimates$candidate_sets[[es]]
    if (is.null(e_cand)) e_cand <- character()
    cm <- .set_family_metrics(t_cand, e_cand)
    c_tpr[[i]] <- cm$tpr
    c_fdp[[i]] <- cm$fdp
    c_prec[[i]] <- 1 - cm$fdp
    c_exact[[i]] <- cm$exact
    c_jac[[i]] <- cm$jaccard

    # beta effects (merged by trait key)
    tb <- .beta_for_signal(truth$beta, t_loc, ts, NULL, NULL)
    eb <- .beta_for_signal(
      estimates$beta, e_loc, es,
      m$est_signals$signal_order[m$est_signals$signal_id == es][1L],
      m$est_signals$representative_snp[m$est_signals$signal_id == es][1L]
    )
    if (!is.null(tb) && nrow(tb) > 0L &&
        !is.null(eb) && nrow(eb) > 0L) {
      idx <- match(tb$trait, eb$trait)
      has <- !is.na(idx)
      if (any(has)) {
        est_b <- eb$beta[idx[has]]
        tru_b <- tb$beta[has]
        d <- est_b - tru_b
        ok <- is.finite(d)
        if (any(ok)) {
          b_bias[[i]] <- mean(d[ok])
          b_rmse[[i]] <- sqrt(mean(d[ok]^2))
          nz <- ok & tru_b != 0
          if (any(nz)) {
            b_sign[[i]] <- mean(sign(est_b[nz]) == sign(tru_b[nz]))
          }
          if ("se" %in% names(eb)) {
            se <- eb$se[idx[has]]
            cov_ok <- ok & is.finite(se) & se > 0
            if (any(cov_ok)) {
              b_cov[[i]] <- mean(abs(d[cov_ok]) <= 1.96 * se[cov_ok])
            }
          }
        }
      }
    }

    # conditional effect eta and representation loss rho
    t_sub <- .rows_for_signal(truth$subset_table, ts, t_loc)
    e_sub <- .rows_for_signal(estimates$subset_table, es, e_loc)
    sm <- .match_subset_maps(t_sub, e_sub)
    if (sm$n > 0L) {
      de <- sm$est_eta - sm$truth_eta
      if (length(de) > 0L) {
        e_bias[[i]] <- mean(de)
        e_rmse[[i]] <- sqrt(mean(de^2))
        nz <- sm$truth_eta != 0
        if (any(nz)) {
          e_sign[[i]] <- mean(sign(sm$est_eta[nz]) ==
                                sign(sm$truth_eta[nz]))
        }
      }
      dr <- sm$est_loss - sm$truth_loss
      if (length(dr) > 0L) {
        r_bias[[i]] <- mean(dr)
        r_rmse[[i]] <- sqrt(mean(dr^2))
        r_mae[[i]] <- mean(abs(dr))
        feas_est <- sm$est_loss <= tol + 1e-10
        feas_true <- sm$truth_loss <= tol + 1e-10
        r_acc[[i]] <- mean(feas_est == feas_true)
      }
    }

    # representative set family (primary tolerance)
    t_rep_tab <- .rows_for_signal(truth$minimum_representative_sets,
                                  ts, t_loc)
    e_rep_tab <- .rows_for_signal(estimates$minimum_representative_sets,
                                  es, e_loc)
    t_reps <- .family_keys_at(t_rep_tab, tol)
    e_reps <- .family_keys_at(e_rep_tab, tol)
    rep_fam <- .set_family_metrics(t_reps, e_reps)
    rep_exact[[i]] <- rep_fam$exact
    t_min <- .min_size_at(t_rep_tab, tol)
    e_min <- .min_size_at(e_rep_tab, tol)
    rep_card[[i]] <- as.numeric(t_min == e_min)
    rep_tie[[i]] <- if (length(t_reps) == 0L) {
      if (length(e_reps) == 0L) 1 else 0
    } else {
      length(intersect(t_reps, e_reps)) / length(t_reps)
    }

    # irreducible module family
    t_mod_tab <- .rows_for_signal(truth$irreducible_modules, ts, t_loc)
    e_mod_tab <- .rows_for_signal(estimates$irreducible_modules, es, e_loc)
    t_mods <- .family_keys_at(t_mod_tab, tol)
    e_mods <- .family_keys_at(e_mod_tab, tol)
    mod_fam <- .set_family_metrics(t_mods, e_mods)
    irr_exact[[i]] <- mod_fam$exact
    t_sets <- .family_sets_at(t_mod_tab, tol)
    e_sets <- .family_sets_at(e_mod_tab, tol)
    joint <- lengths(t_sets) > 1L
    if (any(joint)) {
      joint_keys <- names(t_sets)[joint]
      missed <- 0L
      fragmented <- 0L
      for (jk in joint_keys) {
        T_set <- t_sets[[jk]]
        if (!(jk %in% names(e_sets))) {
          missed <- missed + 1L
        }
        n_sub <- sum(vapply(e_sets, function(E) {
          length(E) < length(T_set) && all(E %in% T_set)
        }, logical(1)))
        if (n_sub >= 2L) fragmented <- fragmented + 1L
      }
      irr_miss[[i]] <- missed / length(joint_keys)
      irr_split[[i]] <- fragmented / length(joint_keys)
    }

    # effect breadth and direction
    t_breadth <- if (!is.null(truth$effect_breadth) &&
                     ts %in% names(truth$effect_breadth)) {
      as.numeric(truth$effect_breadth[[ts]])
    } else {
      length(t_cand)
    }
    breadth[[i]] <- length(e_cand) - t_breadth
    e_dir <- .direction_from_betas(
      if (!is.null(eb) && nrow(eb) > 0L) {
        eb$beta[eb$trait %in% e_cand]
      } else {
        numeric()
      },
      length(e_cand)
    )
    t_dir <- if (!is.null(truth$effect_direction) &&
                 ts %in% names(truth$effect_direction)) {
      truth$effect_direction[[ts]]
    } else {
      NA_character_
    }
    if (!is.na(t_dir)) dir_ok[[i]] <- as.numeric(e_dir == t_dir)

    overall_v[[i]] <- as.numeric(isTRUE(cm$exact == 1) &&
                                   isTRUE(rep_fam$exact == 1) &&
                                   isTRUE(mod_fam$exact == 1))
    if (isTRUE(cm$exact == 1)) {
      cond_v[[i]] <- as.numeric(isTRUE(rep_fam$exact == 1) &&
                                  isTRUE(mod_fam$exact == 1))
    }
  }

  # ---- ASSET / ablation metrics --------------------------------------------------
  a_power <- a_jac <- a_prec <- a_rec <- a_dir <- NA_real_
  a_lead_infl <- a_res_infl <- a_pseudo <- NA_real_
  asset <- estimates$asset
  if (!is.null(asset) && n_t > 0L) {
    sids <- t_sig$signal_id
    # estimate-side ids via the genomic/LD match (never name equality)
    est_sid <- matches$est_signal_id[match(sids, matches$truth_signal_id)]
    est_lid <- matches$est_locus_id[
      match(t_sig$locus_id, matches$truth_locus_id)]
    get0 <- function(lst, id) {
      if (is.null(lst) || is.na(id) || !id %in% names(lst)) {
        NULL
      } else {
        lst[[id]]
      }
    }
    t_breadths <- vapply(sids, function(s) {
      if (!is.null(truth$effect_breadth) &&
          s %in% names(truth$effect_breadth)) {
        as.numeric(truth$effect_breadth[[s]])
      } else {
        tc <- truth$candidate_traits[[s]]
        if (is.null(tc)) 0 else length(tc)
      }
    }, numeric(1))
    widths_of <- function(entry, key = c("signal", "locus")) {
      key <- match.arg(key)
      ids <- if (key == "signal") est_sid else est_lid
      vapply(seq_along(sids), function(k) {
        st <- if (!is.null(entry) && !is.null(entry$sets)) {
          get0(entry$sets, ids[[k]])
        }
        if (is.null(st)) 0 else length(st)
      }, numeric(1))
    }
    cp <- asset$condped
    lead <- asset$lead
    resolved <- asset$resolved
    if (!is.null(cp)) {
      if (!is.null(cp$p)) {
        pv <- vapply(est_sid, function(es) {
          p <- get0(cp$p, es)
          if (is.null(p) || length(p) == 0L || !is.finite(p)) {
            NA_real_
          } else {
            as.numeric(p <= alpha)
          }
        }, numeric(1))
        a_power <- .mean_or_na(pv)
      }
      jj <- pp <- rr <- rep(NA_real_, n_t)
      for (k in seq_len(n_t)) {
        s <- sids[[k]]
        est_set <- if (!is.null(cp$sets)) get0(cp$sets, est_sid[[k]])
        if (is.null(est_set)) est_set <- character()
        t_set <- truth$candidate_traits[[s]]
        if (is.null(t_set)) t_set <- character()
        fm <- .set_family_metrics(t_set, est_set)
        jj[[k]] <- fm$jaccard
        pp[[k]] <- 1 - fm$fdp
        rr[[k]] <- fm$tpr
      }
      a_jac <- .mean_or_na(jj)
      a_prec <- .mean_or_na(pp)
      a_rec <- .mean_or_na(rr)
      if (!is.null(cp$direction) && !is.null(truth$effect_direction)) {
        dv <- vapply(seq_len(n_t), function(k) {
          td <- truth$effect_direction[[sids[[k]]]]
          ed <- get0(cp$direction, est_sid[[k]])
          if (is.null(td) || is.null(ed)) {
            NA_real_
          } else {
            as.numeric(unname(ed) == unname(td))
          }
        }, numeric(1))
        a_dir <- .mean_or_na(dv)
      }
    }
    # fall back to ASSET-side direction representations (derived from
    # the real positive/negative subsets) when the CondPED entry has
    # none; resolved signals first, then the marginal lead
    if (is.na(a_dir) && !is.null(truth$effect_direction)) {
      for (entry in list(list(e = resolved, key = "signal"),
                         list(e = lead, key = "locus"))) {
        ent <- entry$e
        if (is.null(ent) || is.null(ent$direction)) next
        ids <- if (entry$key == "signal") est_sid else est_lid
        dv <- vapply(seq_len(n_t), function(k) {
          td <- truth$effect_direction[[sids[[k]]]]
          ed <- get0(ent$direction, ids[[k]])
          if (is.null(td) || is.null(ed)) {
            NA_real_
          } else {
            as.numeric(unname(ed) == unname(td))
          }
        }, numeric(1))
        a_dir <- .mean_or_na(dv)
        if (!is.na(a_dir)) break
      }
    }
    if (!is.null(lead)) {
      a_lead_infl <- mean(widths_of(lead, "locus") - t_breadths)
    }
    if (!is.null(resolved)) {
      a_res_infl <- mean(widths_of(resolved, "signal") - t_breadths)
    }
    if (!is.null(lead) && !is.null(cp)) {
      a_pseudo <- mean(widths_of(lead, "locus") >= 2 &
                         widths_of(cp, "signal") <= 1)
    }
  }

  row <- data.frame(
    locus_detection_rate = locus_detection_rate,
    signal_count_exact_recovery = count_exact,
    secondary_signal_power = sec_power,
    missed_signal_rate = missed_rate,
    extra_signal_rate = extra_rate,
    signal_coverage = coverage,
    representative_to_causal_r2 = mean_r2,
    representative_to_causal_distance = mean_dist,
    residual_local_false_positive_rate = local_fpr,
    signal_resolution_attempt_rate = attempt_rate,
    conditional_signal_count_exact_recovery = cond_count_exact,
    conditional_secondary_signal_power = cond_sec_power,
    conditional_missed_signal_rate = cond_missed,
    truth_locus_split_rate = split_rate,
    type1_omnibus = type1,
    power_omnibus = power,
    beta_bias = .mean_or_na(b_bias),
    beta_rmse = .mean_or_na(b_rmse),
    beta_coverage = .mean_or_na(b_cov),
    beta_sign_accuracy = .mean_or_na(b_sign),
    candidate_tpr = .mean_or_na(c_tpr),
    candidate_fdp = .mean_or_na(c_fdp),
    candidate_precision = .mean_or_na(c_prec),
    candidate_exact_recovery = .mean_or_na(c_exact),
    candidate_jaccard = .mean_or_na(c_jac),
    effect_breadth_bias = .mean_or_na(breadth),
    direction_recovery = .mean_or_na(dir_ok),
    conditional_effect_bias = .mean_or_na(e_bias),
    conditional_effect_rmse = .mean_or_na(e_rmse),
    conditional_effect_sign_accuracy = .mean_or_na(e_sign),
    representation_loss_bias = .mean_or_na(r_bias),
    representation_loss_rmse = .mean_or_na(r_rmse),
    representation_loss_mae = .mean_or_na(r_mae),
    representation_threshold_accuracy = .mean_or_na(r_acc),
    representation_map_mae = .mean_or_na(r_mae),
    rep_min_cardinality_recovery = .mean_or_na(rep_card),
    rep_exact_family_recovery = .mean_or_na(rep_exact),
    rep_tie_coverage = .mean_or_na(rep_tie),
    irr_family_recovery = .mean_or_na(irr_exact),
    irr_missed_rate = .mean_or_na(irr_miss),
    irr_extra_split_rate = .mean_or_na(irr_split),
    asset_power = a_power,
    asset_trait_set_jaccard = a_jac,
    asset_trait_set_precision = a_prec,
    asset_trait_set_recall = a_rec,
    asset_direction_recovery = a_dir,
    lead_trait_breadth_inflation = a_lead_infl,
    resolved_trait_breadth_inflation = a_res_infl,
    pseudo_multitrait_interpretation_rate = a_pseudo,
    overall_recovery = .mean_or_na(overall_v),
    conditional_recovery = .mean_or_na(cond_v),
    runtime = if (is.null(runtime)) NA_real_ else runtime,
    stringsAsFactors = FALSE
  )
  # scenario descriptors used for grouping
  for (g in c("experiment", "n", "architecture", "locus_pve",
              "correlation", "target_loss", "tolerance")) {
    val <- settings[[g]]
    if (is.null(val) || length(val) != 1L) {
      val <- if (g == "tolerance") tol else NA
    }
    row[[g]] <- if (is.numeric(val)) as.numeric(val) else as.character(val)
  }
  row
}

#' Match truth and estimated signals by genomic position
#'
#' Loci are matched by chromosome plus interval overlap; estimated
#' locus intervals are taken from explicit `start`/`end` columns when
#' present, parsed from `chr:start-end` locus ids otherwise (falling
#' back to the range of representative positions). Signals within
#' overlapping loci are matched one-to-one and greedily: candidate
#' pairs are sorted by r^2 (descending, when `G` is supplied), then
#' representative-to-causal distance (ascending), then marker id, and
#' each signal is used at most once. Matching never uses signal-id
#' string equality.
#'
#' @param truth Truth component with `signals` (and usually `loci`).
#' @param estimates Estimates list with `signals` and optionally
#'   `positions` (named marker -> position vector).
#' @param G Optional n x p genotype matrix with marker-id column
#'   names; when supplied, matching prefers the highest squared
#'   genotype correlation (via `.r2_between`).
#' @return A list with `matches` (data.frame with `truth_signal_id`,
#'   `est_signal_id`, `truth_locus_id`, `est_locus_id`,
#'   `representative_to_causal_r2`,
#'   `representative_to_causal_distance`), `missed` (truth signal
#'   ids), `extra` (estimated signal ids), `locus_matches`
#'   (data.frame with `truth_locus_id`, `est_locus_id`),
#'   `matched_locus_count`, and the normalised `truth_signals`,
#'   `est_signals`, `truth_loci` and `est_loci` tables.
#' @keywords internal
.match_signals <- function(truth, estimates, G = NULL) {
  positions <- estimates$positions
  pos_of <- function(ids) {
    if (is.null(positions) || length(ids) == 0L) {
      return(rep(NA_real_, length(ids)))
    }
    unname(suppressWarnings(as.numeric(positions[ids])))
  }

  # ---- normalised truth signal table -------------------------------------------
  t_sig <- truth$signals
  t_tab <- .empty_signal_table()
  if (!is.null(t_sig) && nrow(t_sig) > 0L) {
    t_chr <- if ("chromosome" %in% names(t_sig)) {
      as.character(t_sig$chromosome)
    } else if (!is.null(truth$loci) &&
               all(c("locus_id", "chromosome") %in% names(truth$loci))) {
      as.character(truth$loci$chromosome[
        match(t_sig$locus_id, truth$loci$locus_id)
      ])
    } else {
      rep(NA_character_, nrow(t_sig))
    }
    t_pos <- if ("position" %in% names(t_sig)) {
      as.numeric(t_sig$position)
    } else {
      rep(NA_real_, nrow(t_sig))
    }
    miss <- is.na(t_pos)
    if (any(miss)) t_pos[miss] <- pos_of(t_sig$representative_snp[miss])
    t_tab <- data.frame(
      signal_id = as.character(t_sig$signal_id),
      locus_id = as.character(t_sig$locus_id),
      representative_snp = as.character(t_sig$representative_snp),
      chromosome = t_chr,
      position = t_pos,
      signal_order = if ("signal_order" %in% names(t_sig)) {
        as.numeric(t_sig$signal_order)
      } else {
        rep(NA_real_, nrow(t_sig))
      },
      stringsAsFactors = FALSE
    )
  }

  # ---- normalised estimated signal table -----------------------------------------
  e_sig <- estimates$signals
  e_tab <- .empty_signal_table()
  if (!is.null(e_sig) && nrow(e_sig) > 0L) {
    e_ids <- if ("signal_id" %in% names(e_sig)) {
      as.character(e_sig$signal_id)
    } else if ("signal_order" %in% names(e_sig)) {
      paste(as.character(e_sig$locus_id), e_sig$signal_order, sep = "::")
    } else {
      as.character(e_sig$representative_snp)
    }
    e_chr <- if ("chromosome" %in% names(e_sig)) {
      as.character(e_sig$chromosome)
    } else {
      rep(NA_character_, nrow(e_sig))
    }
    e_pos <- if ("position" %in% names(e_sig)) {
      as.numeric(e_sig$position)
    } else {
      rep(NA_real_, nrow(e_sig))
    }
    miss <- is.na(e_pos)
    if (any(miss)) e_pos[miss] <- pos_of(e_sig$representative_snp[miss])
    e_tab <- data.frame(
      signal_id = e_ids,
      locus_id = as.character(e_sig$locus_id),
      representative_snp = as.character(e_sig$representative_snp),
      chromosome = e_chr,
      position = e_pos,
      signal_order = if ("signal_order" %in% names(e_sig)) {
        as.numeric(e_sig$signal_order)
      } else {
        rep(NA_real_, nrow(e_sig))
      },
      stringsAsFactors = FALSE
    )
  }

  # ---- locus tables --------------------------------------------------------------
  t_loc <- truth$loci
  if (!is.null(t_loc) && nrow(t_loc) > 0L) {
    t_loc_tab <- data.frame(
      locus_id = as.character(t_loc$locus_id),
      chromosome = as.character(t_loc$chromosome),
      start = as.numeric(t_loc$start),
      end = as.numeric(t_loc$end),
      stringsAsFactors = FALSE
    )
  } else {
    t_loc_tab <- .loci_from_signals(t_tab)
  }
  has_bounds <- !is.null(e_sig) &&
    all(c("start", "end") %in% names(e_sig))
  if (has_bounds && nrow(e_sig) > 0L) {
    e_tab$start <- as.numeric(e_sig$start)
    e_tab$end <- as.numeric(e_sig$end)
  }
  e_loc_tab <- .est_locus_table(e_tab, has_bounds)

  # ---- locus overlap --------------------------------------------------------------
  locus_matches <- data.frame(
    truth_locus_id = character(), est_locus_id = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(t_loc_tab) > 0L && nrow(e_loc_tab) > 0L) {
    lm_rows <- list()
    for (ti in seq_len(nrow(t_loc_tab))) {
      for (ei in seq_len(nrow(e_loc_tab))) {
        tc <- t_loc_tab$chromosome[[ti]]
        ec <- e_loc_tab$chromosome[[ei]]
        if (is.na(tc) || is.na(ec) || tc != ec) next
        if (is.na(t_loc_tab$start[[ti]]) || is.na(e_loc_tab$start[[ei]])) next
        if (e_loc_tab$start[[ei]] <= t_loc_tab$end[[ti]] &&
            e_loc_tab$end[[ei]] >= t_loc_tab$start[[ti]]) {
          lm_rows[[length(lm_rows) + 1L]] <- data.frame(
            truth_locus_id = t_loc_tab$locus_id[[ti]],
            est_locus_id = e_loc_tab$locus_id[[ei]],
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(lm_rows) > 0L) locus_matches <- do.call(rbind, lm_rows)
  }

  # ---- candidate signal pairs -------------------------------------------------------
  matches <- data.frame(
    truth_signal_id = character(), est_signal_id = character(),
    truth_locus_id = character(), est_locus_id = character(),
    representative_to_causal_r2 = numeric(),
    representative_to_causal_distance = numeric(),
    stringsAsFactors = FALSE
  )
  if (nrow(locus_matches) > 0L) {
    cand <- list()
    for (li in seq_len(nrow(locus_matches))) {
      t_idx <- which(t_tab$locus_id == locus_matches$truth_locus_id[[li]])
      e_idx <- which(e_tab$locus_id == locus_matches$est_locus_id[[li]])
      for (tj in t_idx) {
        for (ej in e_idx) {
          r2 <- if (!is.null(G)) {
            .r2_between(G, colnames(G), e_tab$representative_snp[[ej]],
                        t_tab$representative_snp[[tj]])
          } else {
            NA_real_
          }
          cand[[length(cand) + 1L]] <- data.frame(
            truth_signal_id = t_tab$signal_id[[tj]],
            est_signal_id = e_tab$signal_id[[ej]],
            truth_locus_id = t_tab$locus_id[[tj]],
            est_locus_id = e_tab$locus_id[[ej]],
            r2 = r2,
            dist = abs(e_tab$position[[ej]] - t_tab$position[[tj]]),
            est_rep = e_tab$representative_snp[[ej]],
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(cand) > 0L) {
      cand <- do.call(rbind, cand)
      r2_key <- ifelse(is.na(cand$r2), -1, cand$r2)
      d_key <- ifelse(is.na(cand$dist), Inf, cand$dist)
      ord <- order(-r2_key, d_key, cand$est_rep, cand$truth_signal_id)
      t_used <- e_used <- character()
      keep <- logical(nrow(cand))
      for (k in ord) {
        if (cand$truth_signal_id[[k]] %in% t_used ||
            cand$est_signal_id[[k]] %in% e_used) next
        keep[[k]] <- TRUE
        t_used <- c(t_used, cand$truth_signal_id[[k]])
        e_used <- c(e_used, cand$est_signal_id[[k]])
      }
      matches <- data.frame(
        truth_signal_id = cand$truth_signal_id[keep],
        est_signal_id = cand$est_signal_id[keep],
        truth_locus_id = cand$truth_locus_id[keep],
        est_locus_id = cand$est_locus_id[keep],
        representative_to_causal_r2 = cand$r2[keep],
        representative_to_causal_distance = cand$dist[keep],
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    matches = matches,
    missed = setdiff(t_tab$signal_id, matches$truth_signal_id),
    extra = setdiff(e_tab$signal_id, matches$est_signal_id),
    locus_matches = locus_matches,
    matched_locus_count = length(unique(locus_matches$truth_locus_id)),
    truth_signals = t_tab,
    est_signals = e_tab,
    truth_loci = t_loc_tab,
    est_loci = e_loc_tab
  )
}

#' Empty normalised signal table
#' @keywords internal
.empty_signal_table <- function() {
  data.frame(
    signal_id = character(), locus_id = character(),
    representative_snp = character(), chromosome = character(),
    position = numeric(), signal_order = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Derive a locus table from a normalised signal table
#' @param sig_tab Normalised signal table.
#' @keywords internal
.loci_from_signals <- function(sig_tab) {
  if (nrow(sig_tab) == 0L) {
    return(data.frame(
      locus_id = character(), chromosome = character(),
      start = numeric(), end = numeric(), stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(unique(sig_tab$locus_id), function(lid) {
    sub <- sig_tab[sig_tab$locus_id == lid, ]
    chr <- sub$chromosome[!is.na(sub$chromosome)]
    pos <- sub$position[is.finite(sub$position)]
    data.frame(
      locus_id = lid,
      chromosome = if (length(chr) > 0L) chr[[1L]] else NA_character_,
      start = if (length(pos) > 0L) min(pos) else NA_real_,
      end = if (length(pos) > 0L) max(pos) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Build the estimated locus table (explicit bounds, parsed ids, or
#' representative-position ranges)
#' @param e_tab Normalised estimated signal table.
#' @param has_bounds Whether explicit `start`/`end` columns exist.
#' @keywords internal
.est_locus_table <- function(e_tab, has_bounds) {
  if (nrow(e_tab) == 0L) {
    return(data.frame(
      locus_id = character(), chromosome = character(),
      start = numeric(), end = numeric(), stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(unique(e_tab$locus_id), function(lid) {
    sub <- e_tab[e_tab$locus_id == lid, ]
    chr <- sub$chromosome[!is.na(sub$chromosome)]
    chr <- if (length(chr) > 0L) chr[[1L]] else NA_character_
    start <- end <- NA_real_
    if (has_bounds && any(is.finite(sub$start)) &&
        any(is.finite(sub$end))) {
      start <- min(sub$start, na.rm = TRUE)
      end <- max(sub$end, na.rm = TRUE)
    } else {
      parsed <- .parse_locus_id(lid)
      if (!is.null(parsed)) {
        if (is.na(chr)) chr <- parsed$chromosome
        start <- parsed$start
        end <- parsed$end
      }
    }
    if (is.na(start) || is.na(end)) {
      pos <- sub$position[is.finite(sub$position)]
      if (length(pos) > 0L) {
        start <- min(pos)
        end <- max(pos)
      }
    }
    data.frame(
      locus_id = lid, chromosome = chr, start = start, end = end,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Parse a `chr:start-end` locus id
#' @param id Locus id string.
#' @return List with `chromosome`, `start`, `end`, or `NULL` when the
#'   id does not parse.
#' @keywords internal
.parse_locus_id <- function(id) {
  parts <- regmatches(
    id, regexec("^(.*):([0-9][0-9,]*)-([0-9][0-9,]*)$", id)
  )[[1L]]
  if (length(parts) != 4L) return(NULL)
  list(
    chromosome = parts[[2L]],
    start = as.numeric(gsub(",", "", parts[[3L]])),
    end = as.numeric(gsub(",", "", parts[[4L]]))
  )
}

#' Stop on duplicated keys in a table
#'
#' @param tab Data.frame (or `NULL`).
#' @param keys Key columns; intersected with the columns present.
#' @param label Table label for the error message.
#' @keywords internal
.check_key_uniqueness <- function(tab, keys, label) {
  if (is.null(tab) || !is.data.frame(tab) || nrow(tab) == 0L) {
    return(invisible(TRUE))
  }
  keys <- intersect(keys, names(tab))
  if (length(keys) == 0L) return(invisible(TRUE))
  if (anyDuplicated(tab[keys])) {
    .stop_invalid_input(
      "Duplicate keys in %s (key: %s).",
      label, paste(keys, collapse = " + ")
    )
  }
  invisible(TRUE)
}

#' Subset a table to one signal (and locus, when columns exist)
#' @param tab Data.frame or `NULL`.
#' @param signal_id Signal id.
#' @param locus_id Optional locus id.
#' @keywords internal
.rows_for_signal <- function(tab, signal_id, locus_id = NULL) {
  if (is.null(tab) || !is.data.frame(tab) || nrow(tab) == 0L) {
    return(tab)
  }
  keep <- rep(TRUE, nrow(tab))
  if ("locus_id" %in% names(tab) && !is.null(locus_id)) {
    keep <- keep & tab$locus_id == locus_id
  }
  if ("signal_id" %in% names(tab)) {
    keep <- keep & tab$signal_id == signal_id
  }
  tab[keep, , drop = FALSE]
}

#' Extract the beta rows of one signal from a long beta table
#'
#' Rows are selected by `locus_id` plus `signal_id` when present,
#' else by `signal_order`, else by `representative_snp`.
#'
#' @param beta_tab Long beta data.frame or `NULL`.
#' @param locus_id,signal_id,signal_order,representative_snp Signal
#'   identifiers (any may be `NULL`).
#' @keywords internal
.beta_for_signal <- function(beta_tab, locus_id, signal_id,
                             signal_order, representative_snp) {
  if (is.null(beta_tab) || !is.data.frame(beta_tab) ||
      nrow(beta_tab) == 0L) {
    return(NULL)
  }
  keep <- rep(TRUE, nrow(beta_tab))
  if ("locus_id" %in% names(beta_tab) && !is.null(locus_id)) {
    keep <- keep & beta_tab$locus_id == locus_id
  }
  if ("signal_id" %in% names(beta_tab) && !is.null(signal_id)) {
    keep <- keep & beta_tab$signal_id == signal_id
  } else if ("signal_order" %in% names(beta_tab) &&
             !is.null(signal_order) && !is.na(signal_order)) {
    keep <- keep & beta_tab$signal_order == signal_order
  } else if ("representative_snp" %in% names(beta_tab) &&
             !is.null(representative_snp)) {
    keep <- keep & beta_tab$representative_snp == representative_snp
  }
  beta_tab[keep, , drop = FALSE]
}

#' Direction pattern from candidate betas (frozen rule)
#'
#' @param beta_vals Beta values over the candidate traits.
#' @param n_cand Number of candidate traits.
#' @return One of `not_applicable`, `single_trait`, `concordant`,
#'   `antagonistic`, `mixed`.
#' @keywords internal
.direction_from_betas <- function(beta_vals, n_cand) {
  if (n_cand == 0L) return("not_applicable")
  if (n_cand == 1L) return("single_trait")
  s <- sign(beta_vals[is.finite(beta_vals)])
  s <- s[s != 0]
  if (length(s) == 0L || all(s == s[[1L]])) return("concordant")
  if (n_cand == 2L) "antagonistic" else "mixed"
}

#' Mean of the finite part of a vector (`NA` when empty)
#' @param x Numeric vector.
#' @keywords internal
.mean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

#' Set/family metrics under the frozen empty-set rules
#'
#' @param truth_keys,est_keys Character vectors of normalised keys
#'   (traits or trait-set keys).
#' @return List with `tpr`, `fdp`, `exact`, `jaccard`.
#' @keywords internal
.set_family_metrics <- function(truth_keys, est_keys) {
  truth_keys <- unique(truth_keys)
  est_keys <- unique(est_keys)
  n_t <- length(truth_keys)
  n_e <- length(est_keys)
  tp <- length(intersect(truth_keys, est_keys))
  if (n_t == 0L && n_e == 0L) {
    return(list(tpr = NA_real_, fdp = 0, exact = 1, jaccard = 1))
  }
  if (n_e == 0L) {
    return(list(tpr = 0, fdp = 0, exact = 0, jaccard = 0))
  }
  if (n_t == 0L) {
    return(list(tpr = NA_real_, fdp = 1, exact = 0, jaccard = 0))
  }
  list(
    tpr = tp / n_t,
    fdp = (n_e - tp) / max(1, n_e),
    exact = as.numeric(n_t == n_e && tp == n_t),
    jaccard = tp / (n_t + n_e - tp)
  )
}

#' Family trait keys at one tolerance
#' @param tab Extraction table (or `NULL`).
#' @param tol Tolerance.
#' @keywords internal
.family_keys_at <- function(tab, tol) {
  if (is.null(tab) || nrow(tab) == 0L) return(character())
  sort(unique(tab$trait_key[tab$tolerance == tol]))
}

#' Family sets (list of trait vectors, named by key) at one tolerance
#' @keywords internal
.family_sets_at <- function(tab, tol) {
  if (is.null(tab) || nrow(tab) == 0L) return(list())
  keep <- tab$tolerance == tol
  sets <- tab$trait_set[keep]
  names(sets) <- tab$trait_key[keep]
  sets
}

#' Minimum set size at one tolerance (Inf when empty)
#' @keywords internal
.min_size_at <- function(tab, tol) {
  if (is.null(tab) || nrow(tab) == 0L) return(Inf)
  keep <- tab$tolerance == tol
  if (!any(keep)) return(Inf)
  min(tab$set_size[keep])
}

#' Grouped means and Monte Carlo standard errors
#'
#' @param df Replicate metrics data.frame.
#' @param metrics Metric columns to summarise.
#' @param group_cols Grouping columns (possibly empty).
#' @param conf_level Confidence level.
#' @return List with `summary` and `mcse` data.frames.
#' @keywords internal
.summarise_metrics <- function(df, metrics, group_cols, conf_level) {
  metrics <- intersect(metrics, names(df))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  if (length(group_cols) == 0L) {
    groups <- list(seq_len(nrow(df)))
    group_keys <- NULL
  } else {
    key <- interaction(
      lapply(df[group_cols], function(x) {
        x <- as.character(x)
        x[is.na(x)] <- "<NA>"
        x
      }),
      drop = TRUE, lex.order = TRUE
    )
    groups <- split(seq_len(nrow(df)), key)
    group_keys <- lapply(groups, function(idx) {
      df[idx[1L], group_cols, drop = FALSE]
    })
  }
  sum_rows <- list()
  mcse_rows <- list()
  for (gi in seq_along(groups)) {
    idx <- groups[[gi]]
    srow <- if (!is.null(group_keys)) {
      as.list(group_keys[[gi]])
    } else {
      list()
    }
    srow$n <- length(idx)
    mrow <- srow
    for (met in metrics) {
      x <- df[[met]][idx]
      x <- x[!is.na(x)]
      n <- length(x)
      if (n == 0L) {
        srow[[met]] <- NA_real_
        srow[[paste0(met, "_ci_lo")]] <- NA_real_
        srow[[paste0(met, "_ci_hi")]] <- NA_real_
        mrow[[met]] <- NA_real_
        next
      }
      mn <- mean(x)
      is_prop <- all(x %in% c(0, 1))
      se <- if (n > 1L) {
        if (is_prop) sqrt(mn * (1 - mn) / n) else stats::sd(x) / sqrt(n)
      } else {
        NA_real_
      }
      srow[[met]] <- mn
      srow[[paste0(met, "_ci_lo")]] <- if (is.na(se)) NA_real_ else mn - z * se
      srow[[paste0(met, "_ci_hi")]] <- if (is.na(se)) NA_real_ else mn + z * se
      mrow[[met]] <- se
    }
    sum_rows[[gi]] <- as.data.frame(srow, stringsAsFactors = FALSE,
                                    optional = TRUE)
    mcse_rows[[gi]] <- as.data.frame(mrow, stringsAsFactors = FALSE,
                                     optional = TRUE)
  }
  list(
    summary = do.call(rbind, sum_rows),
    mcse = do.call(rbind, mcse_rows)
  )
}

#' Match truth and estimate subset maps for eta/rho metrics
#'
#' Rows are matched by `representing_key`; eta elements are aligned by
#' complement trait names. Only finite values enter the metric pools.
#'
#' @param truth_tab,est_tab Subset tables with `representing_key`,
#'   `conditional_effect` and `representation_loss`.
#' @return List with `n` (matched rows), `truth_eta`, `est_eta`,
#'   `truth_loss`, `est_loss` (aligned numeric vectors).
#' @keywords internal
.match_subset_maps <- function(truth_tab, est_tab) {
  out <- list(n = 0L, truth_eta = numeric(), est_eta = numeric(),
              truth_loss = numeric(), est_loss = numeric())
  if (is.null(truth_tab) || is.null(est_tab) ||
      nrow(truth_tab) == 0L || nrow(est_tab) == 0L) {
    return(out)
  }
  keys <- intersect(truth_tab$representing_key,
                    est_tab$representing_key)
  out$n <- length(keys)
  for (k in keys) {
    rt <- truth_tab[truth_tab$representing_key == k, ][1L, ]
    re <- est_tab[est_tab$representing_key == k, ][1L, ]
    tl <- rt$representation_loss
    el <- re$representation_loss
    if (is.finite(tl) && is.finite(el)) {
      out$truth_loss <- c(out$truth_loss, tl)
      out$est_loss <- c(out$est_loss, el)
    }
    et <- rt$conditional_effect[[1L]]
    ee <- re$conditional_effect[[1L]]
    if (!is.null(et) && !is.null(ee) && length(et) > 0L) {
      ee <- ee[names(et)]
      keep <- is.finite(et) & is.finite(ee)
      out$truth_eta <- c(out$truth_eta, unname(et[keep]))
      out$est_eta <- c(out$est_eta, unname(ee[keep]))
    }
  }
  out
}
