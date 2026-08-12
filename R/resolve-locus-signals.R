# Stepwise locus-level signal resolution core (Stage 6B-3).
#
# Decides HOW a locus's representative signals are selected
# sequentially from conditional scans. No conditional statistic is
# recomputed here: every step calls `.build_conditional_projection()`
# + `.conditional_mt_scan()` (Stage 6B-2) with the fixed global-null
# V. No final joint GLS and no signal-specific beta are produced; the
# stepwise conditional betas are never exposed downstream.

#' Select the next signal from one conditional scan
#'
#' Applies the local screening rule and picks the winner
#' deterministically. `M_remaining` is fixed BEFORE the scan (the
#' number of eligible remaining regional QC markers) and is never
#' redefined by the resulting p-values. Only markers with a usable
#' conditional statistic (finite Q and p-value) can win.
#'
#' @param tab Conditional scan table (marker_id, Q, df, p_value,
#'   status).
#' @param chromosome,position Named marker annotation vectors.
#' @param screening_rule One of `"within_locus_bonferroni"`, `"fixed"`,
#'   `"none"`.
#' @param alpha_signal Level for the bonferroni/none rules.
#' @param fixed_p_threshold Raw-p threshold for the `"fixed"` rule.
#' @param M_remaining Eligible candidate count fixed before the scan.
#'
#' @return A list with `found`, `marker_id`, `p_raw`, `p_adjusted`,
#'   `n_passing`, `M_remaining`.
#' @keywords internal
.select_next_signal <- function(tab, chromosome, position,
                                screening_rule, alpha_signal,
                                fixed_p_threshold, M_remaining) {
  usable <- is.finite(tab$Q) & is.finite(tab$p_value)
  tab <- tab[usable, , drop = FALSE]
  not_found <- list(found = FALSE, marker_id = NA_character_,
                    p_raw = NA_real_, p_adjusted = NA_real_,
                    n_passing = 0L, M_remaining = M_remaining)
  if (nrow(tab) == 0L) return(not_found)

  p_adjusted <- switch(screening_rule,
    within_locus_bonferroni = pmin(1, M_remaining * tab$p_value),
    fixed = tab$p_value,
    none = tab$p_value
  )
  pass <- switch(screening_rule,
    within_locus_bonferroni = p_adjusted <= alpha_signal,
    fixed = tab$p_value <= fixed_p_threshold,
    none = tab$p_value <= alpha_signal
  )
  n_passing <- sum(pass)
  if (n_passing == 0L) return(not_found)

  cand <- tab[pass, , drop = FALSE]
  cand_adj <- p_adjusted[pass]
  # smallest raw conditional p; deterministic tie-break by
  # chromosome -> position -> marker_id (never by input row order)
  ord <- order(cand$p_value,
               as.character(chromosome[cand$marker_id]),
               position[cand$marker_id],
               cand$marker_id)
  winner <- cand$marker_id[ord[1L]]
  list(
    found = TRUE,
    marker_id = winner,
    p_raw = cand$p_value[ord[1L]],
    p_adjusted = cand_adj[ord[1L]],
    n_passing = n_passing,
    M_remaining = M_remaining
  )
}

#' Resolve the representative signals of one locus stepwise
#'
#' Fixed algorithm: the first signal is the locus lead SNP (marginal
#' evidence). Then repeat: condition on all selected signals, rescan
#' ALL remaining regional QC markers, apply the local screening rule,
#' pick the smallest conditional raw p (deterministic tie-break), stop
#' when nobody passes or `max_signals` is reached.
#'
#' @param locus_id Locus identifier.
#' @param lead_snp Lead SNP marker id (first signal).
#' @param members Character vector of the locus's regional QC markers.
#' @param null_fit The single global null fit (with rotation).
#' @param G `n x p` genotype matrix.
#' @param marker_ids Marker order of `G`.
#' @param chromosome,position Named marker annotation vectors.
#' @param lead_marginal_p Marginal omnibus p-value of the lead.
#' @param screening_rule Local screening rule.
#' @param alpha_signal Level for bonferroni/none rules.
#' @param fixed_p_threshold Raw-p threshold for the `"fixed"` rule.
#' @param max_signals Maximum number of signals per locus.
#' @param rank_tol Tolerance forwarded to [`.safe_inverse()`].
#' @param keep_history Logical; keep each step's conditional scan.
#'
#' @return A list with `signals` (data.frame: `locus_id`,
#'   `signal_order`, `representative_snp`,
#'   `conditioning_before_selection` (list-column),
#'   `conditioning_key`, `n_conditioning`, `conditional_p_raw`,
#'   `conditional_p_adjusted`, `screening_rule`, `M_remaining`,
#'   `selection_status`), `selected` (ordered signal marker ids),
#'   `history` (per-step conditional scan tables, or `NULL`),
#'   `status` and `diagnostics` (with `stop_reason`, `n_signals`,
#'   `n_steps`).
#' @keywords internal
.resolve_one_locus_signals <- function(locus_id, lead_snp, members,
                                       null_fit, G,
                                       marker_ids = colnames(G),
                                       chromosome, position,
                                       lead_marginal_p = NA_real_,
                                       screening_rule = c(
                                         "within_locus_bonferroni",
                                         "fixed", "none"),
                                       alpha_signal = 0.05,
                                       fixed_p_threshold = 0.05,
                                       max_signals = 5L,
                                       rank_tol = sqrt(.Machine$double.eps),
                                       keep_history = TRUE) {
  screening_rule <- match.arg(screening_rule)
  .check_prob(alpha_signal, "alpha_signal")
  .check_prob(fixed_p_threshold, "fixed_p_threshold")
  if (!is.numeric(max_signals) || length(max_signals) != 1L ||
      !is.finite(max_signals) || max_signals < 1L) {
    .stop_invalid_input("max_signals must be a positive integer.")
  }
  members <- unique(as.character(members))
  if (!lead_snp %in% members) {
    .stop_invalid_input("lead_snp must be a member of the locus.")
  }
  missing_g <- setdiff(members, marker_ids)
  if (length(missing_g) > 0L) {
    .stop_invalid_input("locus members missing from G: %s.",
                        missing_g[1L])
  }

  mk_row <- function(order, snp, cond_set, n_cond, p_raw, p_adj,
                     M, status) {
    data.frame(
      locus_id = locus_id,
      signal_order = order,
      representative_snp = snp,
      conditioning_before_selection = I(list(cond_set)),
      conditioning_key = .trait_set_key(cond_set),
      n_conditioning = n_cond,
      conditional_p_raw = p_raw,
      conditional_p_adjusted = p_adj,
      screening_rule = screening_rule,
      M_remaining = M,
      selection_status = status,
      stringsAsFactors = FALSE
    )
  }

  selected <- lead_snp
  rows <- list(mk_row(1L, lead_snp, character(), 0L,
                      NA_real_, NA_real_, NA_integer_, "selected"))
  history <- list()
  stop_reason <- "max_signals"

  if (max_signals > 1L) {
    step <- 1L
    repeat {
      candidates <- setdiff(members, selected)
      if (length(candidates) == 0L) {
        stop_reason <- "no_remaining_markers"
        break
      }
      M_remaining <- length(candidates)   # fixed before the scan
      X_C <- G[, match(selected, marker_ids), drop = FALSE]
      proj <- .build_conditional_projection(null_fit, X_C,
                                            rank_tol = rank_tol)
      if (!isTRUE(proj$status$ok)) {
        stop_reason <- "projection_failed"
        break
      }
      cscan <- .conditional_mt_scan(
        proj, G[, match(candidates, marker_ids), drop = FALSE],
        marker_ids = candidates, rank_tol = rank_tol,
        return_effects = FALSE
      )
      if (isTRUE(keep_history)) history[[step]] <- cscan$conditional
      pick <- .select_next_signal(
        cscan$conditional, chromosome, position, screening_rule,
        alpha_signal, fixed_p_threshold, M_remaining
      )
      if (!pick$found) {
        rows[[length(rows) + 1L]] <- mk_row(
          step + 1L, NA_character_, selected, length(selected),
          NA_real_, NA_real_, M_remaining, "no_passing_marker"
        )
        stop_reason <- "no_passing_marker"
        break
      }
      rows[[length(rows) + 1L]] <- mk_row(
        step + 1L, pick$marker_id, selected, length(selected),
        pick$p_raw, pick$p_adjusted, M_remaining, "selected"
      )
      selected <- c(selected, pick$marker_id)
      step <- step + 1L
      if (length(selected) >= max_signals) {
        stop_reason <- "max_signals"
        break
      }
    }
  } else {
    stop_reason <- "max_signals"
  }

  signals <- do.call(rbind, rows)
  rownames(signals) <- NULL
  # lead row carries the marginal omnibus p-value separately
  signals$marginal_p <- NA_real_
  signals$marginal_p[1L] <- lead_marginal_p

  list(
    signals = signals,
    selected = selected,
    history = if (isTRUE(keep_history)) history else NULL,
    status = .new_status(ok = TRUE, code = "ok"),
    diagnostics = list(
      stop_reason = stop_reason,
      n_signals = length(selected),
      n_steps = length(history),
      screening_rule = screening_rule,
      alpha_signal = alpha_signal,
      fixed_p_threshold = fixed_p_threshold,
      max_signals = max_signals
    )
  )
}

#' Resolve locus signals and estimate their joint effects
#'
#' Orchestrates the locus-level signal resolution: for each locus of a
#' [define_associated_loci()] result, the lead SNP is the first
#' representative; the stepwise conditional core
#' [`.resolve_one_locus_signals()`] selects further representatives;
#' the final signal-specific effects come from ONE joint GLS refit of
#' the selected set via [`.fit_joint_signal_effects()`] (the same
#' path used by [estimate_mt_effects()] with conditioning sets).
#' Stepwise temporary conditional betas never appear in the output.
#'
#' @param null_fit The single global null fit (with rotation).
#' @param G `n x p` genotype matrix.
#' @param locus_object A [define_associated_loci()] result (or a list
#'   with `loci` and `membership` data.frames).
#' @param marker_ids Marker order of `G`.
#' @param signal_adjust Local screening rule:
#'   `"within_locus_bonferroni"`, `"fixed"` or `"none"`.
#' @param alpha_signal Level for the bonferroni/none rules.
#' @param fixed_p_threshold Raw-p threshold; required when
#'   `signal_adjust = "fixed"`.
#' @param max_signals Maximum signals per locus.
#' @param rank_tol Tolerance forwarded to [`.safe_inverse()`].
#' @param return_conditional_scan Logical; keep the per-step
#'   conditional scan tables.
#'
#' @return A list with components
#' \describe{
#'   \item{signals}{data.frame, one row per resolved signal:
#'     `locus_id`, `signal_order`, `representative_snp`,
#'     `conditioning_key`, `conditional_p_raw`,
#'     `conditional_p_adjusted`, `screening_rule`, `M_remaining`,
#'     `selection_status`, `marginal_p`.}
#'   \item{beta}{Long data.frame: `locus_id`, `signal_order`,
#'     `representative_snp`, `trait`, `beta`, `se` from the final
#'     joint model.}
#'   \item{covariance}{Named list per `locus_id`: the full joint
#'     `(Km) x (Km)` covariance matrix of that locus's signal effects
#'     (blocks ordered by `signal_order`).}
#'   \item{conditional_scan}{Per-locus stepwise scan history, or
#'     `NULL`.}
#'   \item{status}{Standard status list.}
#'   \item{diagnostics}{List with `n_loci`, `n_signals`,
#'     `stop_reasons`, `n_rank_deficient`.}
#' }
#' @export
resolve_locus_signals <- function(
  null_fit,
  G,
  locus_object,
  marker_ids = colnames(G),
  signal_adjust = c("within_locus_bonferroni", "fixed", "none"),
  alpha_signal = 0.05,
  fixed_p_threshold = NULL,
  max_signals = 10L,
  rank_tol = sqrt(.Machine$double.eps),
  return_conditional_scan = FALSE
) {
  signal_adjust <- match.arg(signal_adjust)
  .check_prob(alpha_signal, "alpha_signal")
  if (signal_adjust == "fixed") {
    if (is.null(fixed_p_threshold)) {
      .stop_invalid_input(
        "fixed_p_threshold is required when signal_adjust = \"fixed\"."
      )
    }
    .check_prob(fixed_p_threshold, "fixed_p_threshold")
  }
  if (!is.numeric(max_signals) || length(max_signals) != 1L ||
      !is.finite(max_signals) || max_signals < 1L) {
    .stop_invalid_input("max_signals must be a positive integer.")
  }
  .validate_null_fit(null_fit)
  if (is.null(null_fit$rotation)) {
    .stop_invalid_input(
      "null_fit does not contain the rotation object; re-fit with return_rotation = TRUE."
    )
  }
  if (is.null(marker_ids)) marker_ids <- paste0("M", seq_len(ncol(G)))
  if (anyDuplicated(marker_ids)) {
    .stop_invalid_input("marker_ids must be unique.")
  }
  if (!is.list(locus_object) || is.null(locus_object$loci) ||
      is.null(locus_object$membership)) {
    .stop_invalid_input(
      "locus_object must provide loci and membership data.frames."
    )
  }
  loci <- locus_object$loci
  membership <- locus_object$membership
  need_loci <- c("locus_id", "lead_snp", "lead_p")
  if (!all(need_loci %in% names(loci))) {
    .stop_invalid_input("loci must contain locus_id, lead_snp, lead_p.")
  }
  if (!all(c("locus_id", "marker_id", "position") %in% names(membership))) {
    .stop_invalid_input(
      "membership must contain locus_id, marker_id and position."
    )
  }
  # marker annotation for the deterministic tie-break: chromosome from
  # the locus table, position from the membership table
  position <- stats::setNames(membership$position,
                              membership$marker_id)
  chr_of_locus <- stats::setNames(as.character(loci$chromosome),
                                  loci$locus_id)
  chromosome <- stats::setNames(
    chr_of_locus[membership$locus_id],
    membership$marker_id
  )

  signal_rows <- list()
  beta_rows <- list()
  cov_list <- list()
  scan_hist <- list()
  stop_reasons <- character()
  locus_status_rows <- list()
  n_rank_deficient <- 0L

  for (li in seq_len(nrow(loci))) {
    lid <- loci$locus_id[li]
    lead <- loci$lead_snp[li]
    members <- membership$marker_id[membership$locus_id == lid]
    res <- .resolve_one_locus_signals(
      locus_id = lid, lead_snp = lead, members = members,
      null_fit = null_fit, G = G, marker_ids = marker_ids,
      chromosome = chromosome, position = position,
      lead_marginal_p = loci$lead_p[li],
      screening_rule = signal_adjust,
      alpha_signal = alpha_signal,
      fixed_p_threshold = if (!is.null(fixed_p_threshold)) {
        fixed_p_threshold
      } else {
        alpha_signal
      },
      max_signals = max_signals, rank_tol = rank_tol,
      keep_history = TRUE
    )
    stop_reasons[lid] <- res$diagnostics$stop_reason

    # final joint refit of the selected representatives
    joint <- .fit_joint_signal_effects(
      null_fit, G, res$selected, marker_ids = marker_ids,
      rank_tol = rank_tol
    )
    if (joint$status$code != "ok") {
      n_rank_deficient <- n_rank_deficient + 1L
    }
    cov_list[[lid]] <- joint$joint_covariance
    locus_status_rows[[li]] <- data.frame(
      locus_id = lid,
      stop_reason = res$diagnostics$stop_reason,
      n_signals = length(res$selected),
      numerical_status = joint$status$code,
      stringsAsFactors = FALSE
    )

    sig <- res$signals
    sig <- sig[!is.na(sig$representative_snp), , drop = FALSE]
    signal_rows[[li]] <- sig
    for (k in seq_along(res$selected)) {
      beta_rows[[length(beta_rows) + 1L]] <- data.frame(
        locus_id = lid,
        signal_id = paste0(lid, "::S", k),
        signal_order = k,
        representative_snp = res$selected[k],
        trait = colnames(joint$beta),
        beta = joint$beta[k, ],
        se = joint$se[k, ],
        stringsAsFactors = FALSE
      )
    }
    if (isTRUE(return_conditional_scan)) scan_hist[[lid]] <- res$history
  }

  signals_df <- if (length(signal_rows) > 0L) {
    do.call(rbind, signal_rows)
  } else {
    data.frame()
  }
  if (nrow(signals_df) > 0L) {
    signals_df$signal_id <- paste0(signals_df$locus_id, "::S",
                                   signals_df$signal_order)
    signals_df <- signals_df[, c("locus_id", "signal_id",
                                 setdiff(names(signals_df),
                                         c("locus_id", "signal_id")))]
  }
  beta_df <- if (length(beta_rows) > 0L) {
    do.call(rbind, beta_rows)
  } else {
    data.frame()
  }

  list(
    signals = signals_df,
    beta = beta_df,
    covariance = cov_list,
    conditional_scan = if (isTRUE(return_conditional_scan)) {
      scan_hist
    } else {
      NULL
    },
    status = .new_status(ok = TRUE, code = "ok"),
    diagnostics = list(
      n_loci = nrow(loci),
      n_signals = nrow(signals_df),
      stop_reasons = stop_reasons,
      locus_status = if (length(locus_status_rows) > 0L) {
        do.call(rbind, locus_status_rows)
      } else {
        data.frame()
      },
      n_rank_deficient = n_rank_deficient
    )
  )
}
