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
