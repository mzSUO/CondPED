#' Candidate trait sets per signal via within-signal Holm testing
#'
#' Defines the candidate trait set of every signal from its FINAL
#' joint signal-specific effects. There is no marginal omnibus
#' re-gate: secondary signals found by conditional resolution enter
#' the Holm step even when they are not marginally significant.
#' Formal mode applies the Holm procedure to the per-trait Wald
#' p-values WITHIN each signal (never pooled across signals or loci).
#'
#' @details
#' Per signal `l`, per-trait Wald statistics
#' \deqn{Q_{il} = \widehat\beta_{il}^2 / \widehat S_{l,ii}}
#' give raw p-values; the candidate set is
#' \deqn{\widehat{\mathcal A}_l = \{i: p_{il}^{Holm} \le \alpha_{trait}\}}
#' and the effect breadth \eqn{\widehat B_l = |\widehat{\mathcal A}_l|}.
#'
#' `candidate_mode`:
#' \describe{
#'   \item{`holm_fwer`}{formal within-signal Holm.}
#'   \item{`all_traits`}{oracle/sensitivity: every analysed trait is a
#'     candidate. Does NOT act as a trait oracle (use `predefined`).}
#'   \item{`predefined`}{simulation trait oracle: candidate sets are
#'     taken from `predefined_sets` (named by `signal_id`); every
#'     signal must have an entry.}
#' }
#'
#' @param effects Final joint signal effects: a
#'   [resolve_locus_signals()] result, an [estimate_mt_effects()]
#'   result, or a data.frame with `trait`, `beta`, `se` and marker or
#'   signal id columns. SNP-level inputs map each marker to a
#'   degenerate one-signal locus.
#' @param candidate_mode `"holm_fwer"`, `"all_traits"` or
#'   `"predefined"`.
#' @param alpha_trait Within-signal Holm level.
#' @param predefined_sets Named list (by `signal_id`) of trait vectors
#'   for `candidate_mode = "predefined"`.
#' @param return_all Logical; `TRUE` keeps every trait row,
#'   `FALSE` keeps only candidate rows. Never changes
#'   `candidate_sets`.
#'
#' @return A list with components
#' \describe{
#'   \item{selected_signals}{Character vector of signal ids.}
#'   \item{trait_table}{data.frame: `locus_id`, `signal_id`,
#'     `representative_snp`, `trait`, `beta`, `se`, `Q`, `p_raw`,
#'     `p_adjusted`, `candidate`, `effect_sign`, `adjust_scope`.}
#'   \item{candidate_sets}{Named list (by `signal_id`) of candidate
#'     trait sets; only signals with at least one candidate.}
#'   \item{signal_table}{data.frame: `locus_id`, `signal_id`,
#'     `n_candidates`, `effect_breadth`, `direction_pattern`,
#'     `signal_status`.}
#'   \item{mode}{The `candidate_mode` used.}
#'   \item{status}{Standard status list.}
#'   \item{diagnostics}{List of counts.}
#' }
#' @export
attribute_traits <- function(
  effects,
  candidate_mode = c("holm_fwer", "all_traits", "predefined"),
  alpha_trait = 0.05,
  predefined_sets = NULL,
  return_all = TRUE
) {
  candidate_mode <- match.arg(candidate_mode)
  .check_prob(alpha_trait, "alpha_trait")
  if (!is.logical(return_all) || length(return_all) != 1L || is.na(return_all)) {
    .stop_invalid_input("return_all must be TRUE or FALSE.")
  }

  eff <- .as_signal_effects(effects)
  signals <- unique(eff$signal_id)
  warnings <- character()

  if (candidate_mode == "predefined") {
    if (is.null(predefined_sets) || !is.list(predefined_sets)) {
      .stop_invalid_input(
        "predefined_sets (a named list of trait vectors by signal_id) is required for candidate_mode = \"predefined\"."
      )
    }
    missing_sets <- setdiff(signals, names(predefined_sets))
    if (length(missing_sets) > 0L) {
      .stop_invalid_input(
        "predefined_sets has no entry for signal \"%s\".",
        missing_sets[1L]
      )
    }
  }

  # ---- per-signal candidate sets ----------------------------------------------
  p_adjusted <- rep(NA_real_, nrow(eff))
  candidate <- rep(FALSE, nrow(eff))
  adjust_scope <- switch(candidate_mode,
    holm_fwer = "within_signal",
    all_traits = "none",
    predefined = "predefined"
  )
  na_p <- is.na(eff$p_value)

  for (sid in signals) {
    rows <- which(eff$signal_id == sid)
    if (candidate_mode == "holm_fwer") {
      ok <- rows[!na_p[rows]]
      if (length(ok) > 0L) {
        adj <- stats::p.adjust(eff$p_value[ok], method = "holm")
        p_adjusted[ok] <- adj
        candidate[ok] <- adj <= alpha_trait
      }
    } else if (candidate_mode == "all_traits") {
      p_adjusted[rows] <- eff$p_value[rows]
      candidate[rows] <- TRUE
    } else {
      preset <- predefined_sets[[sid]]
      known <- eff$trait[rows]
      unknown <- setdiff(preset, known)
      if (length(unknown) > 0L) {
        .stop_invalid_input(
          "predefined set for signal \"%s\" contains unknown trait(s): %s.",
          sid, unknown[1L]
        )
      }
      p_adjusted[rows] <- eff$p_value[rows]
      candidate[rows] <- known %in% preset
    }
  }
  if (candidate_mode == "holm_fwer" && any(na_p)) {
    warnings <- c(warnings, sprintf(
      "%d signal x trait p-value(s) are NA; those rows are kept with candidate = FALSE.",
      sum(na_p)
    ))
  }

  Q <- (eff$beta / eff$se)^2
  Q[!is.finite(Q)] <- NA_real_
  trait_table <- data.frame(
    locus_id = eff$locus_id,
    signal_id = eff$signal_id,
    representative_snp = eff$representative_snp,
    trait = eff$trait,
    beta = eff$beta,
    se = eff$se,
    Q = Q,
    p_raw = eff$p_value,
    p_adjusted = p_adjusted,
    candidate = candidate,
    effect_sign = sign(eff$beta),
    adjust_scope = rep(adjust_scope, nrow(eff)),
    stringsAsFactors = FALSE
  )

  cand_rows <- trait_table$candidate %in% TRUE
  candidate_sets <- split(trait_table$trait[cand_rows],
                          trait_table$signal_id[cand_rows])
  candidate_sets <- lapply(candidate_sets, unname)

  # ---- signal table -------------------------------------------------------------
  n_cand <- vapply(signals, function(sid) {
    sum(trait_table$candidate[trait_table$signal_id == sid])
  }, integer(1))
  sig_info <- eff[match(signals, eff$signal_id), , drop = FALSE]
  direction <- vapply(signals, function(sid) {
    B <- n_cand[[sid]]
    if (B == 0L) return("not_applicable")
    if (B == 1L) return("single_trait")
    sgn <- sign(trait_table$beta[trait_table$signal_id == sid &
                                   trait_table$candidate])
    if (all(sgn == sgn[1L])) return("concordant")
    if (B == 2L) return("antagonistic")
    "mixed"
  }, character(1))
  signal_status <- ifelse(n_cand == 0L, "omnibus_only",
                          ifelse(n_cand == 1L, "trait_restricted",
                                 "multi_trait"))
  signal_table <- data.frame(
    locus_id = sig_info$locus_id,
    signal_id = signals,
    n_candidates = unname(n_cand),
    effect_breadth = unname(n_cand),
    direction_pattern = unname(direction),
    signal_status = signal_status,
    stringsAsFactors = FALSE
  )

  if (!return_all) {
    trait_table <- trait_table[cand_rows, , drop = FALSE]
    rownames(trait_table) <- NULL
  }
  if (length(signals) == 0L) {
    warnings <- c(warnings, "No signals to attribute; returning an empty result.")
  }

  list(
    selected_signals = signals,
    trait_table = trait_table,
    candidate_sets = candidate_sets,
    signal_table = signal_table,
    mode = candidate_mode,
    status = .new_status(
      ok = TRUE,
      code = if (length(signals) == 0L) "empty_selection" else "ok",
      warnings = warnings
    ),
    diagnostics = list(
      n_signals = length(signals),
      n_candidates = sum(cand_rows),
      n_omnibus_only = sum(signal_status == "omnibus_only"),
      n_trait_restricted = sum(signal_status == "trait_restricted"),
      n_multi_trait = sum(signal_status == "multi_trait"),
      n_na_effect_p = sum(na_p),
      alpha_trait = if (candidate_mode == "holm_fwer") {
        alpha_trait
      } else {
        NA_real_
      }
    )
  )
}

#' Coerce effect inputs to a per-signal long data.frame
#'
#' Accepts a [resolve_locus_signals()] result (`$beta` long table), an
#' [estimate_mt_effects()] result (`$effects_long`), or a data.frame.
#' SNP-level inputs map each marker to a degenerate one-signal locus:
#' `locus_id = marker_id`, `signal_id = marker_id::S1`,
#' `representative_snp = marker_id`.
#'
#' @param effects Result object or data.frame.
#' @return data.frame with `locus_id`, `signal_id`,
#'   `representative_snp`, `trait`, `beta`, `se`, `p_value`.
#' @keywords internal
.as_signal_effects <- function(effects) {
  if (is.list(effects) && !is.data.frame(effects)) {
    if (!is.null(effects$status) && isFALSE(effects$status$ok)) {
      .stop_invalid_input(
        "effects result has status$ok == FALSE; cannot attribute traits."
      )
    }
    if (!is.null(effects$beta) && is.data.frame(effects$beta) &&
        "representative_snp" %in% names(effects$beta)) {
      # resolve_locus_signals() output
      b <- effects$beta
      sid <- if ("signal_id" %in% names(b)) {
        b$signal_id
      } else {
        paste0(b$locus_id, "::S", b$signal_order)
      }
      df <- data.frame(
        locus_id = b$locus_id,
        signal_id = sid,
        representative_snp = b$representative_snp,
        trait = as.character(b$trait),
        beta = b$beta,
        se = b$se,
        stringsAsFactors = FALSE
      )
    } else if (!is.null(effects$effects_long)) {
      df <- effects$effects_long
    } else {
      .stop_invalid_input(
        "effects must be a resolve_locus_signals() or estimate_mt_effects() result, or a data.frame."
      )
    }
  } else {
    df <- effects
  }
  if (!is.data.frame(df)) {
    .stop_invalid_input("effects must be a result object or a data.frame.")
  }
  if (!all(c("trait", "beta", "se") %in% names(df))) {
    .stop_invalid_input(
      "effects must provide trait, beta and se columns."
    )
  }
  if ("signal_id" %in% names(df) && "locus_id" %in% names(df) &&
      "representative_snp" %in% names(df)) {
    df$signal_id <- as.character(df$signal_id)
    df$locus_id <- as.character(df$locus_id)
    df$representative_snp <- as.character(df$representative_snp)
  } else if ("marker_id" %in% names(df)) {
    df$marker_id <- as.character(df$marker_id)
    df$locus_id <- df$marker_id
    df$signal_id <- paste0(df$marker_id, "::S1")
    df$representative_snp <- df$marker_id
  } else {
    .stop_invalid_input(
      "effects must carry signal ids (locus_id/signal_id/representative_snp) or marker_id."
    )
  }
  df$trait <- as.character(df$trait)
  if (!is.numeric(df$beta) || !is.numeric(df$se)) {
    .stop_invalid_input("effects beta and se must be numeric.")
  }
  key <- paste(df$signal_id, df$trait, sep = "\r")
  if (anyDuplicated(key)) {
    .stop_invalid_input(
      "effects contains duplicated signal_id x trait rows (e.g. \"%s\").",
      gsub("\r", " / ", unique(key[duplicated(key)])[1L])
    )
  }
  if (!"p_value" %in% names(df)) {
    ok_se <- is.finite(df$se) & df$se > 0 & is.finite(df$beta)
    df$p_value <- NA_real_
    df$p_value[ok_se] <- stats::pchisq(
      (df$beta[ok_se] / df$se[ok_se])^2, df = 1, lower.tail = FALSE
    )
  }
  df[, c("locus_id", "signal_id", "representative_snp", "trait",
         "beta", "se", "p_value")]
}

#' Build the per-signal summary table
#'
#' Joins the attribution signal table with resolution diagnostics into
#' the frozen signal-summary fields. Pure bookkeeping; no statistics.
#'
#' @param attribution An [attribute_traits()] result.
#' @param resolution A [resolve_locus_signals()] result (or `NULL`
#'   for SNP-level/predefined paths, when `locus_info` may carry the
#'   per-locus numerical status).
#' @param locus_info Optional data.frame with `locus_id`,
#'   `n_signals`, `numerical_status` used when `resolution` is NULL.
#' @return data.frame with `locus_id`, `signal_id`,
#'   `representative_snp`, `n_signals_in_locus`, `effect_breadth`,
#'   `direction_pattern`, `signal_status`, `numerical_status`.
#' @keywords internal
.build_signal_summary <- function(attribution, resolution = NULL,
                                  locus_info = NULL) {
  st <- attribution$signal_table
  if (nrow(st) == 0L) {
    return(data.frame(
      locus_id = character(), signal_id = character(),
      representative_snp = character(), n_signals_in_locus = integer(),
      effect_breadth = integer(), direction_pattern = character(),
      signal_status = character(), numerical_status = character(),
      stringsAsFactors = FALSE
    ))
  }
  if (!is.null(resolution)) {
    sig <- resolution$signals
    sig$signal_id <- paste0(sig$locus_id, "::S", sig$signal_order)
    n_per_locus <- tapply(sig$representative_snp, sig$locus_id,
                          length)
    lst <- resolution$diagnostics$locus_status
    out <- st
    out$representative_snp <- sig$representative_snp[
      match(st$signal_id, sig$signal_id)]
    out$n_signals_in_locus <- unname(n_per_locus[
      match(st$locus_id, names(n_per_locus))])
    out$numerical_status <- if (is.data.frame(lst) && nrow(lst) > 0L) {
      lst$numerical_status[match(st$locus_id, lst$locus_id)]
    } else {
      "ok"
    }
  } else {
    out <- st
    tt <- attribution$trait_table
    out$representative_snp <- tt$representative_snp[
      match(st$signal_id, tt$signal_id)]
    if (!is.null(locus_info) && nrow(locus_info) > 0L) {
      out$n_signals_in_locus <- locus_info$n_signals[
        match(st$locus_id, locus_info$locus_id)]
      out$numerical_status <- locus_info$numerical_status[
        match(st$locus_id, locus_info$locus_id)]
    } else {
      out$n_signals_in_locus <- 1L
      out$numerical_status <- "ok"
    }
  }
  out$n_signals_in_locus[is.na(out$n_signals_in_locus)] <- 1L
  out$numerical_status[is.na(out$numerical_status)] <- "ok"
  out[, c("locus_id", "signal_id", "representative_snp",
          "n_signals_in_locus", "effect_breadth", "direction_pattern",
          "signal_status", "numerical_status")]
}
