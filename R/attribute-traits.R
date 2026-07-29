#' Trait attribution with hierarchical error control
#'
#' Attributes the omnibus-selected loci to individual traits, with
#' hierarchical multiplicity correction (Benjamini-Bogomolov style FDR or
#' Holm FWER) layered on top of the first-layer omnibus selection.
#'
#' @details
#' The first layer selects loci from the omnibus scan: the `M` valid
#' (non-`NA`) omnibus p-values are adjusted with `omnibus_method` and
#' compared against the first-layer level (`alpha_omnibus`, or
#' `alpha_split[1]` in `"holm_fwer"` mode). `M` counts tested locus
#' families only -- neither all rows of the omnibus table nor the number
#' of selected loci.
#'
#' The second layer attributes traits within the selected set `S`:
#' \describe{
#'   \item{`bb_fdr`}{Benjamini-Bogomolov style: within each selected locus
#'     family, BH is applied to that family's per-trait p-values at level
#'     `q_within = q_target * |S| / M` (Methods eq. 30). The per-trait
#'     hypotheses of different loci are never pooled into a single BH.}
#'   \item{`holm_fwer`}{the error budget `alpha_total` is split according
#'     to `alpha_split`; all locus x trait hypotheses actually tested in
#'     the second layer are pooled and adjusted with Holm at level
#'     `alpha_split[2]`. For the confirmatory FWER mode of Methods
#'     section 3.7.2, use `omnibus_method = "bonferroni"` so that the
#'     first layer spends `alpha_split[1]` as `alpha_1 / M`.}
#'   \item{`none`}{MT-Posthoc baseline: per-trait p-values are reported
#'     unchanged and `attributed` is decided by the explicit unadjusted
#'     rule `p_raw <= q_target`.}
#' }
#'
#' `dependence = "empirical"` is accepted for interface compatibility, but
#' the interface contract specifies no permutation/bootstrap calibration
#' algorithm, so none is performed: the procedure is identical to
#' `"assumed"` and a warning plus an `assumptions` entry record that the
#' error-rate labels are unverified under the empirical dependence
#' structure.
#'
#' @param omnibus Omnibus scan result returned by [scan_mt_omnibus()], or
#'   its `omnibus` data.frame (columns `marker_id`, `p_value`, ...).
#' @param effects Effect estimates returned by [estimate_mt_effects()], or
#'   its `effects_long` data.frame (columns `marker_id`, `trait`, `beta`,
#'   `se`, `p_value`).
#' @param omnibus_method Multiple-testing adjustment for the first-layer
#'   omnibus p-values: `"BH"`, `"bonferroni"` or `"none"`.
#' @param alpha_omnibus First-layer omnibus significance level (modes
#'   `"bb_fdr"` and `"none"`).
#' @param attribution_mode Second-layer correction mode: `"bb_fdr"`
#'   (Benjamini-Bogomolov style FDR on the selected families), `"holm_fwer"`
#'   (split error budget with Holm on all subsequent hypotheses) or
#'   `"none"`.
#' @param q_target Target FDR level for `"bb_fdr"`; in `"none"` mode the
#'   nominal unadjusted per-trait threshold.
#' @param alpha_total Total error budget for `"holm_fwer"`.
#' @param alpha_split Numeric vector of length 2 giving the split of
#'   `alpha_total` between the first and subsequent layers; must sum to
#'   `alpha_total`.
#' @param dependence Dependence assumption for the BB step: `"assumed"` or
#'   `"empirical"` (see Details; no empirical calibration is performed).
#' @param return_all Logical; return the full per-trait table for the
#'   selected loci, including non-attributed rows. `FALSE` keeps only the
#'   attributed rows and never changes `A`.
#'
#' @return A list with components
#' \describe{
#'   \item{selected_loci}{Character vector of first-layer selected marker
#'     ids.}
#'   \item{trait_table}{`data.frame` with columns `marker_id`, `trait`,
#'     `beta`, `se`, `p_raw`, `p_adjusted`, `attributed`.}
#'   \item{A}{Named list; per-locus associated-trait sets
#'     \eqn{\widehat{\mathcal A}_l} (Methods eq. 18), exactly the rows of
#'     `trait_table` with `attributed == TRUE` split by `marker_id`. Loci
#'     with no attributed trait are absent (omnibus-only/unresolved).}
#'   \item{mode}{The `attribution_mode` used.}
#'   \item{assumptions}{Character vector recording the selection rule, the
#'     dependence/selection assumptions of the BB procedure and the
#'     threshold rule of `"none"` mode.}
#'   \item{status}{Standard status list (`ok`, `code`, `message`,
#'     `warnings`).}
#'   \item{diagnostics}{List with `M`, `n_selected`, `q_within`,
#'     `n_attributed` and additional counts (`n_na_omnibus_p`,
#'     `n_missing_effect_loci`, `n_na_effect_p`, layer levels).}
#' }
#' @export
attribute_traits <- function(
  omnibus,
  effects,
  omnibus_method = c("BH", "bonferroni", "none"),
  alpha_omnibus = 0.05,
  attribution_mode = c("bb_fdr", "holm_fwer", "none"),
  q_target = 0.05,
  alpha_total = 0.05,
  alpha_split = c(0.025, 0.025),
  dependence = c("assumed", "empirical"),
  return_all = TRUE
) {
  omnibus_method <- match.arg(omnibus_method)
  attribution_mode <- match.arg(attribution_mode)
  dependence <- match.arg(dependence)
  .check_prob(alpha_omnibus, "alpha_omnibus")
  .check_prob(q_target, "q_target")
  .check_prob(alpha_total, "alpha_total")
  if (!is.numeric(alpha_split) || length(alpha_split) != 2L ||
      any(!is.finite(alpha_split)) || any(alpha_split <= 0)) {
    .stop_invalid_input(
      "alpha_split must be a numeric vector of length 2 with positive entries."
    )
  }
  if (attribution_mode == "holm_fwer" &&
      abs(sum(alpha_split) - alpha_total) > 1e-8) {
    .stop_invalid_input(
      paste0(
        "alpha_split must sum to alpha_total for attribution_mode = ",
        "\"holm_fwer\" (got sum(alpha_split) = %.6g, alpha_total = %.6g)."
      ),
      sum(alpha_split), alpha_total
    )
  }
  if (!is.logical(return_all) || length(return_all) != 1L || is.na(return_all)) {
    .stop_invalid_input("return_all must be TRUE or FALSE.")
  }

  om <- .as_omnibus_table(omnibus)
  eff <- .as_effects_table(effects)

  warnings <- character()
  if (dependence == "empirical") {
    warnings <- c(warnings, paste0(
      "dependence = \"empirical\": the interface contract specifies no ",
      "permutation/bootstrap calibration algorithm, so none was performed; ",
      "the procedure is identical to dependence = \"assumed\" and the ",
      "error-rate labels are unverified under the empirical dependence ",
      "structure."
    ))
  }

  # ---- layer 1: omnibus selection ------------------------------------------
  # M counts valid (non-NA) tested locus families only.
  valid <- !is.na(om$p_value)
  M <- sum(valid)
  n_na_omnibus_p <- sum(!valid)
  if (M == 0L) {
    .stop_invalid_input(
      "omnibus contains no valid (non-NA) p-values; nothing can be selected."
    )
  }
  p_valid <- om$p_value[valid]
  p_adj1 <- switch(omnibus_method,
    BH = stats::p.adjust(p_valid, method = "BH"),
    bonferroni = stats::p.adjust(p_valid, method = "bonferroni"),
    none = p_valid
  )
  alpha_layer1 <- if (attribution_mode == "holm_fwer") {
    alpha_split[1L]
  } else {
    alpha_omnibus
  }
  selected_loci <- om$marker_id[valid][p_adj1 <= alpha_layer1]
  n_selected <- length(selected_loci)

  # ---- layer 2 bookkeeping ---------------------------------------------------
  key <- paste(eff$marker_id, eff$trait, sep = "\r")
  if (anyDuplicated(key)) {
    dup <- unique(key[duplicated(key)])
    .stop_invalid_input(
      paste0(
        "effects contains duplicated marker_id x trait keys (e.g. \"%s\"); ",
        "each locus x trait hypothesis must appear exactly once."
      ),
      gsub("\r", " / ", dup[1L])
    )
  }
  eff_idx <- split(seq_len(nrow(eff)), eff$marker_id)
  missing_loci <- setdiff(selected_loci, names(eff_idx))
  if (length(missing_loci) > 0L) {
    warnings <- c(warnings, sprintf(
      paste0(
        "%d first-layer selected locus/loci have no rows in effects and ",
        "cannot enter A (e.g. \"%s\")."
      ),
      length(missing_loci), missing_loci[1L]
    ))
  }
  testable_loci <- intersect(selected_loci, names(eff_idx))

  # Collect second-layer rows in omnibus selection order.
  rows <- unlist(eff_idx[testable_loci], use.names = FALSE)
  sub <- eff[rows, , drop = FALSE]
  p_raw <- sub$p_value
  na_p <- is.na(p_raw)
  n_na_effect_p <- sum(na_p)
  if (n_na_effect_p > 0L) {
    warnings <- c(warnings, sprintf(
      paste0(
        "%d locus x trait p-value(s) are NA; these rows are kept with ",
        "attributed = FALSE and excluded from their family's adjustment."
      ),
      n_na_effect_p
    ))
  }

  # ---- layer 2: correction ---------------------------------------------------
  q_within <- NA_real_
  alpha_layer2 <- NA_real_
  p_adjusted <- rep(NA_real_, length(p_raw))
  attributed <- rep(FALSE, length(p_raw))

  if (n_selected > 0L) {
    if (attribution_mode == "bb_fdr") {
      q_within <- q_target * n_selected / M
      fam <- split(seq_len(nrow(sub)), sub$marker_id)
      for (fi in fam) {
        ok <- fi[!na_p[fi]]
        if (length(ok) == 0L) next
        adj <- stats::p.adjust(p_raw[ok], method = "BH")
        p_adjusted[ok] <- adj
        attributed[ok] <- adj <= q_within
      }
    } else if (attribution_mode == "holm_fwer") {
      alpha_layer2 <- alpha_split[2L]
      ok <- which(!na_p)
      if (length(ok) > 0L) {
        adj <- stats::p.adjust(p_raw[ok], method = "holm")
        p_adjusted[ok] <- adj
        attributed[ok] <- adj <= alpha_layer2
      }
    } else { # none: MT-Posthoc baseline, no adjustment of any kind
      p_adjusted <- p_raw
      attributed <- !na_p & p_raw <= q_target
    }
  } else {
    if (attribution_mode == "bb_fdr") q_within <- q_target * 0 / M
    if (attribution_mode == "holm_fwer") alpha_layer2 <- alpha_split[2L]
    warnings <- c(warnings,
                  "First layer selected no loci; returning an empty result.")
  }

  trait_table <- data.frame(
    marker_id = sub$marker_id,
    trait = sub$trait,
    beta = sub$beta,
    se = sub$se,
    p_raw = p_raw,
    p_adjusted = p_adjusted,
    attributed = attributed,
    stringsAsFactors = FALSE
  )

  # A: exactly the attributed rows, split by marker.
  att_rows <- trait_table$attributed %in% TRUE
  A <- split(trait_table$trait[att_rows],
             trait_table$marker_id[att_rows])
  A <- lapply(A, unname)

  if (!return_all) {
    trait_table <- trait_table[att_rows, , drop = FALSE]
    rownames(trait_table) <- NULL
  }

  assumptions <- switch(attribution_mode,
    bb_fdr = c(
      sprintf(paste0(
        "Layer 1: omnibus p-values adjusted by \"%s\" over M = %d valid ",
        "tested locus families at level %.6g; |S| = %d loci selected."
      ), omnibus_method, M, alpha_layer1, n_selected),
      sprintf(paste0(
        "Layer 2 (Benjamini-Bogomolov): within each selected family, BH at ",
        "q_within = q_target * |S| / M = %.6g; per-trait hypotheses are ",
        "never pooled across loci."
      ), q_within),
      paste0(
        "FDR guarantee assumes the layer-1 selection rule is simple/stable ",
        "in the Benjamini-Bogomolov sense and that the selected families ",
        "satisfy the independence/PRDS dependence conditions (Methods ",
        "section 3.7.2); these conditions are assumed, not verified."
      )
    ),
    holm_fwer = c(
      sprintf(paste0(
        "Layer 1: omnibus p-values adjusted by \"%s\" over M = %d valid ",
        "tested locus families at level alpha_split[1] = %.6g; |S| = %d ",
        "loci selected. Use omnibus_method = \"bonferroni\" for the ",
        "confirmatory alpha_1 / M first layer of Methods section 3.7.2."
      ), omnibus_method, M, alpha_layer1, n_selected),
      sprintf(paste0(
        "Layer 2: all %d locus x trait hypotheses actually tested are ",
        "pooled and Holm-adjusted at level alpha_split[2] = %.6g; the ",
        "total budget alpha_total = %.6g is split as (%.6g, %.6g)."
      ), sum(!na_p), alpha_layer2, alpha_total, alpha_split[1L],
      alpha_split[2L])
    ),
    none = c(
      sprintf(paste0(
        "Layer 1: omnibus p-values adjusted by \"%s\" over M = %d valid ",
        "tested locus families at level %.6g; |S| = %d loci selected."
      ), omnibus_method, M, alpha_layer1, n_selected),
      sprintf(paste0(
        "Layer 2 (none / MT-Posthoc baseline): per-trait p-values are ",
        "reported unadjusted (p_adjusted == p_raw); attributed is decided ",
        "by the explicit unadjusted rule p_raw <= q_target = %.6g. No ",
        "error-rate guarantee is claimed."
      ), q_target)
    )
  )
  if (dependence == "empirical") {
    assumptions <- c(assumptions, paste0(
      "dependence = \"empirical\" was requested but no calibration ",
      "algorithm is specified by the contract; results are identical to ",
      "\"assumed\" and dependence-sensitive error rates are unverified."
    ))
  }

  list(
    selected_loci = selected_loci,
    trait_table = trait_table,
    A = A,
    mode = attribution_mode,
    assumptions = assumptions,
    status = .new_status(ok = TRUE, code = "ok", warnings = warnings),
    diagnostics = list(
      M = M,
      n_selected = n_selected,
      q_within = q_within,
      n_attributed = sum(att_rows),
      n_omnibus_rows = nrow(om),
      n_na_omnibus_p = n_na_omnibus_p,
      n_missing_effect_loci = length(missing_loci),
      n_na_effect_p = n_na_effect_p,
      alpha_layer1 = alpha_layer1,
      alpha_layer2 = alpha_layer2
    )
  )
}

#' Coerce the omnibus input to a validated data.frame
#'
#' Accepts either a [scan_mt_omnibus()] result object or its `omnibus`
#' data.frame directly.
#'
#' @param omnibus Result object or data.frame.
#' @return data.frame with at least `marker_id` and `p_value`.
#' @keywords internal
.as_omnibus_table <- function(omnibus) {
  if (is.list(omnibus) && !is.data.frame(omnibus)) {
    if (!is.null(omnibus$status) && isFALSE(omnibus$status$ok)) {
      .stop_invalid_input(
        "omnibus result has status$ok == FALSE; refusing to attribute traits."
      )
    }
    omnibus <- omnibus$omnibus
  }
  if (!is.data.frame(omnibus)) {
    .stop_invalid_input(
      "omnibus must be a scan_mt_omnibus() result or a data.frame."
    )
  }
  need <- c("marker_id", "p_value")
  missing_cols <- setdiff(need, names(omnibus))
  if (length(missing_cols) > 0L) {
    .stop_invalid_input(
      "omnibus data.frame is missing column(s): %s.",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (!is.numeric(omnibus$p_value) ||
      any(!is.na(omnibus$p_value) & (omnibus$p_value < 0 | omnibus$p_value > 1))) {
    .stop_invalid_input("omnibus$p_value must be numeric values in [0, 1].")
  }
  omnibus$marker_id <- as.character(omnibus$marker_id)
  omnibus
}

#' Coerce the effects input to a validated data.frame
#'
#' Accepts either an [estimate_mt_effects()] result object or its
#' `effects_long` data.frame directly.
#'
#' @param effects Result object or data.frame.
#' @return data.frame with at least `marker_id`, `trait`, `beta`, `se`,
#'   `p_value`.
#' @keywords internal
.as_effects_table <- function(effects) {
  if (is.list(effects) && !is.data.frame(effects)) {
    if (!is.null(effects$status) && isFALSE(effects$status$ok)) {
      .stop_invalid_input(
        "effects result has status$ok == FALSE; refusing to attribute traits."
      )
    }
    effects <- effects$effects_long
  }
  if (!is.data.frame(effects)) {
    .stop_invalid_input(
      "effects must be an estimate_mt_effects() result or a data.frame."
    )
  }
  need <- c("marker_id", "trait", "beta", "se", "p_value")
  missing_cols <- setdiff(need, names(effects))
  if (length(missing_cols) > 0L) {
    .stop_invalid_input(
      "effects data.frame is missing column(s): %s.",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (!is.numeric(effects$p_value) ||
      any(!is.na(effects$p_value) &
          (effects$p_value < 0 | effects$p_value > 1))) {
    .stop_invalid_input("effects$p_value must be numeric values in [0, 1].")
  }
  effects$marker_id <- as.character(effects$marker_id)
  effects$trait <- as.character(effects$trait)
  effects
}

#' Check that an argument is a single probability in (0, 1)
#'
#' @param x Value to check.
#' @param name Argument name for the error message.
#' @keywords internal
.check_prob <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    .stop_invalid_input("%s must be a single number in (0, 1).", name)
  }
  invisible(x)
}
