#' Candidate trait sets via within-locus Holm testing
#'
#' Defines the candidate trait set of each omnibus-selected locus. The
#' formal mode applies the Holm procedure to the per-trait Wald
#' p-values WITHIN each locus (never pooled across loci); traits
#' passing at `alpha_trait` form the candidate set
#' \eqn{\widehat{\mathcal A}_l}. `candidate_mode = "all_traits"` is an
#' oracle/sensitivity mode that takes every analysed trait as
#' candidate, still restricted to selected loci.
#'
#' @details
#' Layer 1 selects loci: the `M` valid (non-`NA`) omnibus p-values are
#' adjusted with `omnibus_method` and compared against
#' `alpha_omnibus`. `M` counts tested locus families only. Loci not
#' passing layer 1 never enter the candidate results.
#'
#' Layer 2 (`candidate_mode = "holm_fwer"`): within each selected
#' locus, the locus's own per-trait p-values are Holm-adjusted
#' (`adjust_scope = "within_locus"`) and compared with `alpha_trait`.
#' Locus status: `omnibus_only` (no candidate), `trait_restricted`
#' (exactly one), `multi_trait` (two or more).
#'
#' @param omnibus Omnibus scan result returned by [scan_mt_omnibus()],
#'   or its `omnibus` data.frame (columns `marker_id`, `p_value`, ...).
#' @param effects Effect estimates returned by [estimate_mt_effects()],
#'   or its `effects_long` data.frame (columns `marker_id`, `trait`,
#'   `beta`, `se`, `p_value`).
#' @param omnibus_method Layer-1 adjustment: `"BH"`, `"bonferroni"` or
#'   `"none"`.
#' @param alpha_omnibus Layer-1 omnibus significance level.
#' @param candidate_mode `"holm_fwer"` (formal within-locus Holm) or
#'   `"all_traits"` (oracle: all traits are candidates).
#' @param alpha_trait Within-locus Holm level.
#' @param return_all Logical; `TRUE` keeps every trait row of the
#'   selected loci, `FALSE` keeps only candidate rows. Never changes
#'   `candidate_sets`.
#'
#' @return A list with components
#' \describe{
#'   \item{selected_loci}{Character vector of layer-1 selected marker
#'     ids.}
#'   \item{trait_table}{`data.frame` with columns `marker_id`, `trait`,
#'     `beta`, `se`, `Q`, `p_raw`, `p_adjusted`, `candidate`,
#'     `adjust_scope`.}
#'   \item{candidate_sets}{Named list; per-locus candidate trait sets
#'     (only loci with at least one candidate trait).}
#'   \item{locus_table}{`data.frame` with columns `marker_id`,
#'     `n_candidates`, `locus_status` (`omnibus_only`,
#'     `trait_restricted`, `multi_trait`).}
#'   \item{mode}{The `candidate_mode` used.}
#'   \item{status}{Standard status list.}
#'   \item{diagnostics}{List with `M`, `n_selected`, `n_candidates`,
#'     locus-status counts and anomaly counts.}
#' }
#' @export
attribute_traits <- function(
  omnibus,
  effects,
  omnibus_method = c("BH", "bonferroni", "none"),
  alpha_omnibus = 0.05,
  candidate_mode = c("holm_fwer", "all_traits"),
  alpha_trait = 0.05,
  return_all = TRUE
) {
  omnibus_method <- match.arg(omnibus_method)
  candidate_mode <- match.arg(candidate_mode)
  .check_prob(alpha_omnibus, "alpha_omnibus")
  .check_prob(alpha_trait, "alpha_trait")
  if (!is.logical(return_all) || length(return_all) != 1L || is.na(return_all)) {
    .stop_invalid_input("return_all must be TRUE or FALSE.")
  }

  om <- .as_omnibus_table(omnibus)
  eff <- .as_effects_table(effects)
  warnings <- character()

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
  selected_loci <- om$marker_id[valid][p_adj1 <= alpha_omnibus]
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
        "%d selected locus/loci have no rows in effects and cannot ",
        "receive candidate traits (e.g. \"%s\")."
      ),
      length(missing_loci), missing_loci[1L]
    ))
  }
  testable_loci <- intersect(selected_loci, names(eff_idx))

  rows <- unlist(eff_idx[testable_loci], use.names = FALSE)
  sub <- eff[rows, , drop = FALSE]
  p_raw <- sub$p_value
  na_p <- is.na(p_raw)
  n_na_effect_p <- sum(na_p)

  # ---- layer 2: within-locus Holm or oracle ----------------------------------
  p_adjusted <- rep(NA_real_, length(p_raw))
  candidate <- rep(FALSE, length(p_raw))
  if (candidate_mode == "holm_fwer") {
    if (n_na_effect_p > 0L) {
      warnings <- c(warnings, sprintf(
        paste0(
          "%d locus x trait p-value(s) are NA; these rows are kept with ",
          "candidate = FALSE and excluded from their locus's Holm ",
          "adjustment."
        ),
        n_na_effect_p
      ))
    }
    fam <- split(seq_len(nrow(sub)), sub$marker_id)
    for (fi in fam) {
      ok <- fi[!na_p[fi]]
      if (length(ok) == 0L) next
      adj <- stats::p.adjust(p_raw[ok], method = "holm")
      p_adjusted[ok] <- adj
      candidate[ok] <- adj <= alpha_trait
    }
    adjust_scope <- "within_locus"
  } else {
    # all_traits oracle: no testing, every analysed trait is a candidate
    p_adjusted <- p_raw
    candidate <- rep(TRUE, length(p_raw))
    adjust_scope <- "none"
  }

  Q <- (sub$beta / sub$se)^2
  Q[!is.finite(Q)] <- NA_real_

  trait_table <- data.frame(
    marker_id = sub$marker_id,
    trait = sub$trait,
    beta = sub$beta,
    se = sub$se,
    Q = Q,
    p_raw = p_raw,
    p_adjusted = p_adjusted,
    candidate = candidate,
    adjust_scope = rep(adjust_scope, nrow(sub)),
    stringsAsFactors = FALSE
  )

  # candidate_sets: only loci with at least one candidate trait
  cand_rows <- trait_table$candidate %in% TRUE
  candidate_sets <- split(trait_table$trait[cand_rows],
                          trait_table$marker_id[cand_rows])
  candidate_sets <- lapply(candidate_sets, unname)

  # locus_table: every selected locus with effects rows
  n_cand <- vapply(testable_loci, function(mk) {
    sum(trait_table$candidate[trait_table$marker_id == mk])
  }, integer(1))
  locus_status <- ifelse(n_cand == 0L, "omnibus_only",
                         ifelse(n_cand == 1L, "trait_restricted",
                                "multi_trait"))
  locus_table <- data.frame(
    marker_id = testable_loci,
    n_candidates = unname(n_cand),
    locus_status = locus_status,
    stringsAsFactors = FALSE
  )

  if (!return_all) {
    trait_table <- trait_table[cand_rows, , drop = FALSE]
    rownames(trait_table) <- NULL
  }

  if (n_selected == 0L) {
    warnings <- c(warnings,
                  "Layer 1 selected no loci; returning an empty result.")
  }

  list(
    selected_loci = selected_loci,
    trait_table = trait_table,
    candidate_sets = candidate_sets,
    locus_table = locus_table,
    mode = candidate_mode,
    status = .new_status(ok = TRUE, code = "ok", warnings = warnings),
    diagnostics = list(
      M = M,
      n_selected = n_selected,
      n_candidates = sum(cand_rows),
      n_omnibus_only = sum(locus_status == "omnibus_only"),
      n_trait_restricted = sum(locus_status == "trait_restricted"),
      n_multi_trait = sum(locus_status == "multi_trait"),
      n_omnibus_rows = nrow(om),
      n_na_omnibus_p = n_na_omnibus_p,
      n_missing_effect_loci = length(missing_loci),
      n_na_effect_p = n_na_effect_p,
      alpha_omnibus = alpha_omnibus,
      alpha_trait = if (candidate_mode == "holm_fwer") {
        alpha_trait
      } else {
        NA_real_
      }
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
