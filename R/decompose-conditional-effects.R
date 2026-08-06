#' Conditional representation loss and exact trait-subset summaries
#'
#' The public name is kept for package compatibility; the v1.0 duty of
#' this function is "conditional representation loss and trait subset
#' summaries". For every locus it evaluates the representation loss of
#' candidate representing subsets with the unique Stage-1 loss path
#' [`.compute_subset_loss()`], and in `subset_mode = "all"` extracts
#' the minimum representative trait sets and the irreducible trait
#' modules at each tolerance. No conditional P-values, no profile or
#' partition objects are produced.
#'
#' @details
#' Per-locus computation order (frozen): normalise the candidate set
#' A, evaluate the full quadratic form, enumerate the subsets, compute
#' Gamma/eta/Omega and the three quadratic forms per subset, check the
#' decomposition and range, run the full-table monotonicity check
#' (`"all"` mode only), then extract all minimum representative sets
#' and all irreducible modules per tolerance and build the tolerance
#' path. Only `"all"` yields formal set-valued results; `"singleton"`
#' is the restricted Simulation-III baseline and `"custom"` evaluates
#' user-supplied sets (both add the empty and full sets; neither runs
#' the global monotonicity check or extraction).
#'
#' @param effects Locus effects: an [estimate_mt_effects()] result or
#'   its `effects_long` data.frame (columns `marker_id`, `trait`,
#'   `beta`).
#' @param attribution Candidate trait sets: an [attribute_traits()]
#'   result (its `candidate_sets` or legacy `A` component) or a named
#'   list mapping marker ids to trait vectors.
#' @param contrasts A `"condped_representation_basis"` object from
#'   [derive_conditional_contrasts()] (or a list with `Sigma_P_ref`
#'   and `trait_names`).
#' @param restrict_to_candidates Logical; `TRUE` (default, formal mode)
#'   analyses each locus's candidate trait set; `FALSE` uses all
#'   available traits (oracle/diagnostic only).
#' @param subset_mode `"all"` (enumerate all `2^k` subsets),`
#'   `"singleton"` (empty, full, singletons and their complements) or
#'   `"custom"` (`custom_sets` plus the empty and full sets).
#' @param custom_sets List of trait sets for `subset_mode = "custom"`.
#' @param tolerance Primary tolerance tau in [0, 1).
#' @param sensitivity_tolerance Additional tolerances for the
#'   sensitivity path.
#' @param max_traits_exact Maximum candidate-set size for exact
#'   enumeration; larger candidate sets in `"all"` mode are reported
#'   as `too_many_traits` and never silently demoted.
#' @param return_subset_table Logical; `FALSE` drops the (potentially
#'   large) subset table from the return value.
#'
#' @return A list with components `subset_table` (or `NULL`),
#'   `minimum_representative_sets`, `irreducible_modules`,
#'   `tolerance_path`, `status` and `diagnostics` (with
#'   `max_decomposition_error`, `max_monotonicity_violation`,
#'   `n_rank_deficient`, `n_unstable` and locus counts).
#' @export
decompose_conditional_effects <- function(
  effects,
  attribution,
  contrasts,
  restrict_to_candidates = TRUE,
  subset_mode = c("all", "singleton", "custom"),
  custom_sets = NULL,
  tolerance = 0.10,
  sensitivity_tolerance = c(0.05, 0.10, 0.20),
  max_traits_exact = 10L,
  return_subset_table = TRUE
) {
  subset_mode <- match.arg(subset_mode)
  if (!is.logical(restrict_to_candidates) ||
      length(restrict_to_candidates) != 1L || is.na(restrict_to_candidates)) {
    .stop_invalid_input("restrict_to_candidates must be TRUE or FALSE.")
  }
  .check_prob01 <- function(x, name) {
    ok <- is.numeric(x) && all(is.finite(x)) &&
      all(x >= 0) && all(x < 1)
    if (!ok) .stop_invalid_input("%s must contain values in [0, 1).", name)
  }
  .check_prob01(tolerance, "tolerance")
  .check_prob01(sensitivity_tolerance, "sensitivity_tolerance")
  if (!is.numeric(max_traits_exact) || length(max_traits_exact) != 1L ||
      !is.finite(max_traits_exact) || max_traits_exact < 1) {
    .stop_invalid_input("max_traits_exact must be a positive integer.")
  }
  if (!is.logical(return_subset_table) ||
      length(return_subset_table) != 1L || is.na(return_subset_table)) {
    .stop_invalid_input("return_subset_table must be TRUE or FALSE.")
  }

  # ---- inputs ---------------------------------------------------------------
  if (!is.list(contrasts) || is.null(contrasts$Sigma_P_ref) ||
      is.null(contrasts$trait_names)) {
    .stop_invalid_input(
      paste0(
        "contrasts must be a derive_conditional_contrasts() basis object ",
        "(list with Sigma_P_ref and trait_names)."
      )
    )
  }
  Sigma_ref <- contrasts$Sigma_P_ref
  trait_names <- contrasts$trait_names
  inverse_tol <- if (!is.null(contrasts$basis$inverse_tol)) {
    contrasts$basis$inverse_tol
  } else {
    sqrt(.Machine$double.eps)
  }
  monotonicity_tol <- 1e-10   # frozen default (contract section 5.2)

  beta_list <- .effects_to_beta_list(effects)
  cand_sets <- .attribution_to_candidate_sets(
    attribution, required = restrict_to_candidates
  )
  markers <- names(beta_list)
  analysis_scope <- if (restrict_to_candidates) "candidates" else "all_traits"
  # custom sets are defined on the global trait space; per locus they are
  # intersected with that locus's candidate set
  if (subset_mode == "custom") {
    if (is.null(custom_sets) || length(custom_sets) == 0L) {
      .stop_invalid_input(
        "custom_sets must be a non-empty list when subset_mode = \"custom\"."
      )
    }
    custom_sets <- lapply(custom_sets, .normalize_trait_set,
                          trait_names, "custom_sets")
  }

  # ---- empty-table prototypes ------------------------------------------------
  empty_subset <- data.frame(
    marker_id = character(), set_id = integer(),
    analysis_scope = character(),
    representing_set = I(list()), representing_key = character(),
    complement_set = I(list()), complement_key = character(),
    set_size = integer(),
    full_qform = numeric(), subset_qform = numeric(),
    residual_qform = numeric(), representation_loss = numeric(),
    decomposition_error = numeric(), feasible_primary = logical(),
    rank_Sigma_SS = integer(), rank_Omega = integer(),
    condition_Sigma_SS = numeric(), condition_Omega = numeric(),
    used_pseudoinverse = logical(), status = character(),
    stringsAsFactors = FALSE
  )
  empty_reps <- data.frame(
    marker_id = character(), tolerance = numeric(),
    solution_id = integer(), trait_set = I(list()),
    trait_key = character(), set_size = integer(),
    representation_loss = numeric(), n_tied_solutions = integer(),
    status = character(), stringsAsFactors = FALSE
  )
  empty_mods <- data.frame(
    marker_id = character(), tolerance = numeric(),
    module_id = integer(), trait_set = I(list()),
    trait_key = character(), module_size = integer(),
    complement_set = I(list()), complement_key = character(),
    complement_loss = numeric(), status = character(),
    stringsAsFactors = FALSE
  )
  empty_path <- data.frame(
    marker_id = character(), tolerance = numeric(),
    minimum_set_size = integer(), n_minimum_sets = integer(),
    n_irreducible_modules = integer(), status = character(),
    stringsAsFactors = FALSE
  )

  tolerances <- unique(c(tolerance, sensitivity_tolerance))
  subset_rows <- list()
  reps_list <- list()
  mods_list <- list()
  path_list <- list()
  n_empty_candidates <- 0L
  n_single_candidate <- 0L
  n_too_many <- 0L
  n_unstable_loci <- 0L
  n_rank_def_rows <- 0L
  max_decomp <- 0
  max_mono <- 0

  for (mk in markers) {
    beta_l <- beta_list[[mk]]
    if (restrict_to_candidates) {
      A_l <- if (mk %in% names(cand_sets)) cand_sets[[mk]] else character()
      A_l <- .normalize_trait_set(A_l, trait_names, "candidate set")
    } else {
      A_l <- .normalize_trait_set(names(beta_l), trait_names,
                                  "effects traits")
    }
    k <- length(A_l)

    if (k == 0L) {
      n_empty_candidates <- n_empty_candidates + 1L
      next
    }

    locus_status <- "ok"
    if (subset_mode == "all" && k > max_traits_exact) {
      n_too_many <- n_too_many + 1L
      for (tol in tolerances) {
        path_list[[length(path_list) + 1L]] <- data.frame(
          marker_id = mk, tolerance = tol,
          minimum_set_size = NA_integer_, n_minimum_sets = 0L,
          n_irreducible_modules = 0L, status = "too_many_traits",
          stringsAsFactors = FALSE
        )
      }
      next
    }

    # Evaluate every enumerated subset through the unique loss path.
    Sigma_AA <- Sigma_ref[A_l, A_l, drop = FALSE]
    beta_A <- beta_l[A_l]
    custom_l <- if (subset_mode == "custom") {
      lapply(custom_sets, function(s) A_l[A_l %in% s])
    } else {
      NULL
    }
    enum <- .enumerate_subsets(A_l, subset_mode, custom_l)
    rows <- vector("list", length(enum$sets))
    for (si in seq_along(enum$sets)) {
      out <- .compute_subset_loss(
        beta_A, Sigma_AA, enum$sets[[si]],
        trait_names = A_l, inverse_tol = inverse_tol
      )
      rows[[si]] <- data.frame(
        marker_id = mk, set_id = si,
        analysis_scope = analysis_scope,
        representing_set = I(list(out$representing_set)),
        representing_key = out$representing_key,
        complement_set = I(list(out$complement_set)),
        complement_key = out$complement_key,
        set_size = length(out$representing_set),
        full_qform = out$full_qform,
        subset_qform = out$subset_qform,
        residual_qform = out$residual_qform,
        representation_loss = out$representation_loss,
        decomposition_error = out$decomposition_error,
        feasible_primary = !is.na(out$representation_loss) &&
          out$representation_loss <= tolerance + 1e-10,
        rank_Sigma_SS = out$rank_Sigma_SS,
        rank_Omega = out$rank_Omega,
        condition_Sigma_SS = out$condition_Sigma_SS,
        condition_Omega = out$condition_Omega,
        used_pseudoinverse = out$used_pseudoinverse,
        status = out$status$code,
        stringsAsFactors = FALSE
      )
    }
    tab <- do.call(rbind, rows)
    max_decomp <- max(max_decomp,
                      max(tab$decomposition_error, na.rm = TRUE, 0))
    n_rank_def_rows <- n_rank_def_rows +
      sum(tab$status == "rank_deficient")

    if (k == 1L) {
      # trait_restricted: no combinatorial optimisation (contract 3.3)
      n_single_candidate <- n_single_candidate + 1L
      tab$status <- "not_applicable"
      tab$feasible_primary <- FALSE
      locus_status <- "not_applicable"
    } else if (any(tab$status != "ok")) {
      # any unusable row (rank deficiency / failed decomposition / no
      # signal) makes formal extraction unsafe for this locus
      locus_status <- if (any(tab$status == "not_applicable")) {
        "not_applicable"
      } else {
        "unstable"
      }
    } else if (subset_mode == "all") {
      mono <- .check_loss_monotonicity(tab, monotonicity_tol)
      max_mono <- max(max_mono, mono$max_violation)
      if (mono$n_violations > 0L) {
        locus_status <- "unstable"
        tab$status[tab$status == "ok"] <- "unstable"
      }
    }
    if (locus_status == "unstable") {
      n_unstable_loci <- n_unstable_loci + 1L
    }

    subset_rows[[length(subset_rows) + 1L]] <- tab

    # Formal extraction: complete enumeration, usable locus, k >= 2.
    do_extract <- subset_mode == "all" && k >= 2L && locus_status == "ok"
    if (do_extract) {
      mr <- stats::setNames(vector("list", length(tolerances)),
                            format(tolerances))
      mo <- mr
      for (tol in tolerances) {
        key <- format(tol)
        r <- .extract_minimum_representative_sets(tab, tol)
        m <- .extract_irreducible_modules(tab, tol)
        mr[[key]] <- r
        mo[[key]] <- m
        if (nrow(r) > 0L) {
          reps_list[[length(reps_list) + 1L]] <-
            cbind(marker_id = mk, r, stringsAsFactors = FALSE)
        }
        if (nrow(m) > 0L) {
          mods_list[[length(mods_list) + 1L]] <-
            cbind(marker_id = mk, m, stringsAsFactors = FALSE)
        }
      }
      path <- .build_tolerance_path(tolerances, mr, mo, status = "ok")
      path_list[[length(path_list) + 1L]] <-
        cbind(marker_id = mk, path, stringsAsFactors = FALSE)
    } else {
      for (tol in tolerances) {
        path_list[[length(path_list) + 1L]] <- data.frame(
          marker_id = mk, tolerance = tol,
          minimum_set_size = NA_integer_, n_minimum_sets = 0L,
          n_irreducible_modules = 0L,
          status = if (subset_mode != "all" && locus_status == "ok") {
            "not_applicable"   # extraction only defined for mode "all"
          } else {
            locus_status
          },
          stringsAsFactors = FALSE
        )
      }
    }
  }

  bind_or_empty <- function(lst, empty) {
    if (length(lst) == 0L) return(empty)
    out <- do.call(rbind, lst)
    rownames(out) <- NULL
    out
  }
  subset_table <- bind_or_empty(subset_rows, empty_subset)
  min_reps <- bind_or_empty(reps_list, empty_reps)
  modules <- bind_or_empty(mods_list, empty_mods)
  tol_path <- bind_or_empty(path_list, empty_path)

  n_loci <- length(markers)
  code <- if (n_loci == 0L ||
              (nrow(subset_table) == 0L && n_too_many == 0L)) {
    "empty_selection"
  } else {
    "ok"
  }
  list(
    subset_table = if (return_subset_table) subset_table else NULL,
    minimum_representative_sets = min_reps,
    irreducible_modules = modules,
    tolerance_path = tol_path,
    status = .new_status(
      ok = TRUE, code = code,
      message = if (code == "empty_selection") {
        "No loci with candidate traits to analyse."
      } else {
        ""
      }
    ),
    diagnostics = list(
      n_loci = n_loci,
      n_subset_rows = nrow(subset_table),
      n_empty_candidates = n_empty_candidates,
      n_single_candidate = n_single_candidate,
      n_too_many_traits = n_too_many,
      max_decomposition_error = max_decomp,
      max_monotonicity_violation = max_mono,
      n_rank_deficient = n_rank_def_rows,
      n_unstable = n_unstable_loci
    )
  )
}

#' Convert an effects input to a per-marker named beta list
#'
#' @param effects An [estimate_mt_effects()] result or an
#'   `effects_long` data.frame.
#' @return Named list: marker id -> named numeric effect vector.
#' @keywords internal
.effects_to_beta_list <- function(effects) {
  if (is.list(effects) && !is.data.frame(effects)) {
    if (!is.null(effects$status) && isFALSE(effects$status$ok)) {
      .stop_invalid_input(
        "effects result has status$ok == FALSE; cannot decompose."
      )
    }
    effects <- effects$effects_long
  }
  if (!is.data.frame(effects) ||
      !all(c("marker_id", "trait", "beta") %in% names(effects))) {
    .stop_invalid_input(
      "effects must provide an effects_long data.frame with marker_id, trait and beta."
    )
  }
  if (anyNA(effects$beta)) {
    .stop_invalid_input("effects contains NA beta values.")
  }
  out <- split(
    seq_len(nrow(effects)),
    factor(effects$marker_id, levels = unique(effects$marker_id))
  )
  lapply(out, function(idx) {
    b <- effects$beta[idx]
    names(b) <- as.character(effects$trait[idx])
    b
  })
}

#' Extract candidate trait sets from an attribution input
#'
#' @param attribution An [attribute_traits()] result or a named list
#'   of trait vectors.
#' @param required Logical; when `FALSE` (oracle mode) `NULL` is
#'   accepted and an empty list returned.
#' @return Named list: marker id -> character vector of traits.
#' @keywords internal
.attribution_to_candidate_sets <- function(attribution, required = TRUE) {
  if (is.null(attribution)) {
    if (required) {
      .stop_invalid_input(
        "attribution is required when restrict_to_candidates = TRUE."
      )
    }
    return(list())
  }
  sets <- attribution
  if (is.list(attribution) && !is.data.frame(attribution) &&
      !is.null(names(attribution))) {
    if (!is.null(attribution$status) && isFALSE(attribution$status$ok)) {
      .stop_invalid_input(
        "attribution result has status$ok == FALSE; cannot decompose."
      )
    }
    if (!is.null(attribution$candidate_sets)) {
      sets <- attribution$candidate_sets
    } else if (!is.null(attribution$A)) {
      sets <- attribution$A   # pre-v1.0 field name
    }
  }
  if (!is.list(sets) || is.null(names(sets)) ||
      !all(vapply(sets, function(x) {
        is.character(x) || is.null(x)
      }, logical(1)))) {
    .stop_invalid_input(
      "attribution must provide candidate trait sets (a named list of character vectors)."
    )
  }
  lapply(sets, function(x) if (is.null(x)) character() else x)
}
