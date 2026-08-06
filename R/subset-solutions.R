# Subset enumeration and set-valued solutions for CondPED v1.0 (Stage 2).
#
# These helpers operate on top of the unique loss path
# `.compute_subset_loss()` (R/subset-loss.R); the loss formulas are
# never duplicated here. Extraction functions read exactly one
# per-locus subset table and never recompute losses.

#' Enumerate candidate representing subsets
#'
#' @param A Normalised candidate trait set (character vector).
#' @param mode One of `"all"`, `"singleton"`, `"custom"`.
#' @param custom_sets Optional list of trait sets for `mode = "custom"`.
#'
#' @return A list with components `sets` (list of normalised character
#'   vectors, unique by stable key) and `keys` (their string keys).
#'   `"all"` returns all `2^|A|` subsets ordered by size then trait
#'   order; `"singleton"` returns the empty set, the full set, every
#'   singleton and every singleton complement (duplicates removed, so
#'   `|A| = 2` yields 4 sets); `"custom"` returns the normalised
#'   `custom_sets` plus the empty and full sets.
#' @keywords internal
.enumerate_subsets <- function(A,
                               mode = c("all", "singleton", "custom"),
                               custom_sets = NULL) {
  mode <- match.arg(mode)
  k <- length(A)
  sets <- switch(mode,
    all = {
      out <- list(character())
      if (k > 0L) {
        for (size in seq_len(k)) {
          out <- c(out, utils::combn(A, size, simplify = FALSE))
        }
      }
      out
    },
    singleton = {
      out <- list(character(), A)
      for (i in seq_len(k)) {
        out <- c(out, list(A[i]), list(setdiff(A, A[i])))
      }
      out
    },
    custom = {
      if (is.null(custom_sets) || length(custom_sets) == 0L) {
        .stop_invalid_input(
          "custom_sets must be a non-empty list when subset_mode = \"custom\"."
        )
      }
      out <- lapply(custom_sets, function(s) .normalize_trait_set(s, A, "custom_sets"))
      c(list(character()), out, list(A))
    }
  )
  # Deduplicate by stable key, preserving first occurrence.
  keys <- vapply(sets, .trait_set_key, character(1))
  keep <- !duplicated(keys)
  list(sets = sets[keep], keys = keys[keep])
}

#' Check the quadratic-form decomposition identity
#'
#' Relative error |full - subset - residual| / max(full, qform_tol).
#' This is the same identity `.compute_subset_loss()` reports per
#' subset; the helper exists so validation code and later stages share
#' one definition.
#'
#' @param full_qform,subset_qform,residual_qform Numeric scalars.
#' @param qform_tol Positive energy threshold.
#' @return Numeric scalar relative error.
#' @keywords internal
.check_qform_decomposition <- function(full_qform, subset_qform,
                                       residual_qform, qform_tol = 1e-12) {
  abs(full_qform - subset_qform - residual_qform) /
    max(full_qform, qform_tol)
}

#' Full-table monotonicity check for one locus
#'
#' For every pair S subset of T in a complete per-locus subset table,
#' verifies rho(T) <= rho(S) + monotonicity_tol. Only meaningful for
#' the complete enumeration (`subset_mode = "all"`). Rows with
#' non-usable (NA) losses are skipped.
#'
#' @param tab Per-locus subset table with columns `representing_set`
#'   (list-column), `set_size` and `representation_loss`.
#' @param monotonicity_tol Non-negative tolerance.
#' @return A list with `max_violation` (numeric, 0 when none) and
#'   `n_violations` (integer).
#' @keywords internal
.check_loss_monotonicity <- function(tab, monotonicity_tol = 1e-10) {
  usable <- !is.na(tab$representation_loss)
  tab <- tab[usable, , drop = FALSE]
  max_violation <- 0
  n_violations <- 0L
  n <- nrow(tab)
  if (n < 2L) {
    return(list(max_violation = 0, n_violations = 0L))
  }
  for (i in seq_len(n - 1L)) {
    Si <- tab$representing_set[[i]]
    li <- tab$representation_loss[i]
    for (j in (i + 1L):n) {
      Sj <- tab$representing_set[[j]]
      lj <- tab$representation_loss[j]
      if (all(Si %in% Sj)) {
        viol <- lj - li          # rho(T) - rho(S) must be <= tol
        if (viol > monotonicity_tol) n_violations <- n_violations + 1L
        if (viol > max_violation) max_violation <- viol
      }
      if (all(Sj %in% Si)) {
        viol <- li - lj
        if (viol > monotonicity_tol) n_violations <- n_violations + 1L
        if (viol > max_violation) max_violation <- viol
      }
    }
  }
  list(max_violation = max_violation, n_violations = n_violations)
}

#' Extract all minimum representative trait sets from one subset table
#'
#' Feasible sets satisfy `representation_loss <= tolerance`; the
#' minimum feasible size wins and ALL tied sets are returned (no
#' tie-breaking by loss or lexicographic order).
#'
#' @param tab Per-locus subset table (see `.check_loss_monotonicity`).
#' @param tolerance Numeric scalar in [0, 1).
#'
#' @return A data.frame with columns `tolerance`, `solution_id`,
#'   `trait_set` (list-column), `trait_key`, `set_size`,
#'   `representation_loss`, `n_tied_solutions`, `status`; zero rows
#'   when no usable feasible set exists.
#' @keywords internal
.extract_minimum_representative_sets <- function(tab, tolerance) {
  empty <- data.frame(
    tolerance = numeric(), solution_id = integer(),
    trait_set = I(list()), trait_key = character(),
    set_size = integer(), representation_loss = numeric(),
    n_tied_solutions = integer(), status = character()
  )
  usable <- !is.na(tab$representation_loss)
  # Feasibility includes the exact boundary: constructions with
  # target_loss equal to a sensitivity tolerance are feasible by
  # definition; the 1e-10 slack absorbs floating-point noise only and
  # is far below any truth_margin.
  feasible <- usable &
    tab$representation_loss <= tolerance + 1e-10
  if (!any(feasible)) return(empty)
  min_size <- min(tab$set_size[feasible])
  win <- which(feasible & tab$set_size == min_size)
  data.frame(
    tolerance = tolerance,
    solution_id = seq_along(win),
    trait_set = I(tab$representing_set[win]),
    trait_key = tab$representing_key[win],
    set_size = tab$set_size[win],
    representation_loss = tab$representation_loss[win],
    n_tied_solutions = length(win),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

#' Extract all irreducible trait modules from one subset table
#'
#' A non-empty set T is an irreducible module when removing it makes
#' representation infeasible (`rho(A \ T) > tolerance`) while no
#' proper subset of T has that property (inclusion-minimality). All
#' modules are returned, including singleton, joint and multiple
#' mutually non-containing modules. Requires a complete subset table
#' (`subset_mode = "all"`): the complement loss of every candidate
#' module must be present. NA (unusable) losses never qualify a module.
#'
#' @param tab Complete per-locus subset table.
#' @param tolerance Numeric scalar in [0, 1).
#'
#' @return A data.frame with columns `tolerance`, `module_id`,
#'   `trait_set` (list-column), `trait_key`, `module_size`,
#'   `complement_set` (list-column), `complement_key`,
#'   `complement_loss`, `status`.
#' @keywords internal
.extract_irreducible_modules <- function(tab, tolerance) {
  empty <- data.frame(
    tolerance = numeric(), module_id = integer(),
    trait_set = I(list()), trait_key = character(),
    module_size = integer(),
    complement_set = I(list()), complement_key = character(),
    complement_loss = numeric(), status = character()
  )
  k <- max(tab$set_size)
  if (k == 0L) return(empty)
  A <- tab$representing_set[[which(tab$set_size == k)[1L]]]
  loss_by_key <- tab$representation_loss
  names(loss_by_key) <- tab$representing_key

  # Non-empty candidate modules, smallest first (enables minimality
  # pruning by scanning known modules).
  candidates <- unlist(
    lapply(seq_len(k), function(sz) utils::combn(A, sz, simplify = FALSE)),
    recursive = FALSE
  )
  qualifies <- function(T) {
    comp_key <- .trait_set_key(setdiff(A, T))
    comp_loss <- loss_by_key[[comp_key]]
    # mirror of the feasibility rule: strictly above the tolerance,
    # with the same 1e-10 floating-point slack
    !is.null(comp_loss) && !is.na(comp_loss) &&
      comp_loss > tolerance + 1e-10
  }
  modules <- list()
  for (T in candidates) {
    # inclusion-minimality: no already-found module may be a proper
    # subset of T (candidates are ordered by size, so any qualifying
    # proper subset is already in `modules`)
    if (any(vapply(modules, function(M) all(M %in% T), logical(1)))) next
    if (qualifies(T)) modules <- c(modules, list(T))
  }
  if (length(modules) == 0L) return(empty)
  data.frame(
    tolerance = tolerance,
    module_id = seq_along(modules),
    trait_set = I(modules),
    trait_key = vapply(modules, .trait_set_key, character(1)),
    module_size = lengths(modules),
    complement_set = I(lapply(modules, function(T) setdiff(A, T))),
    complement_key = vapply(modules,
                            function(T) .trait_set_key(setdiff(A, T)),
                            character(1)),
    complement_loss = vapply(
      modules,
      function(T) unname(loss_by_key[[.trait_set_key(setdiff(A, T))]]),
      numeric(1)
    ),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

#' Build the tolerance path for one locus
#'
#' One row per tolerance: size of the minimum representative sets, how
#' many tied solutions, and how many irreducible modules.
#'
#' @param tolerances Numeric vector of tolerances.
#' @param min_reps,modules Extraction results per tolerance (lists
#'   keyed by tolerance) or `NULL` when extraction was not performed.
#' @param status Locus-level status string.
#' @return A data.frame with columns `tolerance`, `minimum_set_size`,
#'   `n_minimum_sets`, `n_irreducible_modules`, `status`.
#' @keywords internal
.build_tolerance_path <- function(tolerances, min_reps, modules,
                                  status = "ok") {
  rows <- lapply(tolerances, function(tol) {
    key <- format(tol)
    mr <- min_reps[[key]]
    mo <- modules[[key]]
    data.frame(
      tolerance = tol,
      minimum_set_size = if (!is.null(mr) && nrow(mr) > 0L) {
        mr$set_size[1L]
      } else {
        NA_integer_
      },
      n_minimum_sets = if (!is.null(mr)) nrow(mr) else 0L,
      n_irreducible_modules = if (!is.null(mo)) nrow(mo) else 0L,
      status = status,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
