# Trait-set normalisation utilities for CondPED v1.0.
# Frozen key rules (interface contract section 4.2):
#   empty set  : <empty>
#   non-empty  : traitA|traitB|traitC
# Sets are always stored as character trait names, deduplicated and
# ordered by the reference trait_names order; comparisons never depend
# on user input order.

#' Normalise a trait set against reference trait names
#'
#' Converts a user-supplied trait set (character names or integer
#' positions) into a deduplicated character vector ordered by
#' `trait_names`. The empty set normalises to `character(0)`.
#'
#' @param set `NULL`, a character vector of trait names, or an integer
#'   vector of positions in `trait_names`.
#' @param trait_names Character vector of unique reference trait names.
#' @param arg Name of the calling argument, used in error messages.
#'
#' @return Character vector of unique trait names in `trait_names`
#'   order.
#' @keywords internal
.normalize_trait_set <- function(set, trait_names, arg = "set") {
  if (!is.character(trait_names) || length(trait_names) == 0L ||
      anyNA(trait_names) || anyDuplicated(trait_names)) {
    .stop_invalid_input(
      "trait_names must be a non-empty character vector of unique names."
    )
  }
  if (is.null(set) || length(set) == 0L) return(character())
  if (is.numeric(set)) {
    if (any(!is.finite(set)) || any(set != round(set)) ||
        any(set < 1) || any(set > length(trait_names))) {
      .stop_invalid_input(
        "%s contains invalid trait positions (must be integers between 1 and %d).",
        arg, length(trait_names)
      )
    }
    set <- trait_names[as.integer(set)]
  }
  if (!is.character(set) || anyNA(set)) {
    .stop_invalid_input(
      "%s must be a character vector of trait names or integer positions.",
      arg
    )
  }
  unknown <- setdiff(set, trait_names)
  if (length(unknown) > 0L) {
    .stop_invalid_input(
      "%s contains unknown trait(s): %s.",
      arg, paste(utils::head(unknown, 3L), collapse = ", ")
    )
  }
  # Deduplicate and order by the reference trait order.
  trait_names[trait_names %in% set]
}

#' Stable string key for a normalised trait set
#'
#' Frozen key rule: the empty set maps to `"<empty>"`, any non-empty
#' set to its members joined by `"|"`. The input must already be
#' normalised (see [`.normalize_trait_set()`]) so that keys never
#' depend on input order.
#'
#' @param set Character vector of trait names (already normalised).
#'
#' @return A single string key.
#' @keywords internal
.trait_set_key <- function(set) {
  if (length(set) == 0L) return("<empty>")
  paste(set, collapse = "|")
}
