# Analysis-mode mapping for the formal simulation error-decomposition
# modes (frozen; "restricted" no longer exists).

#' Map an analysis mode to condped() settings
#'
#' Frozen mapping:
#' \describe{
#'   \item{full}{`signal_mode = "resolve"`,
#'     `candidate_mode = "holm_fwer"`.}
#'   \item{signal_oracle}{`signal_mode = "predefined"`,
#'     `candidate_mode = "holm_fwer"` (true signals, estimated traits).}
#'   \item{signal_trait_oracle}{`signal_mode = "predefined"`,
#'     `candidate_mode = "predefined"` (true signals and true
#'     candidate sets; `all_traits` never substitutes for a trait
#'     oracle).}
#' }
#'
#' @param mode One of `"full"`, `"signal_oracle"`,
#'   `"signal_trait_oracle"`.
#' @return A list with `signal_mode` and `candidate_mode`.
#' @keywords internal
.analysis_mode_mapping <- function(mode = c("full", "signal_oracle",
                                            "signal_trait_oracle")) {
  mode <- match.arg(mode)
  switch(mode,
    full = list(signal_mode = "resolve",
                candidate_mode = "holm_fwer"),
    signal_oracle = list(signal_mode = "predefined",
                         candidate_mode = "holm_fwer"),
    signal_trait_oracle = list(signal_mode = "predefined",
                               candidate_mode = "predefined")
  )
}
