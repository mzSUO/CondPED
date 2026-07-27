#' Estimate locus-level PVE
#'
#' Estimates the proportion of phenotypic variance explained by each locus
#' (and, optionally, the conditional PVE per trait), with an optional noise
#' correction based on \eqn{\widehat\beta^2 - \widehat{s}^2}.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param effects Effect estimates returned by [estimate_mt_effects()].
#' @param genotype_variance Numeric vector of per-locus genotype variances.
#' @param Sigma_P Numeric `m x m` phenotypic covariance matrix used as the
#'   PVE denominator.
#' @param conditional Optional conditional-decomposition result returned by
#'   [decompose_conditional_effects()].
#' @param correction Bias correction to report: `"raw"`,
#'   `"noise_corrected"` or `"both"`.
#' @param truncate_zero Logical; truncate corrected PVE at zero. The
#'   truncated version does not claim strict unbiasedness.
#' @param gate Logical; when `TRUE` conditional PVE for non-attributed
#'   traits is not part of the formal output.
#' @param attribution Optional attribution result returned by
#'   [attribute_traits()], required for gating.
#'
#' @return A list with component `pve_table` (a data.frame with columns
#'   `marker_id`, `trait`, `pve_raw`, `pve_untruncated_bc`, `pve_bc`,
#'   `conditional_pve_raw`, `conditional_pve_untruncated_bc`,
#'   `conditional_pve_bc`, `reportable`), plus `status` and `diagnostics`,
#'   as specified in the interface contract.
#' @export
estimate_locus_pve <- function(
  effects,
  genotype_variance,
  Sigma_P,
  conditional = NULL,
  correction = c("raw", "noise_corrected", "both"),
  truncate_zero = TRUE,
  gate = TRUE,
  attribution = NULL
) {
  .not_implemented("estimate_locus_pve")
}
