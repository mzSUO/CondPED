#' Multi-trait marker effect estimation by GLS
#'
#' Estimates the multi-trait effects of selected loci by generalised least
#' squares, \eqn{\widehat\beta = J^+ U}, reusing the rotation object of the
#' null fit.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param null_fit An object of class `"condped_mt_null"` returned by
#'   [fit_mt_null()].
#' @param G Numeric `n x p` genotype dosage matrix coded 0/1/2.
#' @param loci Vector of locus indices or marker identifiers to estimate.
#' @param marker_ids Character vector of marker identifiers; defaults to
#'   `colnames(G)`.
#' @param rank_tol Relative tolerance used for the generalised inverse of
#'   \eqn{J}.
#' @param return_covariance Logical; return the per-locus `m x m` effect
#'   covariance matrices.
#'
#' @return A list with components `effects_long` (a data.frame with columns
#'   `marker_id`, `trait`, `beta`, `se`, `z`, `p_value`), `beta` (an
#'   n_loci x m matrix), `covariance` (an m x m x n_loci array), `status`
#'   and `diagnostics`, as specified in the interface contract.
#' @export
estimate_mt_effects <- function(
  null_fit,
  G,
  loci,
  marker_ids = colnames(G),
  rank_tol = sqrt(.Machine$double.eps),
  return_covariance = TRUE
) {
  .not_implemented("estimate_mt_effects")
}
