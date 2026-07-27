#' Multi-trait omnibus score scan
#'
#' Runs the joint multi-trait score test for every marker, strictly following
#' equations (9)-(12) of the Methods: \eqn{U_l = Z_l^\top P y},
#' \eqn{J_l = Z_l^\top P Z_l}, \eqn{Q_l = U_l^\top J_l^+ U_l}. The full
#' \eqn{Z_l} and \eqn{P} are never constructed explicitly; block computations
#' use the rotation object of the null fit.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param null_fit An object of class `"condped_mt_null"` returned by
#'   [fit_mt_null()].
#' @param G Numeric `n x p` genotype dosage matrix coded 0/1/2.
#' @param marker_ids Character vector of marker identifiers; defaults to
#'   `colnames(G)`.
#' @param maf_min Minimum minor allele frequency; markers below the threshold
#'   are filtered.
#' @param chunk_size Integer; number of markers processed per block.
#' @param rank_tol Relative tolerance used for the rank and generalised
#'   inverse of \eqn{J_l}.
#' @param return_score Logical; also return the score vectors.
#' @param return_effects Logical; also return marker effect estimates.
#' @param bootstrap_p Logical; compute bootstrap p-values.
#' @param n_boot Integer; number of bootstrap replicates when
#'   `bootstrap_p = TRUE`.
#' @param seed Optional random seed.
#'
#' @return A list with components `omnibus` (a data.frame with columns
#'   `marker_id`, `maf`, `n_eff`, `Q`, `df`, `p_value`, `rank_J`,
#'   `condition_J`, `status`), optional `score` and `effects`, plus `status`
#'   and `diagnostics`, as specified in the interface contract.
#' @export
scan_mt_omnibus <- function(
  null_fit,
  G,
  marker_ids = colnames(G),
  maf_min = 0.05,
  chunk_size = 2000L,
  rank_tol = sqrt(.Machine$double.eps),
  return_score = FALSE,
  return_effects = FALSE,
  bootstrap_p = FALSE,
  n_boot = 0L,
  seed = NULL
) {
  .not_implemented("scan_mt_omnibus")
}
