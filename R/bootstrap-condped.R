#' Parametric bootstrap for the CondPED pipeline
#'
#' Propagates the parameter uncertainty of the fitted null model, effects,
#' conditional deviations and PVE through parametric (or residual)
#' bootstrap. With `reselect = FALSE` the reported loci are held fixed; with
#' `reselect = TRUE` the full pipeline is rerun on each replicate to assess
#' selection stability.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param object A fitted CondPED object (e.g. from [condped()]) or a
#'   `"condped_mt_null"` fit.
#' @param Y Optional numeric `n x m` phenotype matrix.
#' @param W Optional numeric `n x q` fixed-effect design matrix.
#' @param G Optional numeric `n x p` genotype dosage matrix.
#' @param loci Optional vector of loci to bootstrap; required when `object`
#'   does not carry a selection.
#' @param B Integer; number of bootstrap replicates.
#' @param type Bootstrap type: `"parametric"` or `"residual"`.
#' @param refit_null Logical; re-estimate the variance components and the
#'   conditional projection coefficients in every replicate.
#' @param reselect Logical; rerun the full pipeline (including selection) in
#'   every replicate.
#' @param interval Confidence interval type: `"percentile"` or `"bca"`.
#' @param seed Random seed.
#' @param workers Integer; number of parallel workers.
#' @param min_valid Minimum fraction of valid replicates; below this
#'   threshold the result carries status code `"bootstrap_failed"`.
#'
#' @return A list with components `replicates`, `intervals` (data.frame),
#'   optional `selection_frequency`, `status` and `diagnostics` (including
#'   `B`, `n_valid`, `refit_null`, `reselect`), as specified in the
#'   interface contract.
#' @export
bootstrap_condped <- function(
  object,
  Y = NULL,
  W = NULL,
  G = NULL,
  loci = NULL,
  B = 200L,
  type = c("parametric", "residual"),
  refit_null = TRUE,
  reselect = FALSE,
  interval = c("percentile", "bca"),
  seed = 1L,
  workers = 1L,
  min_valid = 0.80
) {
  .not_implemented("bootstrap_condped")
}
