#' Cross-fitted effect and PVE estimation
#'
#' Splits the sample into discovery and estimation folds: the discovery fold
#' runs the omnibus scan and trait attribution, the estimation fold
#' re-estimates effects and PVE only at the loci selected in the discovery
#' fold, mitigating selection bias (winner's curse).
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param Y Numeric `n x m` phenotype matrix.
#' @param W Numeric `n x q` fixed-effect design matrix (must explicitly
#'   contain an intercept or be amendable to one).
#' @param G Numeric `n x p` genotype dosage matrix coded 0/1/2.
#' @param fold_id Optional integer vector of fold assignments; `NULL`
#'   generates folds at random.
#' @param n_folds Integer; number of folds (default 2).
#' @param n_repeats Integer; number of repeated random splits.
#' @param analysis_args List of further arguments passed to the analysis
#'   functions.
#' @param gamma_mode `"global_fixed"` estimates the conditional projection
#'   coefficients once on the full data and keeps them fixed across folds;
#'   `"fold_specific"` re-estimates them within each training fold and
#'   reports the two directions separately (never merged).
#' @param seed Random seed.
#' @param workers Integer; number of parallel workers.
#' @param keep_logs Logical; keep per-direction audit logs.
#'
#' @return A list with components `fold_results`, `aggregate` (data.frame),
#'   `audit_log` (a data.frame with columns `repeat`, `direction`,
#'   `discovery_ids_hash`, `estimation_ids_hash`, `gamma_mode`,
#'   `selected_loci`), `status` and `diagnostics`, as specified in the
#'   interface contract.
#' @export
crossfit_effect_pve <- function(
  Y,
  W,
  G,
  fold_id = NULL,
  n_folds = 2L,
  n_repeats = 10L,
  analysis_args = list(),
  gamma_mode = c("global_fixed", "fold_specific"),
  seed = 1L,
  workers = 1L,
  keep_logs = TRUE
) {
  .not_implemented("crossfit_effect_pve")
}
