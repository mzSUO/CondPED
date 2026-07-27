#' End-to-end CondPED analysis
#'
#' Runs the full CondPED pipeline: input checks, construction or reading of
#' the genomic relationship matrix, global multi-trait null model, omnibus
#' scan, effect estimation, trait attribution, conditional-deviation
#' decomposition, PVE estimation, and optional cross-fitting and bootstrap.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param Y Numeric `n x m` phenotype matrix (rows are individuals).
#' @param W Optional numeric `n x q` fixed-effect design matrix; when `NULL`
#'   an intercept-only design is used.
#' @param G Numeric `n x p` genotype dosage matrix coded 0/1/2.
#' @param K Optional numeric `n x n` genomic relationship matrix; when
#'   `NULL` it is built from `G`.
#' @param chromosome Optional vector of chromosome labels per marker, used
#'   for leave-one-chromosome-out scans. The global conditional projection
#'   coefficients are never redefined by these scans.
#' @param alpha_omnibus First-layer omnibus significance level.
#' @param omnibus_adjust Adjustment method for the omnibus p-values.
#' @param attribution_mode Second-layer correction mode passed to
#'   [attribute_traits()].
#' @param q_target Target FDR level for the attribution step.
#' @param conditional_adjust Adjustment method for the conditional tests.
#' @param pve_correction PVE bias correction passed to
#'   [estimate_locus_pve()].
#' @param crossfit Logical; run [crossfit_effect_pve()].
#' @param bootstrap Logical; run [bootstrap_condped()].
#' @param control List of additional control settings passed to sub-modules.
#' @param seed Random seed.
#'
#' @return An object of class `"condped_fit"`: a list with components
#'   `null_global`, `null_loco`, `omnibus`, `effects`, `attribution`,
#'   `conditional`, `pve`, `crossfit`, `bootstrap`, `output_long`,
#'   `settings`, `status` and `diagnostics`, as specified in the interface
#'   contract.
#' @export
condped <- function(
  Y,
  W = NULL,
  G,
  K = NULL,
  chromosome = NULL,
  alpha_omnibus = 0.05,
  omnibus_adjust = "BH",
  attribution_mode = "bb_fdr",
  q_target = 0.05,
  conditional_adjust = "holm",
  pve_correction = "both",
  crossfit = FALSE,
  bootstrap = FALSE,
  control = list(),
  seed = 1L
) {
  .not_implemented("condped")
}
