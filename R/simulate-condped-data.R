#' Simulate multi-trait data with known ground truth
#'
#' Generates simulated data with a known multi-trait effect matrix, true
#' associated-trait sets, true conditional-deviation sets, theoretical PVE
#' and complete covariance ground truth.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param n Positive integer; number of individuals. Main simulations use
#'   500, 750, 1000 or 1250.
#' @param m Integer; number of traits (default 4).
#' @param p Integer; total number of background markers.
#' @param maf_range Numeric vector of length 2 giving the range of minor
#'   allele frequencies.
#' @param h2 Numeric vector of length `m`; target SNP heritabilities.
#' @param R_G Optional `m x m` genetic correlation matrix; `NULL` uses an
#'   exchangeable correlation of 0.4.
#' @param R_E Optional `m x m` residual correlation matrix; `NULL` uses an
#'   exchangeable correlation of 0.4.
#' @param architecture Genetic architecture of the focal locus/loci; one of
#'   `"null"`, `"single_trait"`, `"shared_same"`, `"shared_opposite"`,
#'   `"dense"`, `"covariance_aligned"`, `"conditional_deviation"`,
#'   `"projection_induced"`.
#' @param locus_pve Scalar or per-active-trait vector of target locus PVE.
#' @param delta True deviation effect in the conditional-deviation
#'   architecture.
#' @param target_trait Integer index of the target trait.
#' @param n_groups Number of families/subpopulations used when
#'   `structured = TRUE`.
#' @param fst Mild between-group allele frequency differentiation.
#' @param structured Logical; when `TRUE` individuals carry family/subgroup
#'   structure so that the genomic relationship matrix does not degenerate
#'   to the identity.
#' @param n_qtl Integer; number of focal QTL.
#' @param max_attempts Maximum number of regeneration attempts when the
#'   background genetic covariance is not positive semi-definite.
#' @param seed Optional random seed.
#' @param return_latent Logical; whether to return latent genetic and
#'   residual values.
#'
#' @return A list with components `Y`, `W`, `G`, `K_bg`, `qtl_index`,
#'   `truth` (effect matrices, true associated-trait sets `A`, true
#'   deviation sets `D`, `eta`, PVE matrices and all covariance ground
#'   truth), `latent`, `generator`, `status` and `diagnostics`, as specified
#'   in the interface contract.
#' @export
simulate_condped_data <- function(
  n = 1000L,
  m = 4L,
  p = 2000L,
  maf_range = c(0.05, 0.50),
  h2 = rep(0.50, m),
  R_G = NULL,
  R_E = NULL,
  architecture = c(
    "null",
    "single_trait",
    "shared_same",
    "shared_opposite",
    "dense",
    "covariance_aligned",
    "conditional_deviation",
    "projection_induced"
  ),
  locus_pve = 0.01,
  delta = 0,
  target_trait = 1L,
  n_groups = 8L,
  fst = 0.05,
  structured = TRUE,
  n_qtl = 1L,
  max_attempts = 100L,
  seed = NULL,
  return_latent = TRUE
) {
  .not_implemented("simulate_condped_data")
}
