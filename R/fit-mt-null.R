#' Fit the multi-trait REML null model
#'
#' Fits the multi-trait null model
#' \eqn{y = X_0 b + u + \varepsilon} with covariance
#' \eqn{V = \Sigma_G \otimes K + \Sigma_E \otimes I_n} by restricted maximum
#' likelihood, and returns the rotation objects and global conditional
#' contrasts reused by the whole pipeline.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param Y Numeric `n x m` phenotype matrix (rows are individuals).
#' @param W Optional numeric `n x q` fixed-effect design matrix; when `NULL`
#'   an intercept-only design is used.
#' @param K Numeric `n x n` genomic relationship matrix, scaled to average
#'   diagonal 1.
#' @param add_intercept Logical; add an intercept column to `W` when it is
#'   not already present.
#' @param optimizer Optimizer used for the restricted log-likelihood;
#'   `"BFGS"` or `"nlminb"`.
#' @param n_starts Integer; number of random starting points.
#' @param h2_start Numeric in (0, 1); heritability used to split the sample
#'   phenotypic covariance into starting values.
#' @param start_values Optional explicit starting values (log-Cholesky
#'   parameterisation).
#' @param eig_tol Tolerance for the eigen decomposition of `K`.
#' @param inverse_tol Relative tolerance forwarded to `.safe_inverse()`.
#' @param chol_floor Lower bound for log-Cholesky diagonal parameters.
#' @param control List of optimizer control settings (e.g. `maxit`,
#'   `reltol`).
#' @param compute_gamma Logical; compute the global conditional projection
#'   coefficients from \eqn{\widehat\Sigma_P = \widehat\Sigma_G +
#'   \widehat\Sigma_E}.
#' @param return_rotation Logical; return the rotation object (eigenvectors,
#'   eigenvalues, rotated data and related components).
#' @param verbose Logical; print optimisation progress.
#'
#' @return An object of class `"condped_mt_null"`: a list with components
#'   `Sigma_G`, `Sigma_E`, `Sigma_P`, `fixed_effects`, `gamma`, `contrasts`,
#'   `conditional_variance`, `logLik`, `npar`, `optimizer`, `convergence`,
#'   `rotation`, `diagnostics` and `status`, as specified in the interface
#'   contract.
#' @export
fit_mt_null <- function(
  Y,
  W = NULL,
  K,
  add_intercept = TRUE,
  optimizer = c("BFGS", "nlminb"),
  n_starts = 5L,
  h2_start = 0.50,
  start_values = NULL,
  eig_tol = 1e-8,
  inverse_tol = sqrt(.Machine$double.eps),
  chol_floor = 1e-8,
  control = list(maxit = 500L, reltol = 1e-8),
  compute_gamma = TRUE,
  return_rotation = TRUE,
  verbose = FALSE
) {
  .not_implemented("fit_mt_null")
}
