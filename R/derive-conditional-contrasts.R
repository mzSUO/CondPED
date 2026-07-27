#' Derive conditional projection contrasts
#'
#' Computes the conditional projection coefficients
#' \eqn{\gamma_{i,-i} = \Sigma_{P,-i,-i}^{-1}\Sigma_{P,-i,i}} for every trait
#' from the (fitted or true) phenotypic covariance matrix, together with the
#' contrast matrix `C` and conditional variances. Orthogonality is verified
#' on the model-covariance scale, never on empirical sample covariances.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param Sigma_P Numeric `m x m` phenotypic covariance matrix.
#' @param trait_names Character vector of trait names; defaults to
#'   `colnames(Sigma_P)`.
#' @param inverse_tol Relative tolerance forwarded to `.safe_inverse()`.
#' @param ridge Non-negative ridge added to the diagonal when submatrices
#'   are ill-conditioned; the used value is recorded.
#' @param standardize Logical; standardise contrasts to unit conditional
#'   variance.
#'
#' @return A list with components `gamma` (per-trait coefficient lists),
#'   `C` (an `m x m` matrix whose i-th column is the contrast vector
#'   \eqn{c_i}), `conditional_variance` (length `m`), `diagnostics` (a
#'   data.frame with columns `trait`, `rank`, `condition_number`,
#'   `ridge_used`, `orthogonality_error`) and `status`, as specified in the
#'   interface contract.
#' @export
derive_conditional_contrasts <- function(
  Sigma_P,
  trait_names = colnames(Sigma_P),
  inverse_tol = sqrt(.Machine$double.eps),
  ridge = 0,
  standardize = FALSE
) {
  .not_implemented("derive_conditional_contrasts")
}
