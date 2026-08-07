#' Representation basis: fixed reference covariance for subset analysis
#'
#' The public name is kept for package compatibility. In v1.0 this
#' function does NOT pre-compute trait contrasts; it validates and
#' freezes the fixed reference phenotypic covariance used by every
#' subset loss evaluation, together with the numerical settings shared
#' by all downstream calls.
#'
#' @param Sigma_P Numeric symmetric `m x m` reference covariance matrix
#'   (typically `null_fit$Sigma_P_ref`).
#' @param trait_names Unique trait names; defaults to
#'   `colnames(Sigma_P)`.
#' @param inverse_tol Eigenvalue tolerance used for the rank, condition
#'   number and pseudo-inverse diagnostics.
#' @param ridge Non-negative scalar added to the diagonal only when the
#'   user explicitly requests it; the default `0` never inflates the
#'   diagonal automatically. The value used is recorded in `basis`.
#' @param standardize Logical; frozen default `FALSE` because the
#'   conditional representation loss is invariant to nonzero diagonal
#'   rescalings of the traits. When `TRUE`, the covariance is converted
#'   to a correlation matrix.
#'
#' @return An object of class `"condped_representation_basis"`: a list
#'   with components `Sigma_P_ref` (the covariance actually used,
#'   after optional ridge/standardization), `trait_names`, `basis`
#'   (list with `inverse_tol`, `ridge`, `standardize`), `status` and
#'   `diagnostics` (list with `rank`, `min_eigenvalue`,
#'   `condition_number`, `used_pseudoinverse`).
#' @export
derive_conditional_contrasts <- function(
  Sigma_P,
  trait_names = colnames(Sigma_P),
  inverse_tol = sqrt(.Machine$double.eps),
  ridge = 0,
  standardize = FALSE
) {
  if (!is.numeric(Sigma_P) || !is.matrix(Sigma_P) ||
      nrow(Sigma_P) != ncol(Sigma_P) || nrow(Sigma_P) == 0L ||
      any(!is.finite(Sigma_P))) {
    .stop_invalid_input(
      "Sigma_P must be a non-empty finite numeric square matrix."
    )
  }
  m <- nrow(Sigma_P)
  if (max(abs(Sigma_P - t(Sigma_P))) > 1e-8 * max(1, max(abs(Sigma_P)))) {
    .stop_invalid_input("Sigma_P must be symmetric.")
  }
  Sigma_P <- (Sigma_P + t(Sigma_P)) / 2
  if (is.null(trait_names)) {
    trait_names <- paste0("trait", seq_len(m))
  }
  if (!is.character(trait_names) || length(trait_names) != m ||
      anyNA(trait_names) || anyDuplicated(trait_names)) {
    .stop_invalid_input(
      "trait_names must be %d unique, non-missing names.", m
    )
  }
  if (!is.numeric(inverse_tol) || length(inverse_tol) != 1L ||
      !is.finite(inverse_tol) || inverse_tol <= 0) {
    .stop_invalid_input("inverse_tol must be a single positive finite number.")
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) ||
      ridge < 0) {
    .stop_invalid_input("ridge must be a single non-negative number.")
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    .stop_invalid_input("standardize must be TRUE or FALSE.")
  }

  warnings <- character()
  Sigma_used <- Sigma_P
  if (ridge > 0) {
    Sigma_used <- Sigma_used + ridge * diag(m)
    warnings <- c(warnings, sprintf(
      "ridge = %.6g was added to the diagonal of Sigma_P (user-specified).",
      ridge
    ))
  }
  if (standardize) {
    d <- sqrt(diag(Sigma_used))
    if (any(d <= 0)) {
      .stop_invalid_input(
        "standardize = TRUE requires positive diagonal entries in Sigma_P."
      )
    }
    Sigma_used <- Sigma_used / (d %o% d)
  }
  dimnames(Sigma_used) <- list(trait_names, trait_names)

  inv <- .safe_inverse(Sigma_used, tol = inverse_tol, symmetric = TRUE)
  eig <- inv$eigenvalues
  code <- if (inv$status == "rank_deficient") {
    "rank_deficient"
  } else if (inv$status == "failed") {
    "unstable"
  } else {
    "ok"
  }
  structure(
    list(
      Sigma_P_ref = Sigma_used,
      trait_names = trait_names,
      basis = list(
        inverse_tol = inverse_tol,
        ridge = ridge,
        standardize = standardize
      ),
      status = .new_status(
        ok = code == "ok",
        code = code,
        message = if (code == "ok") "" else sprintf(
          paste0(
            "The reference covariance is %s at inverse_tol = %.3e; ",
            "subset losses will be computed with the pseudo-inverse and ",
            "flagged per locus."
          ),
          code, inverse_tol
        ),
        warnings = warnings
      ),
      diagnostics = list(
        rank = inv$rank,
        min_eigenvalue = min(eig),
        condition_number = inv$condition_number,
        used_pseudoinverse = inv$used_pseudoinverse
      )
    ),
    class = "condped_representation_basis"
  )
}
