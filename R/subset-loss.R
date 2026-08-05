# The unique conditional representation loss path for CondPED v1.0.
#
# .compute_subset_loss() is the ONLY function that computes Gamma, eta,
# Omega, the three Mahalanobis quadratic forms and the representation
# loss rho (frozen interface contract sections 3.4-3.6). Simulation
# truth and estimation must both call this function; the formulas must
# never be re-implemented elsewhere. Stage 1 scope: no subset
# enumeration, no representative-set or irreducible-module extraction,
# no conditional P-values.

#' Conditional representation loss of one representing subset
#'
#' Given a locus effect vector on a candidate trait set A and a fixed
#' reference covariance matrix, computes the conditional representation
#' loss of the representing subset S: the fraction of the full
#' Mahalanobis quadratic form that the traits in R = A \\ S retain
#' after best-linear representation by S.
#'
#' @details
#' With Sigma blocked by (S, R):
#' \deqn{\Gamma_{R|S} = \Sigma_{SS}^{-1}\Sigma_{SR}}
#' \deqn{\eta_{R|S} = \beta_R - \Gamma_{R|S}^{\top}\beta_S}
#' \deqn{\Omega_{R|S} = \Sigma_{RR} - \Sigma_{RS}\Sigma_{SS}^{-1}\Sigma_{SR}}
#' \deqn{\rho(S) = \frac{\eta^{\top}\Omega^{-1}\eta}
#'   {\beta_A^{\top}\Sigma_{AA}^{-1}\beta_A}}
#'
#' Boundary conventions: `rho(empty) = 1` (no Gamma constructed,
#' `subset_qform = 0`, `residual_qform = full_qform`) and
#' `rho(A) = 0` (`residual_qform = 0`, `subset_qform = full_qform`).
#' When `full_qform <= qform_tol` the locus carries no signal and the
#' result is `not_applicable` (a legitimate state: `status$ok` stays
#' `TRUE`). Only floating-point-scale excursions of
#' the loss outside the unit interval are truncated; grossly out-of-range values
#' are kept and marked `unstable`, never silently truncated.
#'
#' Usability rule (frozen): a result may enter formal subset
#' optimisation only when every component inverse is numerically full
#' rank, the quadratic-form decomposition holds within
#' `decomposition_rel_tol` and the loss lies in the unit interval. If
#' any condition fails (`rank_deficient` or `unstable`),
#' `representation_loss` is `NA` and `status$ok` is `FALSE`; the
#' untruncated value remains available as `raw_representation_loss`
#' together with all other diagnostics. Rank-deficient inputs violate
#' the theoretical `Sigma > 0` assumption, so pseudo-inverse losses
#' are never offered as formal results.
#'
#' @param beta Named numeric vector of locus effects on the candidate
#'   traits; `names(beta)` defines the candidate set A.
#' @param Sigma Numeric symmetric `length(beta) x length(beta)`
#'   reference covariance matrix. If dimnamed, it is reordered to match
#'   `names(beta)`; otherwise the order of `beta` is assumed.
#' @param representing_set Subset of the candidate traits (character
#'   names or integer positions); may be empty.
#' @param trait_names Reference trait order; defaults to `names(beta)`.
#' @param inverse_tol Eigenvalue tolerance for all inverses, ranks and
#'   condition numbers (via [`.safe_inverse()`]).
#' @param qform_tol Energy threshold; `full_qform <= qform_tol` yields
#'   `not_applicable`.
#' @param decomposition_rel_tol Maximum tolerated relative error of the
#'   quadratic-form decomposition before the result is `unstable`.
#' @param condition_number_warning Condition numbers above this value
#'   are recorded as warnings.
#'
#' @return A list with components `representing_set`,
#'   `representing_key`, `complement_set`, `complement_key`, `Gamma`,
#'   `eta`, `Omega`, `full_qform`, `subset_qform`, `residual_qform`,
#'   `representation_loss` (the formal loss; `NA` unless the result is
#'   usable, see Details), `raw_representation_loss` (the untruncated
#'   diagnostic value), `decomposition_error`, `rank_Sigma_SS`,
#'   `rank_Omega`, `condition_Sigma_SS`, `condition_Omega`,
#'   `used_pseudoinverse` and `status` (a standard status list with
#'   code `"ok"`, `"rank_deficient"`, `"unstable"` or
#'   `"not_applicable"`; `ok` is `FALSE` for `rank_deficient` and
#'   `unstable`).
#' @keywords internal
.compute_subset_loss <- function(beta,
                                 Sigma,
                                 representing_set,
                                 trait_names = names(beta),
                                 inverse_tol = sqrt(.Machine$double.eps),
                                 qform_tol = 1e-12,
                                 decomposition_rel_tol = 1e-8,
                                 condition_number_warning = 1e10) {
  # ---- input validation -------------------------------------------------
  if (!is.numeric(beta) || length(beta) == 0L || any(!is.finite(beta))) {
    .stop_invalid_input("beta must be a non-empty finite numeric vector.")
  }
  if (is.null(trait_names) && is.null(names(beta))) {
    .stop_invalid_input(
      "beta must be named, or trait_names must be supplied."
    )
  }
  if (is.null(trait_names)) trait_names <- names(beta)
  if (length(trait_names) != length(beta)) {
    .stop_invalid_input(
      "trait_names must have the same length as beta (%d).", length(beta)
    )
  }
  names(beta) <- trait_names
  A <- .normalize_trait_set(trait_names, trait_names, "trait_names")
  if (!is.numeric(Sigma) || !is.matrix(Sigma) ||
      nrow(Sigma) != length(beta) || ncol(Sigma) != length(beta) ||
      any(!is.finite(Sigma))) {
    .stop_invalid_input(
      "Sigma must be a finite numeric %d x %d matrix.",
      length(beta), length(beta)
    )
  }
  if (!is.null(rownames(Sigma))) {
    if (!setequal(rownames(Sigma), A) ||
        !setequal(colnames(Sigma), A)) {
      .stop_invalid_input(
        "dimnames of Sigma must match the candidate traits in beta."
      )
    }
    Sigma <- Sigma[A, A, drop = FALSE]
  }
  scale_max <- max(abs(Sigma))
  if (max(abs(Sigma - t(Sigma))) > 1e-8 * max(1, scale_max)) {
    .stop_invalid_input("Sigma must be symmetric.")
  }
  Sigma <- (Sigma + t(Sigma)) / 2
  for (nm in c("inverse_tol", "qform_tol", "decomposition_rel_tol",
               "condition_number_warning")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0) {
      .stop_invalid_input("%s must be a single positive finite number.", nm)
    }
  }

  S <- .normalize_trait_set(representing_set, A, "representing_set")
  R <- setdiff(A, S)
  idx_S <- match(S, A)
  idx_R <- match(R, A)

  empty_result <- function(status) {
    list(
      representing_set = S, representing_key = .trait_set_key(S),
      complement_set = R, complement_key = .trait_set_key(R),
      Gamma = NULL, eta = NULL, Omega = NULL,
      full_qform = NA_real_, subset_qform = NA_real_,
      residual_qform = NA_real_, representation_loss = NA_real_,
      raw_representation_loss = NA_real_,
      decomposition_error = NA_real_,
      rank_Sigma_SS = NA_integer_, rank_Omega = NA_integer_,
      condition_Sigma_SS = NA_real_, condition_Omega = NA_real_,
      used_pseudoinverse = NA,
      status = status
    )
  }

  # ---- full quadratic form on the candidate set -------------------------
  inv_AA <- .safe_inverse(Sigma, tol = inverse_tol, symmetric = TRUE)
  if (is.null(inv_AA$inverse)) {
    return(empty_result(.new_status(
      ok = FALSE, code = "unstable",
      message = "Failed to invert the candidate-set covariance matrix."
    )))
  }
  full_qform <- sum(beta * (inv_AA$inverse %*% beta))
  if (!is.finite(full_qform) || full_qform <= qform_tol) {
    res <- empty_result(.new_status(
      ok = TRUE, code = "not_applicable",
      message = sprintf(
        "full_qform (%.3e) is not above qform_tol (%.3e); the locus carries no usable signal.",
        full_qform, qform_tol
      )
    ))
    res$full_qform <- full_qform
    res$used_pseudoinverse <- inv_AA$used_pseudoinverse
    return(res)
  }

  warnings <- character()
  note_condition <- function(inv, label) {
    if (is.finite(inv$condition_number) &&
        inv$condition_number > condition_number_warning) {
      warnings <<- c(warnings, sprintf(
        "%s has condition number %.3e (> %.3e).",
        label, inv$condition_number, condition_number_warning
      ))
    }
  }
  note_condition(inv_AA, "Sigma_AA")

  # ---- boundary sets ------------------------------------------------------
  # Usability rule (frozen): a result may enter formal subset optimisation
  # only when every component inverse is numerically full rank, the
  # quadratic-form decomposition holds within tolerance and the loss lies
  # in [0, 1]. Rank-deficient inputs violate the Sigma > 0 assumption:
  # diagnostics are kept, but representation_loss is NA and status$ok is
  # FALSE so downstream code cannot ingest the result by accident.
  if (length(S) == 0L) {
    usable <- inv_AA$status != "rank_deficient"
    return(list(
      representing_set = S, representing_key = .trait_set_key(S),
      complement_set = R, complement_key = .trait_set_key(R),
      Gamma = NULL, eta = NULL, Omega = NULL,
      full_qform = full_qform, subset_qform = 0,
      residual_qform = full_qform,
      representation_loss = if (usable) 1 else NA_real_,
      raw_representation_loss = 1,
      decomposition_error = 0,
      rank_Sigma_SS = 0L, rank_Omega = NA_integer_,
      condition_Sigma_SS = NA_real_, condition_Omega = NA_real_,
      used_pseudoinverse = inv_AA$used_pseudoinverse,
      status = .new_status(
        ok = usable,
        code = if (usable) "ok" else "rank_deficient",
        message = if (usable) "" else
          "Sigma_AA is numerically rank deficient; diagnostics kept, result not usable for optimisation.",
        warnings = warnings
      )
    ))
  }
  if (length(R) == 0L) {
    usable <- inv_AA$status != "rank_deficient"
    return(list(
      representing_set = S, representing_key = .trait_set_key(S),
      complement_set = R, complement_key = .trait_set_key(R),
      Gamma = NULL, eta = NULL, Omega = NULL,
      full_qform = full_qform, subset_qform = full_qform,
      residual_qform = 0,
      representation_loss = if (usable) 0 else NA_real_,
      raw_representation_loss = 0,
      decomposition_error = 0,
      rank_Sigma_SS = inv_AA$rank, rank_Omega = 0L,
      condition_Sigma_SS = inv_AA$condition_number,
      condition_Omega = NA_real_,
      used_pseudoinverse = inv_AA$used_pseudoinverse,
      status = .new_status(
        ok = usable,
        code = if (usable) "ok" else "rank_deficient",
        message = if (usable) "" else
          "Sigma_AA is numerically rank deficient; diagnostics kept, result not usable for optimisation.",
        warnings = warnings
      )
    ))
  }

  # ---- general subset -----------------------------------------------------
  Sigma_SS <- Sigma[idx_S, idx_S, drop = FALSE]
  Sigma_SR <- Sigma[idx_S, idx_R, drop = FALSE]
  Sigma_RR <- Sigma[idx_R, idx_R, drop = FALSE]

  inv_SS <- .safe_inverse(Sigma_SS, tol = inverse_tol, symmetric = TRUE)
  if (is.null(inv_SS$inverse)) {
    return(empty_result(.new_status(
      ok = FALSE, code = "unstable",
      message = "Failed to invert Sigma_SS."
    )))
  }
  note_condition(inv_SS, "Sigma_SS")

  Gamma <- inv_SS$inverse %*% Sigma_SR
  dimnames(Gamma) <- list(S, R)   # |S| x |R|: maps beta_S to predictions of beta_R
  eta <- beta[idx_R] - as.vector(crossprod(Gamma, beta[idx_S]))
  Omega <- Sigma_RR - crossprod(Sigma_SR, inv_SS$inverse %*% Sigma_SR)
  Omega <- (Omega + t(Omega)) / 2
  dimnames(Omega) <- list(R, R)

  inv_Omega <- .safe_inverse(Omega, tol = inverse_tol, symmetric = TRUE)
  if (is.null(inv_Omega$inverse)) {
    return(empty_result(.new_status(
      ok = FALSE, code = "unstable",
      message = "Failed to invert Omega."
    )))
  }
  note_condition(inv_Omega, "Omega")

  subset_qform <- sum(beta[idx_S] * (inv_SS$inverse %*% beta[idx_S]))
  residual_qform <- sum(eta * (inv_Omega$inverse %*% eta))
  decomposition_error <- abs(full_qform - subset_qform - residual_qform) /
    max(full_qform, qform_tol)
  representation_loss <- residual_qform / full_qform

  # Truncate only floating-point-scale excursions; keep and mark anything
  # larger (contract section 5.3).
  fp_tol <- 1e-8
  if (representation_loss < 0 && representation_loss > -fp_tol) {
    representation_loss <- 0
  }
  if (representation_loss > 1 && representation_loss < 1 + fp_tol) {
    representation_loss <- 1
  }

  rank_deficient <- inv_AA$status == "rank_deficient" ||
    inv_SS$status == "rank_deficient" ||
    inv_Omega$status == "rank_deficient"
  used_pinv <- isTRUE(inv_AA$used_pseudoinverse) ||
    isTRUE(inv_SS$used_pseudoinverse) ||
    isTRUE(inv_Omega$used_pseudoinverse)
  out_of_range <- representation_loss < 0 || representation_loss > 1
  decomp_failed <- !is.finite(decomposition_error) ||
    decomposition_error > decomposition_rel_tol

  code <- "ok"
  message <- ""
  if (rank_deficient) {
    code <- "rank_deficient"
    message <- "One or more component matrices are numerically rank deficient; diagnostics kept, result not usable for optimisation."
  } else if (decomp_failed || out_of_range) {
    code <- "unstable"
    message <- sprintf(
      "Decomposition/loss checks failed (decomposition_error = %.3e, loss = %.6g); diagnostics kept, result not usable for optimisation.",
      decomposition_error, representation_loss
    )
  }
  # Usability rule (frozen): full-rank inverses, decomposition within
  # tolerance and loss inside [0, 1]. Otherwise the formal loss is NA and
  # status$ok is FALSE; the raw value stays available for diagnosis.
  usable <- code == "ok"

  list(
    representing_set = S, representing_key = .trait_set_key(S),
    complement_set = R, complement_key = .trait_set_key(R),
    Gamma = Gamma, eta = eta, Omega = Omega,
    full_qform = full_qform,
    subset_qform = subset_qform,
    residual_qform = residual_qform,
    representation_loss = if (usable) representation_loss else NA_real_,
    raw_representation_loss = representation_loss,
    decomposition_error = decomposition_error,
    rank_Sigma_SS = inv_SS$rank,
    rank_Omega = inv_Omega$rank,
    condition_Sigma_SS = inv_SS$condition_number,
    condition_Omega = inv_Omega$condition_number,
    used_pseudoinverse = used_pinv,
    status = .new_status(ok = usable, code = code, message = message,
                         warnings = warnings)
  )
}
