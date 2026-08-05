# Internal utilities for CondPED v0.3.
# These helpers are not exported; see the interface contract
# (docs/design/CondPED_v0.3_Interface_Contract.md), section 1 and 5.

#' Create a standard status object
#'
#' Every public CondPED function returns a status list with fields
#' `ok`, `code`, `message` and `warnings`. This constructor centralises
#' the allowed status codes defined by the interface contract.
#'
#' @param ok Logical scalar; `TRUE` when the call succeeded.
#' @param code Status code, one of `"ok"`, `"invalid_input"`,
#'   `"non_convergence"`, `"rank_deficient"`, `"ill_conditioned"`,
#'   `"bootstrap_failed"`, `"unstable"`, `"empty_selection"`,
#'   `"not_applicable"`, `"too_many_traits"`.
#' @param message Human readable message (character scalar).
#' @param warnings Character vector of collected warnings.
#'
#' @return A list with components `ok`, `code`, `message`, `warnings`.
#' @keywords internal
.new_status <- function(ok = TRUE,
                        code = c("ok", "invalid_input", "non_convergence",
                                 "rank_deficient", "ill_conditioned",
                                 "bootstrap_failed", "unstable",
                                 "empty_selection", "not_applicable",
                                 "too_many_traits"),
                        message = "",
                        warnings = character()) {
  code <- match.arg(code)
  stopifnot(
    "ok must be a single logical value" = is.logical(ok) && length(ok) == 1L && !is.na(ok),
    "message must be a single character string" = is.character(message) && length(message) == 1L,
    "warnings must be a character vector" = is.character(warnings)
  )
  list(ok = ok, code = code, message = message, warnings = warnings)
}

#' Stop with a contracted invalid-input error
#'
#' All user-input validation failures raise an error of class
#' `condped_invalid_input` (contract section 1.3), so callers can
#' distinguish bad input from internal numerical failures.
#'
#' @param ... Passed to [sprintf()] to build the message.
#' @keywords internal
.stop_invalid_input <- function(...) {
  stop(
    structure(
      class = c("condped_invalid_input", "error", "condition"),
      list(message = sprintf(...), call = NULL)
    )
  )
}

#' Validate dimensions of core data matrices
#'
#' Checks that the core data objects are numeric matrices with consistent
#' dimensions: `Y` is n x m, `W` is n x q, `G` is n x p and `K` is n x n.
#' Any subset of the objects may be supplied; `NULL` entries are skipped.
#'
#' @param Y Optional n x m phenotype matrix.
#' @param W Optional n x q fixed-effect design matrix.
#' @param G Optional n x p genotype dosage matrix.
#' @param K Optional n x n genomic relationship matrix.
#'
#' @return Invisibly returns a list with the shared sample size `n` and the
#'   observed dimensions `m`, `q`, `p` (entries `NULL` when the corresponding
#'   matrix was not supplied). Stops with an informative error otherwise.
#' @keywords internal
.validate_dimensions <- function(Y = NULL, W = NULL, G = NULL, K = NULL) {
  mats <- list(Y = Y, W = W, G = G, K = K)
  for (nm in names(mats)) {
    x <- mats[[nm]]
    if (is.null(x)) next
    if (!is.matrix(x) || !is.numeric(x)) {
      .stop_invalid_input("%s must be a numeric matrix.", nm)
    }
    if (anyNA(x)) {
      .stop_invalid_input("%s must not contain missing values.", nm)
    }
    if (any(!is.finite(x))) {
      .stop_invalid_input("%s must not contain non-finite values (NA/Inf/NaN).",
                          nm)
    }
  }
  n_candidates <- c(
    if (!is.null(Y)) nrow(Y),
    if (!is.null(W)) nrow(W),
    if (!is.null(G)) nrow(G),
    if (!is.null(K)) nrow(K)
  )
  if (length(n_candidates) > 0L && length(unique(n_candidates)) > 1L) {
    .stop_invalid_input(
      "Y, W, G and K must share the same number of rows (individuals)."
    )
  }
  if (!is.null(K) && nrow(K) != ncol(K)) {
    .stop_invalid_input("K must be a square n x n matrix.")
  }
  invisible(list(
    n = if (length(n_candidates) > 0L) n_candidates[[1L]] else NULL,
    m = if (!is.null(Y)) ncol(Y) else NULL,
    q = if (!is.null(W)) ncol(W) else NULL,
    p = if (!is.null(G)) ncol(G) else NULL
  ))
}

#' Numerically safe matrix inverse with tolerance control
#'
#' Computes the inverse (or Moore-Penrose generalised inverse) of a square
#' matrix without ever calling an untoleranced [solve()]. The numerical rank
#' is always determined first, from the eigenvalues (symmetric case) or
#' singular values (non-symmetric case) at relative tolerance `tol`; only a
#' numerically positive definite matrix (full rank with every eigenvalue
#' significantly positive) takes the Cholesky fast path. Everything else
#' falls back to a tolerance-controlled eigen/SVD pseudo-inverse. The
#' Moore-Penrose inverse of a zero matrix is the zero matrix itself,
#' reported with rank 0 and status `"rank_deficient"`.
#'
#' @param A Square numeric matrix.
#' @param tol Relative tolerance used to decide whether an eigenvalue is
#'   numerically zero. Defaults to `sqrt(.Machine$double.eps)`.
#' @param symmetric Logical; set `TRUE` (default) when `A` is known to be
#'   symmetric, enabling the Cholesky fast path and a symmetric eigen
#'   decomposition. When `FALSE`, an SVD-based pseudo-inverse is used.
#' @param allow_pseudoinverse Logical; when `FALSE`, a numerically
#'   rank-deficient matrix returns `inverse = NULL` with status
#'   `"rank_deficient"` instead of the Moore-Penrose generalised inverse.
#'
#' @return A list with components:
#'   \describe{
#'     \item{inverse}{The inverse or pseudo-inverse matrix; `NULL` on failure
#'       or when `allow_pseudoinverse = FALSE` meets a rank-deficient matrix.}
#'     \item{rank}{Integer; numerical rank at tolerance `tol`.}
#'     \item{eigenvalues}{Numeric vector of eigenvalues (symmetric case) or
#'       singular values (non-symmetric case).}
#'     \item{condition_number}{Ratio of largest to smallest nonzero singular
#'       value/eigenvalue magnitude; `Inf` for rank-deficient input.}
#'     \item{method}{`"chol"` or `"eigen_pinv"`.}
#'     \item{used_pseudoinverse}{Logical; `TRUE` when the returned inverse
#'       was built by the eigen/SVD pseudo-inverse path rather than
#'       Cholesky; `NA` on failure.}
#'     \item{status}{`"ok"`, `"rank_deficient"` or `"failed"`.}
#'   }
#' @keywords internal
.safe_inverse <- function(A,
                          tol = sqrt(.Machine$double.eps),
                          symmetric = TRUE,
                          allow_pseudoinverse = TRUE) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) {
    stop("A must be a square matrix.", call. = FALSE)
  }
  if (!is.numeric(A) || anyNA(A)) {
    stop("A must be numeric and free of missing values.", call. = FALSE)
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a single positive finite number.", call. = FALSE)
  }
  n <- nrow(A)

  failed <- list(
    inverse = NULL,
    rank = NA_integer_,
    eigenvalues = rep(NA_real_, n),
    condition_number = NA_real_,
    method = "eigen_pinv",
    used_pseudoinverse = NA,
    status = "failed"
  )
  # Rank-deficient input with pseudo-inverses disabled: report the rank
  # honestly but return no inverse.
  not_allowed <- function(values) {
    list(
      inverse = NULL,
      rank = sum(abs(values) > tol * max(abs(values))),
      eigenvalues = values,
      condition_number = Inf,
      method = "eigen_pinv",
      used_pseudoinverse = FALSE,
      status = "rank_deficient"
    )
  }

  if (symmetric) {
    # Always determine the numerical rank from the eigenvalues first:
    # chol() alone can succeed on mathematically positive definite but
    # numerically near-singular matrices (e.g. diag(c(1, 1e-20))).
    decomp <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
    if (is.null(decomp)) return(failed)
    values <- decomp$values
    vectors <- decomp$vectors
    max_abs <- max(abs(values))
    if (max_abs == 0) {
      # Zero spectrum: the Moore-Penrose inverse of a zero matrix is the
      # zero matrix of the same dimension.
      if (!allow_pseudoinverse) {
        return(list(
          inverse = NULL,
          rank = 0L,
          eigenvalues = values,
          condition_number = Inf,
          method = "eigen_pinv",
          used_pseudoinverse = FALSE,
          status = "rank_deficient"
        ))
      }
      return(list(
        inverse = matrix(0, n, n),
        rank = 0L,
        eigenvalues = values,
        condition_number = Inf,
        method = "eigen_pinv",
        used_pseudoinverse = TRUE,
        status = "rank_deficient"
      ))
    }
    cutoff <- tol * max_abs
    keep <- abs(values) > cutoff
    rank <- sum(keep)

    if (rank < n && !allow_pseudoinverse) return(not_allowed(values))

    if (rank == n && min(values) > cutoff) {
      # Numerically positive definite: Cholesky fast path.
      ch <- tryCatch(chol(A), error = function(e) NULL)
      if (!is.null(ch)) {
        return(list(
          inverse = chol2inv(ch),
          rank = n,
          eigenvalues = values,
          condition_number = max(values) / min(values),
          method = "chol",
          used_pseudoinverse = FALSE,
          status = "ok"
        ))
      }
      # Cholesky failed despite eigenvalue clearance: fall through to the
      # eigen pseudo-inverse (still full rank, so status stays "ok").
    }

    d_inv <- numeric(n)
    d_inv[keep] <- 1 / values[keep]
    inverse <- vectors %*% (d_inv * t(vectors))
    return(list(
      inverse = inverse,
      rank = rank,
      eigenvalues = values,
      condition_number = if (rank < n) Inf else max_abs / min(abs(values)),
      method = "eigen_pinv",
      used_pseudoinverse = TRUE,
      status = if (rank < n) "rank_deficient" else "ok"
    ))
  }

  # Non-symmetric path: SVD-based Moore-Penrose pseudo-inverse.
  decomp <- tryCatch(svd(A), error = function(e) NULL)
  if (is.null(decomp)) return(failed)
  d <- decomp$d
  if (max(d) == 0) {
    # Zero matrix: Moore-Penrose inverse is the zero matrix.
    if (!allow_pseudoinverse) {
      return(list(
        inverse = NULL,
        rank = 0L,
        eigenvalues = d,
        condition_number = Inf,
        method = "eigen_pinv",
        used_pseudoinverse = FALSE,
        status = "rank_deficient"
      ))
    }
    return(list(
      inverse = matrix(0, n, n),
      rank = 0L,
      eigenvalues = d,
      condition_number = Inf,
      method = "eigen_pinv",
      used_pseudoinverse = TRUE,
      status = "rank_deficient"
    ))
  }
  cutoff <- tol * max(d)
  keep <- d > cutoff
  rank <- sum(keep)
  if (rank < n && !allow_pseudoinverse) return(not_allowed(d))
  d_inv <- numeric(n)
  d_inv[keep] <- 1 / d[keep]
  inverse <- decomp$v %*% (d_inv * t(decomp$u))
  list(
    inverse = inverse,
    rank = rank,
    eigenvalues = d,
    condition_number = if (rank < n) Inf else max(d) / min(d),
    method = "eigen_pinv",
    used_pseudoinverse = TRUE,
    status = if (rank < n) "rank_deficient" else "ok"
  )
}

#' Genomic relationship matrix from standardised genotypes
#'
#' Computes `K = Z Z' / p` from a standardised genotype matrix `Z`
#' (columns centred and scaled) and rescales it to average diagonal 1.
#'
#' @param Z `n x p` standardised genotype matrix.
#' @return Symmetric `n x n` relationship matrix with `mean(diag(K)) = 1`.
#' @keywords internal
.make_grm <- function(Z) {
  Z <- as.matrix(Z)
  K <- tcrossprod(Z) / ncol(Z)
  K <- (K + t(K)) / 2
  K / mean(diag(K))
}

#' Signal that a public function is not implemented yet
#'
#' Used by the S0 API skeleton: every public stub throws an error of class
#' `condped_not_implemented` so that unfinished code paths can never be
#' mistaken for real results.
#'
#' @param name Character scalar; name of the public function.
#'
#' @return Never returns; always signals an error.
#' @keywords internal
.not_implemented <- function(name) {
  cond <- structure(
    class = c("condped_not_implemented", "error", "condition"),
    list(
      message = sprintf(
        "%s() belongs to the CondPED v0.3 API skeleton (S0) and is not implemented yet.",
        name
      ),
      call = NULL
    )
  )
  stop(cond)
}
