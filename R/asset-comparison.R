# Stage 6C: ASSET comparator adapter.
#
# `.run_asset_comparison()` is the single boundary between CondPED
# per-signal trait-wise estimates and the external ASSET package (an
# optional, suggested-only dependency). The adapter owns ALL
# normalisation: downstream evaluators only ever see the frozen output
# schema (`asset_p`, `asset_best_subset`, `asset_positive_subset`,
# `asset_negative_subset`, `status`) and never touch ASSET internals.
#
# Status codes (character scalar in `$status`):
#   "ok"                  backend ran and its result was normalised
#   "asset_not_available" ASSET is not installed (backend = NULL); never
#                         an error, so comparison pipelines stay runnable
#   "asset_failed"        backend errored or returned an unrecognisable
#                         structure; p is NA, subsets empty
#   "invalid_input"       malformed adapter inputs; thrown via
#                         `.stop_invalid_input()` (contract section 1.3)
#
# The adapter never fabricates rho/Rep/Irr from ASSET output: it reports
# exactly the best-subset p-value and the trait subset, plus a split of
# that subset by the sign of the CondPED effect estimates.

#' Frozen empty result for the ASSET adapter
#'
#' @param status Status code to stamp on the result.
#' @return List with the frozen ASSET comparison schema.
#' @keywords internal
.asset_empty_result <- function(status) {
  list(
    asset_p = NA_real_,
    asset_best_subset = character(),
    asset_positive_subset = character(),
    asset_negative_subset = character(),
    status = status
  )
}

#' Validate the per-signal inputs of the ASSET adapter
#'
#' All failures raise `condped_invalid_input`. Reorders `beta` (and
#' `se`/`z` when supplied) to `trait_names` and derives `z = beta / se`
#' when only `se` is given.
#'
#' @param beta,se,z,Sigma_Z,trait_names,two_sided See
#'   [`.run_asset_comparison()`].
#' @return List with validated `beta`, `se`, `z`, `Sigma_Z`,
#'   `trait_names`, `two_sided` (all aligned to `trait_names`).
#' @keywords internal
.validate_asset_inputs <- function(beta, se, z, Sigma_Z, trait_names,
                                   two_sided) {
  if (!is.character(trait_names) || length(trait_names) < 1L ||
      anyNA(trait_names) || anyDuplicated(trait_names)) {
    .stop_invalid_input("trait_names must be unique, non-missing names.")
  }
  m <- length(trait_names)

  if (!is.numeric(beta) || length(beta) != m || any(!is.finite(beta)) ||
      is.null(names(beta)) || anyNA(names(beta))) {
    .stop_invalid_input(
      "beta must be a named finite numeric vector with one entry per trait."
    )
  }
  if (!setequal(names(beta), trait_names)) {
    .stop_invalid_input("names(beta) must match trait_names.")
  }
  beta <- beta[trait_names]

  if (!is.null(se)) {
    if (!is.numeric(se) || length(se) != m || any(!is.finite(se)) ||
        any(se <= 0) || is.null(names(se)) ||
        !setequal(names(se), trait_names)) {
      .stop_invalid_input(
        "se must be a named, positive, finite numeric vector matching trait_names."
      )
    }
    se <- se[trait_names]
  }
  if (!is.null(z)) {
    if (!is.numeric(z) || length(z) != m || any(!is.finite(z)) ||
        is.null(names(z)) || !setequal(names(z), trait_names)) {
      .stop_invalid_input(
        "z must be a named finite numeric vector matching trait_names."
      )
    }
    z <- z[trait_names]
  } else {
    if (is.null(se)) {
      .stop_invalid_input(
        "Either z or se must be supplied (z is derived as beta / se)."
      )
    }
    z <- beta / se
  }

  if (!is.numeric(Sigma_Z) || !is.matrix(Sigma_Z) ||
      nrow(Sigma_Z) != m || ncol(Sigma_Z) != m ||
      any(!is.finite(Sigma_Z))) {
    .stop_invalid_input(
      "Sigma_Z must be a finite numeric %d x %d matrix.", m, m
    )
  }
  if (is.null(dimnames(Sigma_Z))) {
    dimnames(Sigma_Z) <- list(trait_names, trait_names)
  }
  if (!setequal(rownames(Sigma_Z), trait_names) ||
      !setequal(colnames(Sigma_Z), trait_names)) {
    .stop_invalid_input("dimnames of Sigma_Z must match trait_names.")
  }
  Sigma_Z <- Sigma_Z[trait_names, trait_names, drop = FALSE]
  scale_max <- max(1, max(abs(Sigma_Z)))
  if (max(abs(Sigma_Z - t(Sigma_Z))) > 1e-8 * scale_max) {
    .stop_invalid_input("Sigma_Z must be symmetric.")
  }
  Sigma_Z <- (Sigma_Z + t(Sigma_Z)) / 2
  if (max(abs(diag(Sigma_Z) - 1)) > 1e-6) {
    .stop_invalid_input("Sigma_Z must be a correlation matrix (diag == 1).")
  }
  eig <- tryCatch(
    eigen(Sigma_Z, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(eig) || min(eig) < -1e-8 * scale_max) {
    .stop_invalid_input("Sigma_Z must be (numerically) positive semidefinite.")
  }

  if (!is.logical(two_sided) || length(two_sided) != 1L ||
      is.na(two_sided)) {
    .stop_invalid_input("two_sided must be TRUE or FALSE.")
  }

  list(beta = beta, se = se, z = z, Sigma_Z = Sigma_Z,
       trait_names = trait_names, two_sided = two_sided)
}

#' Call the real ASSET package defensively
#'
#' Only reached when `backend = NULL` and ASSET is installed. Uses the
#' official `ASSET::fast_asset()` with the frozen mapping
#' trait_names -> traits.lab, beta -> beta.hat, se -> sigma.hat,
#' sample_size -> Neff, Sigma_Z -> cor, and blocks built by
#' `ASSET::create_blocks(Sigma_Z)`. `scr_pthr` is fixed at the official
#' default 0.05 (recorded in the output settings).
#'
#' @param beta,se,z,Sigma_Z,trait_names,two_sided Validated inputs.
#' @param sample_size Validated per-trait sample sizes.
#' @return Raw `fast_asset()` result (opaque to the caller).
#' @keywords internal
.call_asset_package <- function(beta, se, z, Sigma_Z, trait_names,
                                two_sided, sample_size) {
  m <- length(trait_names)
  neff <- if (length(sample_size) == 1L) {
    rep(sample_size, m)
  } else {
    unname(sample_size[trait_names])
  }
  block <- ASSET::create_blocks(Sigma_Z)
  ASSET::fast_asset(
    snp = "signal",
    traits.lab = trait_names,
    beta.hat = matrix(beta, nrow = 1L,
                      dimnames = list("signal", trait_names)),
    sigma.hat = matrix(se, nrow = 1L,
                       dimnames = list("signal", trait_names)),
    Neff = neff,
    cor = Sigma_Z,
    block = block,
    scr_pthr = 0.05
  )
}

#' Normalise a raw ASSET-like result into the frozen schema
#'
#' Understands two raw shapes:
#' \enumerate{
#'   \item a list carrying a p-value under one of `p`, `p_value`,
#'     `p.value`, `pval`, `best_p` and the best trait subset under one of
#'     `subset`, `best_subset`, `traits`, `trait_subset`,
#'     `selected_traits`;
#'   \item a data.frame with one row per searched subset, a p-value
#'     column (`p.value`/`p_value`/`pval`/`p`) and one column per trait
#'     marking membership (non-`NA`, non-`0`); the row with the smallest
#'     p-value wins.
#' }
#' The positive/negative split is derived from the sign of the CondPED
#' effect estimates (`beta`), never from ASSET internals. Anything else
#' raises an error (mapped to status `"asset_failed"` by the caller).
#'
#' @param raw Raw backend result.
#' @param beta Named numeric effect vector (aligned to `trait_names`).
#' @param trait_names Trait names.
#' @return List with the frozen schema and status `"ok"`.
#' @keywords internal
.normalize_asset_result <- function(raw, beta, trait_names,
                                    two_sided = TRUE) {
  if (is.null(raw)) stop("backend returned NULL.")

  # Official ASSET::fast_asset() shape: take ASSET's own sign-split
  # subsets; asset_best_subset = union(positive, negative) — never
  # interpreted as a CondPED representative set.
  if (is.list(raw) && !is.data.frame(raw) &&
      any(c("Subset.2sided", "Subset.1sided") %in% names(raw))) {
    sub <- if (isTRUE(two_sided) && !is.null(raw$Subset.2sided)) {
      raw$Subset.2sided
    } else {
      raw$Subset.1sided
    }
    if (is.null(sub) || is.null(sub$pval) ||
        !is.finite(as.numeric(sub$pval)[1L])) {
      stop("fast_asset result carries no finite p-value.")
    }
    take <- function(ph) {
      if (is.null(ph)) return(character())
      colnames(ph)[as.logical(ph[1L, ])]
    }
    pos <- take(sub$pheno.1)
    neg <- take(sub$pheno.2)
    best <- trait_names[trait_names %in% union(pos, neg)]
    return(list(
      asset_p = as.numeric(sub$pval)[1L],
      asset_best_subset = best,
      asset_positive_subset = pos,
      asset_negative_subset = neg,
      status = "ok"
    ))
  }

  pick_field <- function(x, candidates) {
    for (nm in candidates) {
      if (!is.null(x[[nm]])) return(x[[nm]])
    }
    NULL
  }

  p_raw <- NULL
  subset_raw <- NULL
  if (is.data.frame(raw)) {
    pcol <- intersect(
      c("p.value", "p_value", "pval", "p", "P.value", "P"),
      names(raw)
    )
    if (length(pcol) == 0L) {
      stop("ASSET data.frame result has no recognisable p-value column.")
    }
    pcol <- pcol[1L]
    pvec <- as.numeric(raw[[pcol]])
    if (!any(is.finite(pvec))) stop("no finite p-value in ASSET result.")
    best_row <- which.min(ifelse(is.finite(pvec), pvec, Inf))
    p_raw <- pvec[best_row]
    trait_cols <- intersect(trait_names, names(raw))
    in_set <- vapply(trait_cols, function(tc) {
      v <- raw[[tc]][best_row]
      if (is.logical(v)) return(isTRUE(v))
      if (is.numeric(v)) return(!is.na(v) && v != 0)
      !is.na(v)
    }, logical(1))
    subset_raw <- trait_cols[in_set]
  } else if (is.list(raw)) {
    p_raw <- pick_field(raw, c("p", "p_value", "p.value", "pval",
                               "best_p"))
    subset_raw <- pick_field(raw, c("subset", "best_subset", "traits",
                                    "trait_subset", "selected_traits"))
  } else {
    stop("unrecognised ASSET result of class: ",
         paste(class(raw), collapse = ", "))
  }

  if (is.null(p_raw) || length(p_raw) != 1L || !is.finite(as.numeric(p_raw))) {
    stop("backend result carries no single finite p-value.")
  }
  asset_p <- as.numeric(p_raw)
  if (asset_p < 0 || asset_p > 1) {
    stop("backend p-value outside [0, 1].")
  }

  if (is.null(subset_raw)) {
    best <- character()
  } else {
    best <- unique(as.character(subset_raw))
    unknown <- setdiff(best, trait_names)
    if (length(unknown) > 0L) {
      stop("backend subset contains unknown trait(s): ", unknown[1L])
    }
    best <- trait_names[trait_names %in% best]  # canonical trait order
  }

  sgn <- sign(beta[best])
  list(
    asset_p = asset_p,
    asset_best_subset = best,
    asset_positive_subset = best[sgn > 0],
    asset_negative_subset = best[sgn < 0],
    status = "ok"
  )
}

#' Run one signal's trait-wise effects through ASSET
#'
#' Adapter between CondPED per-signal trait-wise estimates and the
#' external ASSET package. ASSET is an optional dependency: when
#' `backend = NULL` and the package is not installed, the adapter returns
#' the frozen schema with status `"asset_not_available"` instead of
#' erroring, so comparison pipelines degrade gracefully. The trait
#' z-score correlation matrix `Sigma_Z` is REQUIRED and is always passed
#' to the backend unchanged — correlated traits are fully supported and
#' the correlation is never diagonalised away.
#'
#' @param beta Named numeric vector of trait effects for ONE signal
#'   (names must match `trait_names`).
#' @param se Optional named numeric vector of standard errors
#'   (positive, finite).
#' @param z Optional named numeric vector of z-scores; derived as
#'   `beta / se` when absent (then `se` is required).
#' @param Sigma_Z `m x m` correlation matrix of the trait z-scores
#'   (symmetric, unit diagonal, numerically PSD; dimnames must match
#'   `trait_names` or be `NULL`).
#' @param trait_names Trait names (character, unique).
#' @param sample_size Total analysed sample size per trait (positive
#'   scalar, or one value per trait); mapped to `Neff` in the official
#'   `ASSET::fast_asset()` call. Required when `backend = NULL` (the
#'   real ASSET path).
#' @param two_sided Logical; forwarded to the backend (two-sided subset
#'   search when `TRUE`).
#' @param backend Test hook: a
#'   `function(beta, se, z, Sigma_Z, trait_names, two_sided)` returning a
#'   raw ASSET-like result (see [`.normalize_asset_result()`] for the
#'   accepted shapes). When `NULL`, the adapter looks for the installed
#'   ASSET package and calls its Z-score subset-search entry point
#'   defensively; any failure yields status `"asset_failed"`.
#'
#' @return A list with the frozen comparison schema:
#' \describe{
#'   \item{asset_p}{Best-subset p-value (`NA_real_` unless status is
#'     `"ok"`).}
#'   \item{asset_best_subset}{Union of ASSET's positive and negative
#'     subsets (canonical `trait_names` order); never interpreted as a
#'     CondPED representative set.}
#'   \item{asset_positive_subset}{Best-subset traits with positive
#'     effects (ASSET's own sign split on the fast_asset path; the
#'     CondPED beta signs for injected backends).}
#'   \item{asset_negative_subset}{Best-subset traits with negative
#'     effects.}
#'   \item{status}{`"ok"`, `"asset_not_available"` or `"asset_failed"`.
#'     Malformed inputs raise `condped_invalid_input` instead.}
#'   \item{settings}{List recording `scr_pthr` (fixed at the official
#'     default 0.05), `two_sided` and `Neff`.}
#' }
#' @keywords internal
.run_asset_comparison <- function(beta, se = NULL, z = NULL, Sigma_Z,
                                  trait_names, sample_size = NULL,
                                  two_sided = TRUE,
                                  backend = NULL) {
  if (!is.null(backend) && !is.function(backend)) {
    .stop_invalid_input("backend must be NULL or a function.")
  }
  v <- .validate_asset_inputs(beta, se, z, Sigma_Z, trait_names,
                              two_sided)

  if (is.null(backend)) {
    if (!requireNamespace("ASSET", quietly = TRUE)) {
      return(c(.asset_empty_result("asset_not_available"),
               list(settings = list(scr_pthr = 0.05,
                                    two_sided = isTRUE(two_sided),
                                    Neff = sample_size))))
    }
    # the official fast_asset() needs sigma.hat and Neff
    if (is.null(v$se)) {
      .stop_invalid_input(
        "se is required for the real ASSET::fast_asset() path (sigma.hat)."
      )
    }
    if (is.null(sample_size) || !is.numeric(sample_size) ||
        any(!is.finite(sample_size)) || any(sample_size <= 0) ||
        !(length(sample_size) %in% c(1L, length(v$trait_names)))) {
      .stop_invalid_input(
        paste0("sample_size (total analysed sample size per trait) is ",
               "required for the real ASSET path: a positive scalar or ",
               "one value per trait.")
      )
    }
    backend <- function(beta, se, z, Sigma_Z, trait_names, two_sided) {
      .call_asset_package(beta, se, z, Sigma_Z, trait_names, two_sided,
                          sample_size)
    }
  }

  raw <- tryCatch(
    backend(beta = v$beta, se = v$se, z = v$z, Sigma_Z = v$Sigma_Z,
            trait_names = v$trait_names, two_sided = v$two_sided),
    error = function(e) e
  )
  if (inherits(raw, "error")) {
    return(c(.asset_empty_result("asset_failed"),
             list(settings = list(scr_pthr = 0.05,
                                  two_sided = isTRUE(two_sided),
                                  Neff = sample_size))))
  }

  norm <- tryCatch(
    .normalize_asset_result(raw, v$beta, v$trait_names, v$two_sided),
    error = function(e) e
  )
  if (inherits(norm, "error")) {
    return(c(.asset_empty_result("asset_failed"),
             list(settings = list(scr_pthr = 0.05,
                                  two_sided = isTRUE(two_sided),
                                  Neff = sample_size))))
  }
  c(norm, list(settings = list(scr_pthr = 0.05,
                               two_sided = isTRUE(v$two_sided),
                               Neff = sample_size)))
}
