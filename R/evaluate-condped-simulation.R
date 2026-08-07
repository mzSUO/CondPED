#' Evaluate CondPED simulation replicates (v1.0)
#'
#' Computes per-replicate recovery metrics and grouped summaries with
#' Monte Carlo standard errors. All set comparisons use normalised
#' trait keys; families are compared as whole sets of keys (never just
#' the first solution). The empty-set rules of the frozen contract are
#' applied verbatim. No legacy D/profile/persistent metrics exist.
#'
#' @details
#' ## Replicate file format
#'
#' Each entry of `files` is an RDS containing one replicate result:
#' \describe{
#'   \item{truth}{The `truth` component of a
#'     [simulate_condped_data()] call.}
#'   \item{estimates}{A list with `causal_marker`, `omnibus_summary`
#'     (data.frame with `marker_id` and `p_value`), `candidate_sets`
#'     (named list marker -> traits), `minimum_representative_sets`
#'     and `irreducible_modules` (data.frames as returned by
#'     [decompose_conditional_effects()]).}
#'   \item{settings}{List with the scenario descriptors named in
#'     `group_by` plus `tolerance` and `alpha_omnibus`.}
#'   \item{status}{Standard status list.}
#'   \item{runtime}{Numeric seconds.}
#' }
#'
#' ## Metric definitions (per replicate, primary tolerance)
#'
#' \describe{
#'   \item{type1_omnibus}{Fraction of non-causal markers with
#'     `p_value <= alpha_omnibus`.}
#'   \item{power_omnibus}{Indicator that the causal marker passes;
#'     `NA` for the null architecture.}
#'   \item{candidate_tpr, candidate_fdp, candidate_exact_recovery,
#'     candidate_jaccard}{TPR / FDP / exact recovery / Jaccard of the
#'     estimated candidate set against `truth$candidate_traits`.}
#'   \item{representative_exact_recovery, representative_jaccard,
#'     representative_size_bias, tie_coverage}{Family exact recovery,
#'     family Jaccard, set-size bias (estimated minus truth minimum
#'     size) and `tie_coverage`: the fraction of truth tied minimum
#'     sets whose key appears in the estimated family.}
#'   \item{irreducible_tpr, irreducible_fdp,
#'     irreducible_exact_recovery, irreducible_jaccard}{TPR / FDP /
#'     family exact recovery / family Jaccard over module keys.}
#'   \item{joint_module_miss_rate}{Fraction of truth modules of size
#'     greater than 1 whose key is not exactly recovered (`NA` when
#'     the truth has no joint module).}
#'   \item{over_fragmentation_rate}{Fraction of truth joint modules
#'     for which the estimated family contains two or more strict
#'     sub-modules (`NA` when the truth has no joint module).}
#'   \item{overall_recovery}{Candidate exact AND representative family
#'     exact AND module family exact.}
#'   \item{conditional_recovery}{Representative family exact AND
#'     module family exact, evaluated only when the candidate set was
#'     exactly recovered (`NA` otherwise).}
#' }
#'
#' Empty-set rules (frozen): truth empty + estimate empty -> exact = 1,
#' Jaccard = 1, TPR = NA, FDP = 0; truth non-empty + estimate empty ->
#' TPR = 0, FDP = 0, Jaccard = 0, exact = 0; truth empty + estimate
#' non-empty -> TPR = NA, FDP = 1, Jaccard = 0, exact = 0.
#'
#' @param files Character vector of RDS file paths, one per replicate.
#' @param metrics Metrics to report; defaults to the full frozen list.
#' @param include_unstable Logical; keep replicates whose status code
#'   is `unstable` in the metric summaries (they are always listed in
#'   `failures`).
#' @param group_by Scenario columns to group summaries by; columns
#'   absent from the replicate settings are dropped with a warning.
#' @param conf_level Confidence level for the summary intervals.
#'
#' @return A list with components `replicate_metrics` (one row per
#'   evaluated replicate), `summary` (grouped means), `mcse` (grouped
#'   Monte Carlo standard errors), `failures` (read failures and
#'   non-ok replicates), `settings` and `status`.
#' @export
evaluate_condped_simulation <- function(
  files,
  metrics = c(
    "type1_omnibus",
    "power_omnibus",
    "candidate_tpr",
    "candidate_fdp",
    "candidate_exact_recovery",
    "candidate_jaccard",
    "representative_exact_recovery",
    "representative_jaccard",
    "representative_size_bias",
    "tie_coverage",
    "irreducible_tpr",
    "irreducible_fdp",
    "irreducible_exact_recovery",
    "irreducible_jaccard",
    "joint_module_miss_rate",
    "over_fragmentation_rate",
    "overall_recovery",
    "conditional_recovery",
    "runtime"
  ),
  include_unstable = TRUE,
  group_by = c(
    "experiment",
    "n",
    "architecture",
    "locus_pve",
    "correlation",
    "target_loss",
    "tolerance"
  ),
  conf_level = 0.95
) {
  metric_all <- c(
    "type1_omnibus", "power_omnibus",
    "candidate_tpr", "candidate_fdp", "candidate_exact_recovery",
    "candidate_jaccard", "representative_exact_recovery",
    "representative_jaccard", "representative_size_bias",
    "tie_coverage", "irreducible_tpr", "irreducible_fdp",
    "irreducible_exact_recovery", "irreducible_jaccard",
    "joint_module_miss_rate", "over_fragmentation_rate",
    "overall_recovery", "conditional_recovery", "runtime"
  )
  metrics <- match.arg(metrics, several.ok = TRUE)
  if (!is.logical(include_unstable) || length(include_unstable) != 1L ||
      is.na(include_unstable)) {
    .stop_invalid_input("include_unstable must be TRUE or FALSE.")
  }
  .check_prob(conf_level, "conf_level")
  if (!is.character(files) || length(files) == 0L) {
    .stop_invalid_input("files must be a non-empty character vector of paths.")
  }

  warnings <- character()
  rep_rows <- list()
  fail_rows <- list()

  for (i in seq_along(files)) {
    f <- files[[i]]
    obj <- tryCatch(readRDS(f), error = function(e) e)
    if (inherits(obj, "error")) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        file = f, reason = paste("read error:", conditionMessage(obj)),
        included = FALSE, stringsAsFactors = FALSE
      )
      next
    }
    if (is.null(obj$truth) || is.null(obj$estimates)) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        file = f, reason = "missing truth or estimates component",
        included = FALSE, stringsAsFactors = FALSE
      )
      next
    }
    ok <- isTRUE(obj$status$ok)
    code <- if (!is.null(obj$status$code)) obj$status$code else "ok"
    unstable <- !ok || code != "ok"
    included <- !unstable || isTRUE(include_unstable)
    if (unstable) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        file = f, reason = sprintf("replicate status: %s", code),
        included = included, stringsAsFactors = FALSE
      )
    }
    if (!included) next
    row <- tryCatch(
      .evaluate_one_replicate(obj$truth, obj$estimates,
                              settings = obj$settings,
                              runtime = obj$runtime),
      error = function(e) e
    )
    if (inherits(row, "error")) {
      fail_rows[[length(fail_rows) + 1L]] <- data.frame(
        file = f,
        reason = paste("evaluation error:", conditionMessage(row)),
        included = FALSE, stringsAsFactors = FALSE
      )
      next
    }
    row$file <- f
    rep_rows[[length(rep_rows) + 1L]] <- row
  }

  replicate_metrics <- if (length(rep_rows) > 0L) {
    out <- do.call(rbind, rep_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame()
  }
  failures <- if (length(fail_rows) > 0L) {
    out <- do.call(rbind, fail_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame(file = character(), reason = character(),
               included = logical(), stringsAsFactors = FALSE)
  }

  # ---- grouped summary -------------------------------------------------------
  # A group_by column is honoured only when at least one replicate
  # actually carries it; all-missing columns are dropped with a warning.
  present <- intersect(group_by, names(replicate_metrics))
  informative <- present[vapply(present, function(g) {
    nrow(replicate_metrics) > 0L && !all(is.na(replicate_metrics[[g]]))
  }, logical(1))]
  dropped <- setdiff(group_by, informative)
  if (length(dropped) > 0L && nrow(replicate_metrics) > 0L) {
    msg <- sprintf(
      "group_by column(s) not present in replicate settings and dropped: %s.",
      paste(dropped, collapse = ", ")
    )
    warnings <- c(warnings, msg)
    warning(msg, call. = FALSE)
  }
  group_cols <- informative
  summary_df <- data.frame()
  mcse_df <- data.frame()
  if (nrow(replicate_metrics) > 0L) {
    summ <- .summarise_metrics(replicate_metrics, metrics, group_cols,
                               conf_level)
    summary_df <- summ$summary
    mcse_df <- summ$mcse
  }

  list(
    replicate_metrics = replicate_metrics,
    summary = summary_df,
    mcse = mcse_df,
    failures = failures,
    settings = list(
      metrics = metrics,
      include_unstable = include_unstable,
      group_by = group_by,
      conf_level = conf_level,
      n_files = length(files)
    ),
    status = .new_status(
      ok = TRUE,
      code = if (length(rep_rows) == 0L) "empty_selection" else "ok",
      message = if (length(rep_rows) == 0L) {
        "No evaluable replicates."
      } else {
        ""
      },
      warnings = warnings
    )
  )
}

#' Per-replicate metrics against the population truth
#'
#' @param truth `simulate_condped_data()` truth component.
#' @param estimates Estimates list (see
#'   [evaluate_condped_simulation()]).
#' @param settings Replicate settings (scenario descriptors).
#' @param runtime Numeric seconds.
#' @return A one-row data.frame of metrics plus scenario columns.
#' @keywords internal
.evaluate_one_replicate <- function(truth, estimates,
                                    settings = list(),
                                    runtime = NA_real_) {
  tol <- if (!is.null(settings$tolerance)) settings$tolerance else 0.10
  alpha <- if (!is.null(settings$alpha_omnibus)) {
    settings$alpha_omnibus
  } else {
    0.05
  }
  causal <- estimates$causal_marker
  if (is.null(causal)) {
    causal <- rownames(truth$beta)[1L]
  }

  # ---- omnibus ---------------------------------------------------------------
  om <- estimates$omnibus_summary
  is_null_scene <- length(truth$candidate_traits) == 0L
  type1 <- NA_real_
  power <- NA_real_
  if (!is.null(om) && nrow(om) > 0L) {
    tested <- !is.na(om$p_value)
    null_rows <- tested & om$marker_id != causal
    if (any(null_rows)) {
      type1 <- mean(om$p_value[null_rows] <= alpha)
    }
    if (!is_null_scene) {
      p_causal <- om$p_value[om$marker_id == causal & tested]
      if (length(p_causal) > 0L) power <- as.numeric(p_causal[1L] <= alpha)
    }
  }

  # ---- candidate set -----------------------------------------------------------
  est_cand <- estimates$candidate_sets[[causal]]
  if (is.null(est_cand)) est_cand <- character()
  cand <- .set_family_metrics(truth$candidate_traits, est_cand)

  # ---- representative set family (primary tolerance) ----------------------------
  truth_reps <- .family_keys_at(truth$minimum_representative_sets, tol)
  est_reps <- .family_keys_at(estimates$minimum_representative_sets, tol)
  rep_fam <- .set_family_metrics(truth_reps, est_reps)
  truth_min_size <- .min_size_at(truth$minimum_representative_sets, tol)
  est_min_size <- .min_size_at(estimates$minimum_representative_sets, tol)
  size_bias <- if (is.finite(truth_min_size) && is.finite(est_min_size)) {
    est_min_size - truth_min_size
  } else {
    NA_real_
  }
  tie_coverage <- if (length(truth_reps) == 0L) {
    if (length(est_reps) == 0L) 1 else 0
  } else {
    length(intersect(truth_reps, est_reps)) / length(truth_reps)
  }

  # ---- irreducible module family -------------------------------------------------
  truth_mods <- .family_keys_at(truth$irreducible_modules, tol)
  est_mods <- .family_keys_at(estimates$irreducible_modules, tol)
  mod_fam <- .set_family_metrics(truth_mods, est_mods)

  truth_mod_sets <- .family_sets_at(truth$irreducible_modules, tol)
  est_mod_sets <- .family_sets_at(estimates$irreducible_modules, tol)
  est_mod_keys <- names(est_mod_sets)
  joint <- lengths(truth_mod_sets) > 1L
  joint_miss <- NA_real_
  over_frag <- NA_real_
  if (any(joint)) {
    joint_keys <- names(truth_mod_sets)[joint]
    missed <- 0L
    fragmented <- 0L
    for (jk in joint_keys) {
      T_set <- truth_mod_sets[[jk]]
      if (!(jk %in% est_mod_keys)) {
        missed <- missed + 1L
      }
      # fragmentation: two or more estimated modules that are strict
      # subsets of the truth joint module
      n_sub <- sum(vapply(est_mod_sets, function(E) {
        length(E) < length(T_set) && all(E %in% T_set)
      }, logical(1)))
      if (n_sub >= 2L) fragmented <- fragmented + 1L
    }
    joint_miss <- missed / length(joint_keys)
    over_frag <- fragmented / length(joint_keys)
  }

  cand_exact <- cand$exact
  rep_exact <- rep_fam$exact
  mod_exact <- mod_fam$exact
  overall <- as.numeric(isTRUE(cand_exact == 1) && isTRUE(rep_exact == 1) &&
                          isTRUE(mod_exact == 1))
  conditional <- if (isTRUE(cand_exact == 1)) {
    as.numeric(isTRUE(rep_exact == 1) && isTRUE(mod_exact == 1))
  } else {
    NA_real_
  }

  row <- data.frame(
    type1_omnibus = type1,
    power_omnibus = power,
    candidate_tpr = cand$tpr,
    candidate_fdp = cand$fdp,
    candidate_exact_recovery = cand$exact,
    candidate_jaccard = cand$jaccard,
    representative_exact_recovery = rep_fam$exact,
    representative_jaccard = rep_fam$jaccard,
    representative_size_bias = size_bias,
    tie_coverage = tie_coverage,
    irreducible_tpr = mod_fam$tpr,
    irreducible_fdp = mod_fam$fdp,
    irreducible_exact_recovery = mod_fam$exact,
    irreducible_jaccard = mod_fam$jaccard,
    joint_module_miss_rate = joint_miss,
    over_fragmentation_rate = over_frag,
    overall_recovery = overall,
    conditional_recovery = conditional,
    runtime = if (is.null(runtime)) NA_real_ else runtime,
    stringsAsFactors = FALSE
  )
  # scenario descriptors used for grouping
  for (g in c("experiment", "n", "architecture", "locus_pve",
              "correlation", "target_loss", "tolerance")) {
    val <- settings[[g]]
    if (is.null(val) || length(val) != 1L) {
      val <- if (g == "tolerance") tol else NA
    }
    row[[g]] <- if (is.numeric(val)) as.numeric(val) else as.character(val)
  }
  row
}

#' Set/family metrics under the frozen empty-set rules
#'
#' @param truth_keys,est_keys Character vectors of normalised keys
#'   (traits or trait-set keys).
#' @return List with `tpr`, `fdp`, `exact`, `jaccard`.
#' @keywords internal
.set_family_metrics <- function(truth_keys, est_keys) {
  truth_keys <- unique(truth_keys)
  est_keys <- unique(est_keys)
  n_t <- length(truth_keys)
  n_e <- length(est_keys)
  tp <- length(intersect(truth_keys, est_keys))
  if (n_t == 0L && n_e == 0L) {
    return(list(tpr = NA_real_, fdp = 0, exact = 1, jaccard = 1))
  }
  if (n_e == 0L) {
    return(list(tpr = 0, fdp = 0, exact = 0, jaccard = 0))
  }
  if (n_t == 0L) {
    return(list(tpr = NA_real_, fdp = 1, exact = 0, jaccard = 0))
  }
  list(
    tpr = tp / n_t,
    fdp = (n_e - tp) / max(1, n_e),
    exact = as.numeric(n_t == n_e && tp == n_t),
    jaccard = tp / (n_t + n_e - tp)
  )
}

#' Family trait keys at one tolerance
#' @param tab Extraction table (or `NULL`).
#' @param tol Tolerance.
#' @keywords internal
.family_keys_at <- function(tab, tol) {
  if (is.null(tab) || nrow(tab) == 0L) return(character())
  sort(unique(tab$trait_key[tab$tolerance == tol]))
}

#' Family sets (list of trait vectors, named by key) at one tolerance
#' @keywords internal
.family_sets_at <- function(tab, tol) {
  if (is.null(tab) || nrow(tab) == 0L) return(list())
  keep <- tab$tolerance == tol
  sets <- tab$trait_set[keep]
  names(sets) <- tab$trait_key[keep]
  sets
}

#' Minimum set size at one tolerance (Inf when empty)
#' @keywords internal
.min_size_at <- function(tab, tol) {
  if (is.null(tab) || nrow(tab) == 0L) return(Inf)
  keep <- tab$tolerance == tol
  if (!any(keep)) return(Inf)
  min(tab$set_size[keep])
}

#' Grouped means and Monte Carlo standard errors
#'
#' @param df Replicate metrics data.frame.
#' @param metrics Metric columns to summarise.
#' @param group_cols Grouping columns (possibly empty).
#' @param conf_level Confidence level.
#' @return List with `summary` and `mcse` data.frames.
#' @keywords internal
.summarise_metrics <- function(df, metrics, group_cols, conf_level) {
  metrics <- intersect(metrics, names(df))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  if (length(group_cols) == 0L) {
    groups <- list(seq_len(nrow(df)))
    group_keys <- NULL
  } else {
    key <- interaction(
      lapply(df[group_cols], function(x) {
        x <- as.character(x)
        x[is.na(x)] <- "<NA>"
        x
      }),
      drop = TRUE, lex.order = TRUE
    )
    groups <- split(seq_len(nrow(df)), key)
    group_keys <- lapply(groups, function(idx) {
      df[idx[1L], group_cols, drop = FALSE]
    })
  }
  sum_rows <- list()
  mcse_rows <- list()
  for (gi in seq_along(groups)) {
    idx <- groups[[gi]]
    srow <- if (!is.null(group_keys)) {
      as.list(group_keys[[gi]])
    } else {
      list()
    }
    srow$n <- nrow(idx)
    mrow <- srow
    for (met in metrics) {
      x <- df[[met]][idx]
      x <- x[!is.na(x)]
      n <- length(x)
      if (n == 0L) {
        srow[[met]] <- NA_real_
        mrow[[met]] <- NA_real_
        next
      }
      mn <- mean(x)
      is_prop <- all(x %in% c(0, 1))
      se <- if (n > 1L) {
        if (is_prop) sqrt(mn * (1 - mn) / n) else stats::sd(x) / sqrt(n)
      } else {
        NA_real_
      }
      srow[[met]] <- mn
      srow[[paste0(met, "_ci_lo")]] <- if (is.na(se)) NA_real_ else mn - z * se
      srow[[paste0(met, "_ci_hi")]] <- if (is.na(se)) NA_real_ else mn + z * se
      mrow[[met]] <- se
    }
    sum_rows[[gi]] <- as.data.frame(srow, stringsAsFactors = FALSE,
                                    optional = TRUE)
    mcse_rows[[gi]] <- as.data.frame(mrow, stringsAsFactors = FALSE,
                                     optional = TRUE)
  }
  list(
    summary = do.call(rbind, sum_rows),
    mcse = do.call(rbind, mcse_rows)
  )
}
