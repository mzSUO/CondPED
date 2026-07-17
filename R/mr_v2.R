# ==============================================================================
# mr.R — Layer 3: MR-supported directional evidence for CondPED
# ==============================================================================
# Scientific scope
# ----------------
# This file provides directional evidence for trait pairs. It does NOT prove the
# exclusion restriction and does NOT identify a biological mechanism by itself.
#
# Main workflow
# -------------
#   1. select_instruments_pairwise_filter()
#      Exposure relevance + F statistic + LD pruning + pairwise conditional
#      outcome-effect filtering.
#   2. estimate_causal_ld_gls()
#      LD-aware generalized IVW/GLS estimate, optionally with individual-level
#      bootstrap uncertainty.
#   3. steiger_direction_diagnostic(), cochran_q_diagnostic(),
#      weighted_median_mr()
#      Direction and robustness diagnostics.
#   4. bidirectional_mr()
#      Runs both directions for one trait pair.
#   5. run_all_trait_pairs_mr()
#      Runs all unordered trait pairs in an m-trait dataset.
#
# Required marginal-effect format
# -------------------------------
# data.frame with columns:
#   TRAIT, SNPID, A, SE, P_Value
#
# Matrix requirements
# -------------------
#   Y_residual: n x m matrix with trait names as column names.
#   X_all:      n x p matrix with SNP IDs as column names.
# ==============================================================================


# ------------------------------------------------------------------------------
# Validation and numerical helpers
# ------------------------------------------------------------------------------

.validate_marginal_effects <- function(marginal_effects) {
  required <- c("TRAIT", "SNPID", "A", "SE", "P_Value")
  missing_cols <- setdiff(required, names(marginal_effects))
  if (length(missing_cols) > 0L) {
    stop(
      "marginal_effects 缺少必要列: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  out <- marginal_effects[, required, drop = FALSE]
  out$TRAIT <- as.character(out$TRAIT)
  out$SNPID <- as.character(out$SNPID)
  out$A <- suppressWarnings(as.numeric(out$A))
  out$SE <- suppressWarnings(as.numeric(out$SE))
  out$P_Value <- suppressWarnings(as.numeric(out$P_Value))

  duplicated_key <- duplicated(out[c("TRAIT", "SNPID")])
  if (any(duplicated_key)) {
    bad <- out[duplicated_key, c("TRAIT", "SNPID"), drop = FALSE]
    stop(
      "marginal_effects 中存在重复的 TRAIT-SNPID 记录，例如: ",
      paste0(bad$TRAIT[1L], "/", bad$SNPID[1L]),
      call. = FALSE
    )
  }
  out
}

.safe_inverse <- function(M, ridge = 1e-8) {
  M <- as.matrix(M)
  if (nrow(M) != ncol(M)) stop("待求逆矩阵必须为方阵。", call. = FALSE)
  M <- (M + t(M)) / 2

  ans <- tryCatch(solve(M), error = function(e) NULL)
  if (!is.null(ans) && all(is.finite(ans))) {
    return(ans)
  }

  scale_value <- mean(diag(M), na.rm = TRUE)
  if (!is.finite(scale_value) || scale_value <= 0) scale_value <- 1
  M_ridge <- M + diag(ridge * scale_value, nrow(M))

  ans <- tryCatch(solve(M_ridge), error = function(e) NULL)
  if (is.null(ans) || !all(is.finite(ans))) {
    stop("矩阵在 ridge 正则化后仍无法稳定求逆。", call. = FALSE)
  }
  ans
}

compute_f_stat <- function(A, SE) {
  A <- suppressWarnings(as.numeric(A))
  SE <- suppressWarnings(as.numeric(SE))
  out <- (A / SE)^2
  out[!is.finite(out) | SE <= 0] <- NA_real_
  out
}

.fast_lm_effect <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  n <- length(y)

  if (n < 4L || stats::var(x) < 1e-12) {
    return(c(beta = NA_real_, se = NA_real_, pval = NA_real_))
  }

  xc <- x - mean(x)
  yc <- y - mean(y)
  sxx <- sum(xc^2)
  beta <- sum(xc * yc) / sxx
  intercept <- mean(y) - beta * mean(x)
  residual <- y - intercept - beta * x
  sigma2 <- sum(residual^2) / max(n - 2L, 1L)
  se <- sqrt(sigma2 / sxx)
  t_value <- beta / se
  pval <- 2 * stats::pt(-abs(t_value), df = n - 2L)

  c(beta = beta, se = se, pval = pval)
}

.project_outcome_pairwise <- function(Y_residual, exposure, outcome, V = NULL) {
  Y_residual <- as.matrix(Y_residual)
  if (is.null(colnames(Y_residual))) {
    stop("Y_residual 必须具有性状列名。", call. = FALSE)
  }
  if (!all(c(exposure, outcome) %in% colnames(Y_residual))) {
    stop("exposure 或 outcome 不在 Y_residual 列名中。", call. = FALSE)
  }
  if (identical(exposure, outcome)) {
    stop("exposure 与 outcome 不能相同。", call. = FALSE)
  }

  if (is.null(V)) V <- stats::cov(Y_residual, use = "pairwise.complete.obs")
  ie <- match(exposure, colnames(Y_residual))
  io <- match(outcome, colnames(Y_residual))

  var_exp <- V[ie, ie]
  cov_out_exp <- V[io, ie]
  if (!is.finite(var_exp) || var_exp <= 1e-12) {
    stop("暴露性状方差过小，无法构造成对条件结局。", call. = FALSE)
  }

  gamma <- cov_out_exp / var_exp
  y_cond <- Y_residual[, io] - gamma * Y_residual[, ie]

  list(
    y_cond = as.numeric(y_cond),
    gamma = as.numeric(gamma),
    exposure = exposure,
    outcome = outcome
  )
}

.project_outcome_full <- function(Y_residual, outcome, V = NULL) {
  Y_residual <- as.matrix(Y_residual)
  if (is.null(colnames(Y_residual))) {
    stop("Y_residual 必须具有性状列名。", call. = FALSE)
  }
  if (!outcome %in% colnames(Y_residual)) {
    stop("outcome 不在 Y_residual 列名中。", call. = FALSE)
  }
  if (ncol(Y_residual) < 2L) {
    stop("全条件投影至少需要两个性状。", call. = FALSE)
  }

  if (is.null(V)) V <- stats::cov(Y_residual, use = "pairwise.complete.obs")
  io <- match(outcome, colnames(Y_residual))
  cond_idx <- setdiff(seq_len(ncol(Y_residual)), io)
  V_cond <- V[cond_idx, cond_idx, drop = FALSE]
  C_cond_out <- V[cond_idx, io, drop = FALSE]
  gamma <- .safe_inverse(V_cond) %*% C_cond_out
  y_cond <- Y_residual[, io] - Y_residual[, cond_idx, drop = FALSE] %*% gamma

  list(
    y_cond = as.numeric(y_cond),
    gamma = as.numeric(gamma),
    exposure = paste(colnames(Y_residual)[cond_idx], collapse = "+"),
    outcome = outcome
  )
}


# ------------------------------------------------------------------------------
# LD pruning
# ------------------------------------------------------------------------------

ld_prune_ivs <- function(loci,
                         ld_matrix,
                         r2_threshold = 0.1,
                         priority = NULL) {
  loci <- unique(as.character(loci))
  if (length(loci) <= 1L) {
    return(loci)
  }
  if (is.null(rownames(ld_matrix)) || is.null(colnames(ld_matrix))) {
    stop("ld_matrix 必须具有 SNP 行名和列名。", call. = FALSE)
  }

  shared <- intersect(loci, intersect(rownames(ld_matrix), colnames(ld_matrix)))
  if (length(shared) <= 1L) {
    return(shared)
  }

  if (is.null(priority)) {
    ordered <- shared
  } else {
    priority <- priority[shared]
    priority[!is.finite(priority)] <- -Inf
    ordered <- names(sort(priority, decreasing = TRUE))
  }

  keep <- character(0L)
  for (locus in ordered) {
    if (length(keep) == 0L) {
      keep <- locus
      next
    }
    r_values <- as.numeric(ld_matrix[locus, keep, drop = TRUE])
    max_r2 <- max(r_values^2, na.rm = TRUE)
    if (!is.finite(max_r2) || max_r2 < r2_threshold) {
      keep <- c(keep, locus)
    }
  }
  keep
}


# ------------------------------------------------------------------------------
# Step 1: direction-specific IV selection
# ------------------------------------------------------------------------------

#' Select candidate instruments with pairwise conditional outcome filtering
#'
#' The conditional filter is an empirical pleiotropy-reduction step. It does not
#' prove or guarantee the MR exclusion restriction.
#'
#' @param marginal_effects data.frame(TRAIT, SNPID, A, SE, P_Value).
#' @param Y_residual n x m adjusted phenotype matrix.
#' @param X_loci n x p genotype matrix with SNP IDs as column names.
#' @param exposure Exposure trait name.
#' @param outcome Outcome trait name.
#' @param alpha1 Marginal relevance threshold.
#' @param alpha_filter Pairwise conditional outcome-effect filter threshold.
#' @param F_threshold Minimum SNP-specific F statistic.
#' @param r2_threshold LD-pruning threshold on r^2.
#' @param ld_matrix Optional SNP correlation matrix R (not R^2).
#' @param V Optional phenotypic covariance matrix.
#' @param screen_mode "pairwise", "full", or "none". "full" is retained only
#'   for method comparison; "pairwise" is the recommended CondPED mode.
#'
#' @return Structured list with ivs, status, and a per-candidate audit table.
#' @export
select_instruments_pairwise_filter <- function(
    marginal_effects,
    Y_residual,
    X_loci,
    exposure,
    outcome,
    alpha1,
    alpha_filter = 0.05,
    F_threshold = 10,
    r2_threshold = 0.1,
    ld_matrix = NULL,
    V = NULL,
    screen_mode = c("pairwise", "full", "none"),
    min_iv = 3L) {
  screen_mode <- match.arg(screen_mode)
  marginal_effects <- .validate_marginal_effects(marginal_effects)
  Y_residual <- as.matrix(Y_residual)
  X_loci <- as.matrix(X_loci)

  if (is.null(colnames(Y_residual))) stop("Y_residual 必须具有性状列名。", call. = FALSE)
  if (is.null(colnames(X_loci))) stop("X_loci 必须具有 SNP 列名。", call. = FALSE)
  if (!all(c(exposure, outcome) %in% colnames(Y_residual))) {
    stop("exposure/outcome 未在 Y_residual 中找到。", call. = FALSE)
  }
  if (!is.numeric(alpha1) || length(alpha1) != 1L || !is.finite(alpha1) || alpha1 <= 0 || alpha1 >= 1) {
    stop("alpha1 必须为 (0,1) 内的单个数值。", call. = FALSE)
  }
  if (!is.numeric(alpha_filter) || length(alpha_filter) != 1L || !is.finite(alpha_filter) || alpha_filter <= 0 || alpha_filter >= 1) {
    stop("alpha_filter 必须为 (0,1) 内的单个数值。", call. = FALSE)
  }

  exp_df <- marginal_effects[marginal_effects$TRAIT == exposure, , drop = FALSE]
  exp_df <- exp_df[exp_df$SNPID %in% colnames(X_loci), , drop = FALSE]
  exp_df$relevance_pass <- !is.na(exp_df$P_Value) & exp_df$P_Value < alpha1
  exp_df$F_stat <- compute_f_stat(exp_df$A, exp_df$SE)
  exp_df$strength_pass <- !is.na(exp_df$F_stat) & exp_df$F_stat > F_threshold

  audit <- exp_df[exp_df$relevance_pass, c("SNPID", "A", "SE", "P_Value", "F_stat", "strength_pass"), drop = FALSE]
  names(audit)[names(audit) == "A"] <- "beta_exposure"
  names(audit)[names(audit) == "SE"] <- "se_exposure"
  names(audit)[names(audit) == "P_Value"] <- "p_exposure"

  if (nrow(audit) == 0L) {
    return(list(
      ivs = character(0L), n_iv = 0L, status = "no_relevant_instruments",
      exposure = exposure, outcome = outcome, screen_mode = screen_mode,
      audit = audit, projection = NULL
    ))
  }

  projection <- NULL
  if (screen_mode == "pairwise") {
    projection <- .project_outcome_pairwise(Y_residual, exposure, outcome, V)
  } else if (screen_mode == "full") {
    projection <- .project_outcome_full(Y_residual, outcome, V)
  }

  audit$p_cond_outcome <- NA_real_
  audit$conditional_filter_pass <- TRUE

  if (screen_mode != "none") {
    for (i in seq_len(nrow(audit))) {
      snp <- audit$SNPID[i]
      fit <- .fast_lm_effect(X_loci[, snp], projection$y_cond)
      audit$p_cond_outcome[i] <- fit["pval"]
      audit$conditional_filter_pass[i] <- is.finite(fit["pval"]) && fit["pval"] > alpha_filter
    }
  }

  candidates <- audit$SNPID[
    audit$strength_pass & audit$conditional_filter_pass
  ]

  if (length(candidates) == 0L) {
    audit$ld_pass <- FALSE
    audit$selected <- FALSE
    return(list(
      ivs = character(0L), n_iv = 0L, status = "no_instruments_after_filtering",
      exposure = exposure, outcome = outcome, screen_mode = screen_mode,
      audit = audit, projection = projection
    ))
  }

  if (!is.null(ld_matrix) && length(candidates) > 1L) {
    missing_ld <- setdiff(candidates, intersect(rownames(ld_matrix), colnames(ld_matrix)))
    if (length(missing_ld) > 0L) {
      stop(
        "ld_matrix 缺少候选 SNP: ", paste(missing_ld, collapse = ", "),
        call. = FALSE
      )
    }
    priority <- setNames(audit$F_stat, audit$SNPID)
    selected <- ld_prune_ivs(
      loci = candidates,
      ld_matrix = ld_matrix,
      r2_threshold = r2_threshold,
      priority = priority
    )
  } else {
    selected <- candidates
  }

  audit$ld_pass <- audit$SNPID %in% selected
  audit$selected <- audit$SNPID %in% selected
  status <- if (length(selected) < min_iv) "insufficient_iv" else "ok"

  list(
    ivs = selected,
    n_iv = length(selected),
    status = status,
    exposure = exposure,
    outcome = outcome,
    screen_mode = screen_mode,
    audit = audit,
    projection = projection,
    thresholds = list(
      alpha1 = alpha1,
      alpha_filter = alpha_filter,
      F_threshold = F_threshold,
      r2_threshold = r2_threshold,
      min_iv = min_iv
    )
  )
}

# Backward-compatible alias. The return value is now a structured list.
#' @export
select_instruments_pairwise <- select_instruments_pairwise_filter


# ------------------------------------------------------------------------------
# Step 2: LD-aware causal estimate and diagnostics
# ------------------------------------------------------------------------------

.build_outcome_covariance <- function(se_out, R = NULL, ridge = 1e-8) {
  se_out <- as.numeric(se_out)
  k <- length(se_out)
  if (any(!is.finite(se_out)) || any(se_out <= 0)) {
    stop("se_out 必须全部为有限正数。", call. = FALSE)
  }
  if (is.null(R)) R <- diag(k)
  R <- as.matrix(R)
  if (!all(dim(R) == c(k, k))) stop("R 维度与工具变量数不一致。", call. = FALSE)

  D <- diag(se_out, k)
  Omega <- D %*% R %*% D
  Omega <- (Omega + t(Omega)) / 2
  scale_value <- mean(diag(Omega))
  if (!is.finite(scale_value) || scale_value <= 0) scale_value <- 1
  Omega + diag(ridge * scale_value, k)
}

.estimate_gls_core <- function(beta_exp, beta_out, se_out, R = NULL) {
  beta_exp <- as.numeric(beta_exp)
  beta_out <- as.numeric(beta_out)
  se_out <- as.numeric(se_out)
  k <- length(beta_exp)

  Omega <- .build_outcome_covariance(se_out, R)
  W <- .safe_inverse(Omega)
  denominator <- as.numeric(crossprod(beta_exp, W %*% beta_exp))
  if (!is.finite(denominator) || denominator <= 1e-14) {
    stop("GLS 分母过小，无法稳定估计因果效应。", call. = FALSE)
  }
  numerator <- as.numeric(crossprod(beta_exp, W %*% beta_out))
  gamma <- numerator / denominator
  se <- sqrt(1 / denominator)

  list(gamma = gamma, se = se, W = W, Omega = Omega)
}

#' Cochran Q heterogeneity diagnostic
#' @export
cochran_q_diagnostic <- function(beta_exp, beta_out, se_out, gamma, R = NULL) {
  k <- length(beta_exp)
  if (k < 3L || any(!is.finite(c(beta_exp, beta_out, se_out, gamma)))) {
    return(list(Q = NA_real_, df = NA_integer_, pval = NA_real_, heterogeneous = NA))
  }
  Omega <- .build_outcome_covariance(se_out, R)
  W <- .safe_inverse(Omega)
  residual <- as.numeric(beta_out - gamma * beta_exp)
  Q <- as.numeric(crossprod(residual, W %*% residual))
  df <- k - 1L
  pval <- stats::pchisq(Q, df = df, lower.tail = FALSE)
  list(Q = Q, df = df, pval = pval, heterogeneous = is.finite(pval) && pval < 0.05)
}

.weighted_median <- function(values, weights) {
  ok <- is.finite(values) & is.finite(weights) & weights > 0
  values <- values[ok]
  weights <- weights[ok]
  if (length(values) == 0L) {
    return(NA_real_)
  }
  ord <- order(values)
  values <- values[ord]
  weights <- weights[ord] / sum(weights)
  values[which(cumsum(weights) >= 0.5)[1L]]
}

#' Weighted-median MR sensitivity estimate
#' @export
weighted_median_mr <- function(beta_exp, beta_out, se_out) {
  ok <- is.finite(beta_exp) & is.finite(beta_out) & is.finite(se_out) &
    abs(beta_exp) > 1e-12 & se_out > 0
  if (sum(ok) < 3L) {
    return(list(gamma = NA_real_, n_iv = sum(ok), status = "insufficient_iv"))
  }
  ratio <- beta_out[ok] / beta_exp[ok]
  weights <- beta_exp[ok]^2 / se_out[ok]^2
  list(
    gamma = .weighted_median(ratio, weights),
    n_iv = sum(ok),
    status = "ok"
  )
}

#' Steiger direction-consistency diagnostic
#'
#' Compares variance explained in exposure and outcome for each IV. This is a
#' direction diagnostic, not an exclusion-restriction test.
#'
#' @export
steiger_direction_diagnostic <- function(beta_exp,
                                         beta_out,
                                         var_g,
                                         var_exp,
                                         var_out,
                                         consistency_threshold = 0.5) {
  beta_exp <- as.numeric(beta_exp)
  beta_out <- as.numeric(beta_out)
  var_g <- as.numeric(var_g)

  ok <- is.finite(beta_exp) & is.finite(beta_out) & is.finite(var_g) & var_g > 0
  if (sum(ok) == 0L || !is.finite(var_exp) || !is.finite(var_out) ||
    var_exp <= 0 || var_out <= 0) {
    return(list(
      consistent = NA,
      proportion_consistent = NA_real_,
      n_iv = 0L,
      per_iv = data.frame()
    ))
  }

  r2_exp <- var_g[ok] * beta_exp[ok]^2 / var_exp
  r2_out <- var_g[ok] * beta_out[ok]^2 / var_out
  per_iv_consistent <- r2_exp > r2_out
  proportion <- mean(per_iv_consistent)

  list(
    consistent = is.finite(proportion) && proportion >= consistency_threshold,
    proportion_consistent = proportion,
    n_iv = sum(ok),
    per_iv = data.frame(
      r2_exposure = r2_exp,
      r2_outcome = r2_out,
      consistent = per_iv_consistent
    )
  )
}

#' LD-aware generalized IVW/GLS causal estimate
#'
#' @param theta_exp Named exposure-effect vector.
#' @param theta_out Named outcome-effect vector.
#' @param se_out Named outcome-effect SE vector.
#' @param ld_matrix Optional SNP correlation matrix R.
#' @param Y Optional individual-level phenotype matrix for bootstrap.
#' @param X_iv Optional individual-level IV genotype matrix for bootstrap.
#' @param exp_col,out_col Columns in Y.
#' @param n_boot Number of individual-level bootstrap replicates. Set 0 to use
#'   the model-based GLS SE only.
#'
#' @return Causal estimate with Q and weighted-median diagnostics.
#' @export
estimate_causal_ld_gls <- function(theta_exp,
                                   theta_out,
                                   se_out,
                                   ld_matrix = NULL,
                                   Y = NULL,
                                   X_iv = NULL,
                                   exp_col = 1L,
                                   out_col = 2L,
                                   n_boot = 200L,
                                   seed = NULL) {
  shared <- Reduce(intersect, list(names(theta_exp), names(theta_out), names(se_out)))
  if (length(shared) < 3L) stop("至少需要 3 个共有 IV。", call. = FALSE)

  theta_exp <- as.numeric(theta_exp[shared])
  theta_out <- as.numeric(theta_out[shared])
  se_out <- as.numeric(se_out[shared])
  names(theta_exp) <- names(theta_out) <- names(se_out) <- shared

  ok <- is.finite(theta_exp) & is.finite(theta_out) & is.finite(se_out) & se_out > 0
  theta_exp <- theta_exp[ok]
  theta_out <- theta_out[ok]
  se_out <- se_out[ok]
  shared <- names(theta_exp)
  k <- length(shared)
  if (k < 3L) stop("去除缺失值后少于 3 个 IV。", call. = FALSE)

  R <- NULL
  if (!is.null(ld_matrix)) {
    if (is.null(rownames(ld_matrix)) || !all(shared %in% rownames(ld_matrix))) {
      stop("ld_matrix 不包含全部 IV 的行名。", call. = FALSE)
    }
    R <- ld_matrix[shared, shared, drop = FALSE]
  }

  core <- .estimate_gls_core(theta_exp, theta_out, se_out, R)
  gamma_hat <- core$gamma
  se_model <- core$se
  se_final <- se_model
  method <- "model_based_gls"
  bootstrap_values <- numeric(0L)

  can_bootstrap <- n_boot > 1L && !is.null(Y) && !is.null(X_iv)
  if (can_bootstrap) {
    Y <- as.matrix(Y)
    X_iv <- as.matrix(X_iv)
    if (nrow(Y) != nrow(X_iv)) stop("Y 与 X_iv 行数不一致。", call. = FALSE)
    if (!all(shared %in% colnames(X_iv))) stop("X_iv 缺少部分 IV 列。", call. = FALSE)
    if (max(exp_col, out_col) > ncol(Y)) stop("exp_col/out_col 超出 Y 列数。", call. = FALSE)
    X_iv <- X_iv[, shared, drop = FALSE]

    if (!is.null(seed)) set.seed(seed)
    n <- nrow(Y)
    bootstrap_values <- rep(NA_real_, n_boot)

    for (b in seq_len(n_boot)) {
      idx <- sample.int(n, n, replace = TRUE)
      Xb <- X_iv[idx, , drop = FALSE]
      y_exp <- Y[idx, exp_col]
      y_out <- Y[idx, out_col]

      eff_exp <- t(vapply(seq_len(k), function(j) {
        .fast_lm_effect(Xb[, j], y_exp)
      }, numeric(3L)))
      eff_out <- t(vapply(seq_len(k), function(j) {
        .fast_lm_effect(Xb[, j], y_out)
      }, numeric(3L)))

      bx <- eff_exp[, "beta"]
      by <- eff_out[, "beta"]
      sy <- eff_out[, "se"]
      ok_b <- is.finite(bx) & is.finite(by) & is.finite(sy) & sy > 0
      if (sum(ok_b) < 3L) next

      R_b <- NULL
      if (!is.null(R)) {
        R_b <- stats::cor(Xb[, ok_b, drop = FALSE], use = "pairwise.complete.obs")
        if (sum(ok_b) == 1L) R_b <- matrix(1, 1, 1)
      }

      bootstrap_values[b] <- tryCatch(
        .estimate_gls_core(bx[ok_b], by[ok_b], sy[ok_b], R_b)$gamma,
        error = function(e) NA_real_
      )
    }

    valid_boot <- bootstrap_values[is.finite(bootstrap_values)]
    if (length(valid_boot) >= max(30L, ceiling(0.5 * n_boot))) {
      se_boot <- stats::sd(valid_boot)
      if (is.finite(se_boot) && se_boot > 0) {
        se_final <- se_boot
        method <- "individual_bootstrap_gls"
      }
    }
  }

  z <- gamma_hat / se_final
  pval <- 2 * stats::pnorm(-abs(z))
  q_diag <- cochran_q_diagnostic(theta_exp, theta_out, se_out, gamma_hat, R)
  wm <- weighted_median_mr(theta_exp, theta_out, se_out)

  list(
    gamma = gamma_hat,
    se = se_final,
    se_model = se_model,
    z = z,
    pval = pval,
    n_iv = k,
    method = method,
    ivs = shared,
    q = q_diag,
    weighted_median = wm,
    bootstrap_values = bootstrap_values,
    status = "ok"
  )
}

# Backward-compatible wrapper. New code should call estimate_causal_ld_gls().
#' @export
estimate_causal_gls <- function(theta_exp,
                                theta_out,
                                ld_matrix = NULL,
                                Y = NULL,
                                X_iv = NULL,
                                exp_col = 1L,
                                out_col = 2L,
                                n_boot = 200L,
                                se_out = NULL,
                                seed = NULL) {
  if (is.null(se_out)) {
    warning(
      "estimate_causal_gls() 未提供 se_out；将使用等权近似。建议改用 estimate_causal_ld_gls()。",
      call. = FALSE
    )
    se_out <- setNames(rep(1, length(theta_out)), names(theta_out))
  }
  estimate_causal_ld_gls(
    theta_exp = theta_exp,
    theta_out = theta_out,
    se_out = se_out,
    ld_matrix = ld_matrix,
    Y = Y,
    X_iv = X_iv,
    exp_col = exp_col,
    out_col = out_col,
    n_boot = n_boot,
    seed = seed
  )
}


# ------------------------------------------------------------------------------
# Step 3: one-direction and bidirectional MR
# ------------------------------------------------------------------------------

.empty_direction_result <- function(exposure, outcome, status, selection = NULL) {
  list(
    exposure = exposure,
    outcome = outcome,
    gamma = NA_real_,
    se = NA_real_,
    se_model = NA_real_,
    z = NA_real_,
    pval = NA_real_,
    sig = FALSE,
    n_iv = if (is.null(selection)) 0L else selection$n_iv,
    method = NA_character_,
    ivs = if (is.null(selection)) character(0L) else selection$ivs,
    status = status,
    selection = selection,
    steiger = list(consistent = NA, proportion_consistent = NA_real_, n_iv = 0L),
    q = list(Q = NA_real_, df = NA_integer_, pval = NA_real_, heterogeneous = NA),
    weighted_median = list(gamma = NA_real_, n_iv = 0L, status = "not_run"),
    robust_direction_consistent = NA
  )
}

.run_mr_direction <- function(marginal_effects,
                              exposure,
                              outcome,
                              Y_residual,
                              Y,
                              X_all,
                              ld_matrix,
                              alpha1,
                              alpha_filter,
                              F_threshold,
                              r2_threshold,
                              alpha_mr,
                              n_boot,
                              V,
                              screen_mode,
                              min_iv,
                              steiger_threshold,
                              robust_tolerance,
                              seed = NULL) {
  selection <- select_instruments_pairwise_filter(
    marginal_effects = marginal_effects,
    Y_residual = Y_residual,
    X_loci = X_all,
    exposure = exposure,
    outcome = outcome,
    alpha1 = alpha1,
    alpha_filter = alpha_filter,
    F_threshold = F_threshold,
    r2_threshold = r2_threshold,
    ld_matrix = ld_matrix,
    V = V,
    screen_mode = screen_mode,
    min_iv = min_iv
  )

  if (!identical(selection$status, "ok")) {
    return(.empty_direction_result(exposure, outcome, selection$status, selection))
  }

  me <- .validate_marginal_effects(marginal_effects)
  ivs <- selection$ivs
  exp_df <- me[me$TRAIT == exposure & me$SNPID %in% ivs, , drop = FALSE]
  out_df <- me[me$TRAIT == outcome & me$SNPID %in% ivs, , drop = FALSE]

  beta_exp <- setNames(exp_df$A, exp_df$SNPID)[ivs]
  beta_out <- setNames(out_df$A, out_df$SNPID)[ivs]
  se_out <- setNames(out_df$SE, out_df$SNPID)[ivs]
  complete <- is.finite(beta_exp) & is.finite(beta_out) & is.finite(se_out) & se_out > 0
  ivs_used <- ivs[complete]

  if (length(ivs_used) < min_iv) {
    selection$ivs <- ivs_used
    selection$n_iv <- length(ivs_used)
    selection$status <- "insufficient_complete_effects"
    return(.empty_direction_result(exposure, outcome, selection$status, selection))
  }

  beta_exp <- beta_exp[ivs_used]
  beta_out <- beta_out[ivs_used]
  se_out <- se_out[ivs_used]

  trait_names <- colnames(Y_residual)
  exp_col <- match(exposure, trait_names)
  out_col <- match(outcome, trait_names)
  Y_for_boot <- if (is.null(Y)) NULL else as.matrix(Y)[, c(exp_col, out_col), drop = FALSE]
  X_iv <- X_all[, ivs_used, drop = FALSE]

  estimate <- tryCatch(
    estimate_causal_ld_gls(
      theta_exp = beta_exp,
      theta_out = beta_out,
      se_out = se_out,
      ld_matrix = ld_matrix,
      Y = Y_for_boot,
      X_iv = X_iv,
      exp_col = 1L,
      out_col = 2L,
      n_boot = n_boot,
      seed = seed
    ),
    error = function(e) e
  )

  if (inherits(estimate, "error")) {
    result <- .empty_direction_result(exposure, outcome, "estimation_failed", selection)
    result$error_message <- conditionMessage(estimate)
    return(result)
  }

  var_g <- apply(X_iv, 2L, stats::var, na.rm = TRUE)
  steiger <- steiger_direction_diagnostic(
    beta_exp = beta_exp,
    beta_out = beta_out,
    var_g = var_g,
    var_exp = stats::var(Y_residual[, exposure], na.rm = TRUE),
    var_out = stats::var(Y_residual[, outcome], na.rm = TRUE),
    consistency_threshold = steiger_threshold
  )

  wm_gamma <- estimate$weighted_median$gamma
  robust_direction_consistent <- if (is.finite(wm_gamma) && is.finite(estimate$gamma)) {
    same_sign <- sign(wm_gamma) == sign(estimate$gamma)
    relative_difference <- abs(wm_gamma - estimate$gamma) / max(abs(estimate$gamma), 1e-8)
    same_sign && relative_difference <= robust_tolerance
  } else {
    NA
  }

  c(
    list(exposure = exposure, outcome = outcome),
    estimate,
    list(
      sig = is.finite(estimate$pval) && estimate$pval < alpha_mr,
      selection = selection,
      steiger = steiger,
      robust_direction_consistent = robust_direction_consistent,
      alpha_mr = alpha_mr
    )
  )
}

#' Bidirectional MR-supported directional evidence for one trait pair
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{forward}: MR result for the first directional test.
#'   \item \code{reverse}: MR result for the reverse directional test.
#' }
#' #' @export
bidirectional_mr <- function(marginal_effects,
                             traits,
                             Y_residual,
                             Y = NULL,
                             X_all,
                             ld_matrix = NULL,
                             alpha1,
                             alpha_filter = 0.05,
                             F_threshold = 10,
                             r2_threshold = 0.1,
                             alpha_mr = 0.05,
                             n_boot = 200L,
                             V = NULL,
                             screen_mode = c("pairwise", "full", "none"),
                             min_iv = 3L,
                             steiger_threshold = 0.5,
                             robust_tolerance = 1.0,
                             seed = NULL) {
  screen_mode <- match.arg(screen_mode)
  if (length(traits) != 2L || anyDuplicated(traits)) {
    stop("traits 必须是两个不同性状名称。", call. = FALSE)
  }
  Y_residual <- as.matrix(Y_residual)
  X_all <- as.matrix(X_all)
  if (is.null(colnames(Y_residual))) stop("Y_residual 必须具有性状列名。", call. = FALSE)
  if (is.null(colnames(X_all))) stop("X_all 必须具有 SNP 列名。", call. = FALSE)
  if (!all(traits %in% colnames(Y_residual))) stop("traits 不在 Y_residual 中。", call. = FALSE)
  if (nrow(Y_residual) != nrow(X_all)) stop("Y_residual 与 X_all 行数不一致。", call. = FALSE)
  if (!is.null(Y) && nrow(as.matrix(Y)) != nrow(X_all)) stop("Y 与 X_all 行数不一致。", call. = FALSE)
  if (is.null(V)) V <- stats::cov(Y_residual, use = "pairwise.complete.obs")

  trA <- traits[1L]
  trB <- traits[2L]
  seed_AB <- if (is.null(seed)) NULL else as.integer(seed)
  seed_BA <- if (is.null(seed)) NULL else as.integer(seed) + 1L

  AB <- .run_mr_direction(
    marginal_effects, trA, trB, Y_residual, Y, X_all, ld_matrix,
    alpha1, alpha_filter, F_threshold, r2_threshold, alpha_mr,
    n_boot, V, screen_mode, min_iv, steiger_threshold,
    robust_tolerance, seed_AB
  )
  BA <- .run_mr_direction(
    marginal_effects, trB, trA, Y_residual, Y, X_all, ld_matrix,
    alpha1, alpha_filter, F_threshold, r2_threshold, alpha_mr,
    n_boot, V, screen_mode, min_iv, steiger_threshold,
    robust_tolerance, seed_BA
  )

  list(
    traits = traits,
    AB = AB,
    BA = BA,
    iv_AB = AB$ivs,
    iv_BA = BA$ivs,
    settings = list(
      alpha1 = alpha1,
      alpha_filter = alpha_filter,
      F_threshold = F_threshold,
      r2_threshold = r2_threshold,
      alpha_mr = alpha_mr,
      n_boot = n_boot,
      screen_mode = screen_mode,
      min_iv = min_iv,
      steiger_threshold = steiger_threshold,
      robust_tolerance = robust_tolerance
    )
  )
}

#' Run bidirectional MR for all unordered trait pairs
#'
#' For m traits, this evaluates choose(m, 2) pairs and m(m-1) directions.
#' @export
run_all_trait_pairs_mr <- function(marginal_effects,
                                   traits,
                                   Y_residual,
                                   Y = NULL,
                                   X_all,
                                   ld_matrix = NULL,
                                   alpha1,
                                   alpha_filter = 0.05,
                                   F_threshold = 10,
                                   r2_threshold = 0.1,
                                   alpha_mr = NULL,
                                   adjust_mr = TRUE,
                                   n_boot = 200L,
                                   V = NULL,
                                   screen_mode = c("pairwise", "full", "none"),
                                   min_iv = 3L,
                                   steiger_threshold = 0.5,
                                   robust_tolerance = 1.0,
                                   seed = NULL) {
  screen_mode <- match.arg(screen_mode)
  traits <- unique(as.character(traits))
  if (length(traits) < 2L) stop("至少需要两个性状。", call. = FALSE)
  if (is.null(alpha_mr)) {
    alpha_mr <- if (isTRUE(adjust_mr)) 0.05 / (length(traits) * (length(traits) - 1L)) else 0.05
  }
  if (is.null(V)) V <- stats::cov(as.matrix(Y_residual), use = "pairwise.complete.obs")

  pairs <- utils::combn(traits, 2L, simplify = FALSE)
  out <- vector("list", length(pairs))

  for (i in seq_along(pairs)) {
    pair_seed <- if (is.null(seed)) NULL else as.integer(seed) + 2L * (i - 1L)
    out[[i]] <- bidirectional_mr(
      marginal_effects = marginal_effects,
      traits = pairs[[i]],
      Y_residual = Y_residual,
      Y = Y,
      X_all = X_all,
      ld_matrix = ld_matrix,
      alpha1 = alpha1,
      alpha_filter = alpha_filter,
      F_threshold = F_threshold,
      r2_threshold = r2_threshold,
      alpha_mr = alpha_mr,
      n_boot = n_boot,
      V = V,
      screen_mode = screen_mode,
      min_iv = min_iv,
      steiger_threshold = steiger_threshold,
      robust_tolerance = robust_tolerance,
      seed = pair_seed
    )
  }

  names(out) <- vapply(pairs, function(x) paste(x, collapse = "__"), character(1L))
  attr(out, "alpha_mr") <- alpha_mr
  attr(out, "traits") <- traits
  class(out) <- c("condped_pairwise_mr_list", "list")
  out
}


# ------------------------------------------------------------------------------
# Optional MVMR extension
# ------------------------------------------------------------------------------

#' Multivariable MR using generalized least squares
#'
#' This is an optional extension when co-exposures are pre-specified from a
#' defensible causal model. It is not used by the default CondPED pipeline.
#' @export
mvmr_estimate <- function(marginal_effects,
                          exposures,
                          outcome,
                          alpha1,
                          F_threshold = 10,
                          r2_threshold = 0.1,
                          ld_matrix = NULL,
                          focal = NULL) {
  me <- .validate_marginal_effects(marginal_effects)
  exposures <- unique(as.character(exposures))
  if (length(exposures) < 2L) stop("MVMR 至少需要两个暴露。", call. = FALSE)
  if (outcome %in% exposures) stop("outcome 不能同时属于 exposures。", call. = FALSE)
  if (is.null(focal)) focal <- exposures[1L]
  if (!focal %in% exposures) stop("focal 必须属于 exposures。", call. = FALSE)

  pull <- function(trait) {
    rows <- me[me$TRAIT == trait, , drop = FALSE]
    list(
      beta = setNames(rows$A, rows$SNPID),
      se = setNames(rows$SE, rows$SNPID),
      p = setNames(rows$P_Value, rows$SNPID)
    )
  }

  inst <- character(0L)
  for (e in exposures) {
    pe <- pull(e)
    inst <- union(inst, names(pe$p)[is.finite(pe$p) & pe$p < alpha1])
  }
  if (length(inst) < length(exposures) + 1L) {
    return(list(status = "insufficient_iv", gamma = NA_real_, se = NA_real_, pval = NA_real_, n_iv = length(inst)))
  }

  strong <- setNames(rep(FALSE, length(inst)), inst)
  max_f <- setNames(rep(-Inf, length(inst)), inst)
  for (e in exposures) {
    pe <- pull(e)
    f <- compute_f_stat(pe$beta[inst], pe$se[inst])
    strong <- strong | (!is.na(f) & f > F_threshold)
    f[!is.finite(f)] <- -Inf
    max_f <- pmax(max_f, f)
  }
  inst <- inst[strong]

  if (!is.null(ld_matrix) && length(inst) > 1L) {
    inst <- ld_prune_ivs(inst, ld_matrix, r2_threshold, priority = max_f)
  }
  if (length(inst) < length(exposures) + 1L) {
    return(list(status = "insufficient_iv", gamma = NA_real_, se = NA_real_, pval = NA_real_, n_iv = length(inst)))
  }

  Theta <- vapply(exposures, function(e) pull(e)$beta[inst], numeric(length(inst)))
  colnames(Theta) <- exposures
  beta_out <- pull(outcome)$beta[inst]
  se_out <- pull(outcome)$se[inst]
  ok <- stats::complete.cases(Theta, beta_out, se_out) & se_out > 0
  Theta <- Theta[ok, , drop = FALSE]
  beta_out <- beta_out[ok]
  se_out <- se_out[ok]
  inst <- inst[ok]
  k <- nrow(Theta)

  if (k < length(exposures) + 1L) {
    return(list(status = "insufficient_complete_effects", gamma = NA_real_, se = NA_real_, pval = NA_real_, n_iv = k))
  }

  R <- if (is.null(ld_matrix)) diag(k) else ld_matrix[inst, inst, drop = FALSE]
  Omega <- .build_outcome_covariance(se_out, R)
  W <- .safe_inverse(Omega)
  XtWX <- crossprod(Theta, W %*% Theta)
  XtWy <- crossprod(Theta, W %*% beta_out)
  cov_gamma <- .safe_inverse(XtWX)
  gamma_all <- as.numeric(cov_gamma %*% XtWy)
  names(gamma_all) <- exposures
  se_all <- sqrt(diag(cov_gamma))
  names(se_all) <- exposures

  gamma_focal <- gamma_all[focal]
  se_focal <- se_all[focal]
  z <- gamma_focal / se_focal

  list(
    status = "ok",
    gamma = as.numeric(gamma_focal),
    se = as.numeric(se_focal),
    z = as.numeric(z),
    pval = 2 * stats::pnorm(-abs(z)),
    n_iv = k,
    ivs = inst,
    gamma_all = gamma_all,
    se_all = se_all,
    exposures = exposures,
    focal = focal,
    outcome = outcome
  )
}
