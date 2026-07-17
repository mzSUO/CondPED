# ==============================================================================
# mr.R — Layer 3: Mendelian Randomization for Causal Inference
# ==============================================================================
# Implements Section 2.3.1c + 2.3.2 of the CondPED paper:
#
#   Step 1  select_instruments()          IV selection (four criteria)
#   Step 2  estimate_causal_gls()         GLS causal effect estimation (Eq. 7)
#   Step 3  bidirectional_mr()            Bidirectional MR + significance test
#   Step 4  test_mediation_consistency()  Reverse conditional consistency (Sec. 2.3.2)
#
# Input/Output interface:
#
#   select_instruments(marginal_effects, cond_effects, exposure, outcome, ...)
#     marginal_effects : layer1_marginal_scan()$scan_result reformatted as:
#                        data.frame(TRAIT, SNPID, A, SE, P_Value)
#                        where A / SE / P_Value are numeric
#     cond_effects     : fit_conditional_model()
#                        data.frame(trait, locus, theta_cond, se_cond,
#                                   pval_cond, sig_cond)
#     exposure/outcome : trait name strings
#     -> character vector of valid IV SNP IDs
#
#   estimate_causal_gls(theta_exp, theta_out, ld_matrix, Y, X_iv, ...)
#     theta_exp / theta_out : named numeric vectors (names = SNPID)
#     ld_matrix             : k x k matrix of |r| values; NULL -> identity
#     Y                     : n x m phenotype matrix (for bootstrap SE)
#     X_iv                  : n x k genotype matrix of IVs (for bootstrap SE)
#     -> list(gamma, se, z, pval, n_iv, method)
#
#   bidirectional_mr(marginal_effects, cond_effects, traits, Y, X_all, ...)
#     traits  : character(2) c("TraitA", "TraitB")
#     X_all   : full n x p genotype matrix; columns subsetted to IVs internally
#     -> list(AB = list(gamma, se, z, pval, sig, n_iv),
#             BA = list(...),
#             iv_AB, iv_BA)
#
#   test_mediation_consistency(theta_B_marg, theta_B_cond)
#     theta_B_marg : named numeric: IV marginal effects on outcome B
#     theta_B_cond : named numeric: IV conditional effects on B | A
#     -> list(pval, consistent, delta_median, delta)
# ==============================================================================


# ------------------------------------------------------------------------------
# Internal helper: F-statistic
# ------------------------------------------------------------------------------

#' @keywords internal
compute_f_stat <- function(A, SE) {
  A_num <- suppressWarnings(as.numeric(A))
  SE_num <- suppressWarnings(as.numeric(SE))
  (A_num / SE_num)^2
}


# ------------------------------------------------------------------------------
# Internal helper: greedy LD pruning
# ------------------------------------------------------------------------------

#' @keywords internal
ld_prune_ivs <- function(loci, ld_mat, r2_threshold = 0.1) {
  if (length(loci) <= 1L) {
    return(loci)
  }
  keep <- character(0L)
  for (l in loci) {
    if (length(keep) == 0L) {
      keep <- l
      next
    }
    max_r <- max(abs(ld_mat[l, keep, drop = FALSE]))
    if (max_r < r2_threshold) keep <- c(keep, l)
  }
  keep
}


# ------------------------------------------------------------------------------
# Step 1: IV selection (Section 2.3.1a, four criteria)
# ------------------------------------------------------------------------------

#' Select valid instrumental variables for MR
#'

#' IV 筛选：排他性筛选支持成对 / 全条件两种投影模式
#'
#' @param screen_mode "pairwise"（结局只对暴露投影 B|A，挡得住混杂工具）或
#'   "full"（结局对所有其他性状投影 B|A,C,...，混杂被投影掉、污染工具会漏进来）。
#'   m=2 时两者等价；m>=3 时才有区别。包里固定 "pairwise"。
#' @inheritParams 其余同原版（marginal_effects, Y_residual, X_loci,
#'   exposure, outcome, alpha1, alpha2, F_threshold, r2_threshold, ld_matrix, V）
#' @return 有效 IV 的 SNPID 字符向量
#' @export
select_instruments_pairwise <- function(marginal_effects,
                                        Y_residual,
                                        X_loci,
                                        exposure,
                                        outcome,
                                        alpha1,
                                        alpha2 = 0.05,
                                        F_threshold = 10,
                                        r2_threshold = 0.1,
                                        ld_matrix = NULL,
                                        V = NULL,
                                        screen_mode = c("pairwise", "full")) {
  screen_mode <- match.arg(screen_mode)
  Y_residual <- as.matrix(Y_residual)
  if (is.null(colnames(Y_residual))) {
    colnames(Y_residual) <- paste0("Trait", seq_len(ncol(Y_residual)))
  }
  if (is.null(V)) V <- stats::cov(Y_residual)
  m <- ncol(Y_residual)
  ie <- if (is.character(exposure)) match(exposure, colnames(Y_residual)) else as.integer(exposure)
  io <- if (is.character(outcome)) match(outcome, colnames(Y_residual)) else as.integer(outcome)
  if (is.na(ie) || is.na(io)) stop("exposure / outcome 未在 Y_residual 列名中找到。")

  # 准则1：相关性 —— 对暴露边际显著
  exp_sig <- marginal_effects[
    marginal_effects$TRAIT == exposure &
      suppressWarnings(as.numeric(marginal_effects$P_Value)) < alpha1,
    "SNPID"
  ]
  if (length(exp_sig) == 0L) {
    warning(sprintf("暴露 '%s' 无边际显著位点，无法构造 IV。", exposure))
    return(character(0L))
  }

  # 准则2：排他性 —— 按 screen_mode 决定结局对谁投影
  #   pairwise: cond = {exposure}            （挡混杂工具）
  #   full    : cond = 所有其他性状           （混杂被投影掉，污染工具漏进来）
  cond_idx <- if (screen_mode == "pairwise") ie else setdiff(seq_len(m), io)
  y_cond <- .project_outcome_on(Y_residual, io, cond_idx, V)

  cand_geno <- intersect(exp_sig, colnames(X_loci))
  excl_ok <- character(0L)
  for (snp in cand_geno) {
    x <- X_loci[, snp]
    if (var(x) < 1e-10) next
    fit <- lm(y_cond ~ x)
    cf <- summary(fit)$coefficients
    if (nrow(cf) < 2L) next
    if (cf[2L, 4L] > alpha2) excl_ok <- c(excl_ok, snp) # 条件效应不显著 → 保留
  }
  candidates <- intersect(exp_sig, excl_ok)
  if (length(candidates) == 0L) {
    warning(sprintf(
      "方向 %s -> %s（%s）：无位点通过排他性筛选。",
      exposure, outcome, screen_mode
    ))
    return(character(0L))
  }

  # 准则3：工具强度 F > F_threshold
  exp_df <- marginal_effects[
    marginal_effects$TRAIT == exposure & marginal_effects$SNPID %in% candidates,
  ]
  f_stats <- setNames(compute_f_stat(exp_df$A, exp_df$SE), exp_df$SNPID)
  f_pass <- names(f_stats)[!is.na(f_stats) & f_stats > F_threshold]
  candidates <- intersect(candidates, f_pass)
  if (length(candidates) == 0L) {
    warning(sprintf("方向 %s -> %s：无位点通过 F > %g。", exposure, outcome, F_threshold))
    return(character(0L))
  }

  # 准则4：LD 修剪
  if (!is.null(ld_matrix) && length(candidates) > 1L) {
    shared <- intersect(candidates, rownames(ld_matrix))
    if (length(shared) > 1L) {
      candidates <- ld_prune_ivs(shared, ld_matrix[shared, shared, drop = FALSE], r2_threshold)
    }
  }

  if (length(candidates) < 3L) {
    warning(sprintf(
      "方向 %s -> %s：筛选后仅 %d 个 IV（需≥3）。",
      exposure, outcome, length(candidates)
    ))
    return(character(0L))
  }
  candidates
}
# ------------------------------------------------------------------------------
# Step 2: GLS causal effect estimation (Eq. 7) + Bootstrap SE
# ------------------------------------------------------------------------------

#' GLS causal effect estimation with Bootstrap standard error
#'
#' Implements Eq. 7 of the paper:
#' \deqn{
#'   \hat{\gamma}_{A \to B} =
#'   \frac{\hat{\boldsymbol{\theta}}_A^\top \mathbf{R}^{-1} \hat{\boldsymbol{\theta}}_B}
#'        {\hat{\boldsymbol{\theta}}_A^\top \mathbf{R}^{-1} \hat{\boldsymbol{\theta}}_A}
#' }
#'
#' After LD pruning (|r| < 0.1), R approximates the identity matrix and GLS
#' reduces to weighted least squares (WLS) with weights w_i = 1 / Var(theta_i).
#' Standard error is estimated via individual-level bootstrap (n_boot resamples)
#' when Y and X_iv are provided, or by the delta method otherwise.
#'
#' @param theta_exp Named numeric vector. IV marginal effects on exposure.
#'   Names must be SNPID strings matching those in theta_out.
#' @param theta_out Named numeric vector. IV marginal effects on outcome.
#' @param ld_matrix Optional k x k matrix of |r| values for the IV set.
#'   NULL uses the identity matrix (appropriate after LD pruning).
#' @param Y n x m individual-level phenotype matrix. Required for bootstrap SE.
#' @param X_iv n x k genotype matrix for the IV set (columns = IVs used).
#' @param exp_col Integer. Column index of exposure trait in Y. Default 1.
#' @param out_col Integer. Column index of outcome trait in Y. Default 2.
#' @param n_boot Integer. Number of bootstrap resamples. Default 1000.
#'
#' @return Named list:
#' \describe{
#'   \item{gamma}{Point estimate of causal effect}
#'   \item{se}{Standard error (bootstrap or delta method)}
#'   \item{z}{Wald test statistic}
#'   \item{pval}{Two-sided p-value}
#'   \item{n_iv}{Number of IVs used}
#'   \item{method}{"bootstrap" or "delta"}
#' }
#'
#' @export
estimate_causal_gls <- function(theta_exp,
                                theta_out,
                                ld_matrix = NULL,
                                Y = NULL,
                                X_iv = NULL,
                                exp_col = 1L,
                                out_col = 2L,
                                n_boot = 1000L) {
  k <- length(theta_exp)
  stopifnot(length(theta_out) == k, k >= 3L)

  # GLS point estimate
  if (is.null(ld_matrix)) {
    R_inv <- diag(k)
  } else {
    loci <- names(theta_exp)
    sub <- ld_matrix[loci, loci, drop = FALSE]
    R_inv <- tryCatch(
      solve(sub),
      error = function(e) {
        warning("LD matrix is singular; falling back to identity matrix.")
        diag(k)
      }
    )
  }

  numerator <- as.numeric(t(theta_exp) %*% R_inv %*% theta_out)
  denominator <- as.numeric(t(theta_exp) %*% R_inv %*% theta_exp)
  gamma_hat <- numerator / denominator

  # Standard error
  method <- "delta"

  if (!is.null(Y) && !is.null(X_iv)) {
    method <- "bootstrap"
    n <- nrow(Y)
    gammas <- numeric(n_boot)

    for (b in seq_len(n_boot)) {
      idx <- sample.int(n, n, replace = TRUE)
      Yb <- Y[idx, , drop = FALSE]
      Xb <- X_iv[idx, , drop = FALSE]

      # Re-estimate IV effects via OLS in bootstrap sample
      tA_b <- apply(Xb, 2, function(x) {
        coef(lm.fit(cbind(1, x), Yb[, exp_col]))[2L]
      })
      tB_b <- apply(Xb, 2, function(x) {
        coef(lm.fit(cbind(1, x), Yb[, out_col]))[2L]
      })

      denom_b <- sum(tA_b^2)
      gammas[b] <- if (abs(denom_b) > 1e-12) sum(tA_b * tB_b) / denom_b else NA_real_
    }

    se_gamma <- sd(gammas, na.rm = TRUE)
  } else {
    # Delta method (Fieller 1954 approximation)
    sigma2 <- sum((theta_out - gamma_hat * theta_exp)^2) / max(k - 1L, 1L)
    se_gamma <- if (abs(denominator) > 1e-12) sqrt(sigma2 / denominator) else NA_real_
  }

  z <- gamma_hat / se_gamma
  pval <- 2 * pnorm(-abs(z))

  list(
    gamma  = gamma_hat,
    se     = se_gamma,
    z      = z,
    pval   = pval,
    n_iv   = k,
    method = method
  )
}

# ------------------------------------------------------------------------------
# Step 3: Bidirectional MR (pairwise-conditional IV screening version)
# ------------------------------------------------------------------------------

#' Bidirectional Mendelian Randomization (pairwise conditional screening)
#'
#' Runs MR in both directions (A->B and B->A) for a pair of traits. For each
#' direction, IVs are selected by select_instruments_pairwise(), whose exclusion
#' criterion uses a DIRECTION-SPECIFIC pairwise conditional projection
#' (outcome | exposure) computed on the fly. This blocks confounded instruments
#' that leak through a common upstream trait in m >= 3 settings — unlike a full
#' conditional projection, which would project the confounder out and let those
#' instruments pass.
#'
#' NOTE: the old `cond_effects` argument is removed. Pairwise conditional effects
#' are direction-specific and cannot be precomputed once and shared; they are now
#' computed inside each direction from `Y_residual`.
#'
#' @param marginal_effects Data frame with columns TRAIT, SNPID, A, SE, P_Value.
#'   Build it from layer1_marginal_scan()$scan_result by renaming
#'   (snp_id->SNPID, trait->TRAIT, beta->A, se->SE, p_value->P_Value) BEFORE
#'   calling this function.
#' @param traits Character vector length 2: c("TraitA", "TraitB").
#' @param Y_residual n x m residual phenotype matrix (column names = trait names),
#'   same convention as compute_conditional_phenotype(). Used to compute the
#'   pairwise conditional phenotype for IV exclusion screening. Required.
#' @param Y n x m individual-level phenotype matrix for bootstrap SE. NULL ->
#'   delta-method SE. May be the same matrix as Y_residual in simulations.
#' @param X_all n x p full genotype matrix (column names = SNPID). Passed to the
#'   IV selector (subset to candidates internally) and used for bootstrap SE.
#' @param ld_matrix Optional LD matrix (|r|). NULL skips LD pruning / uses identity.
#' @param alpha1 Genome-wide relevance threshold (use l1$threshold, not raw 0.05).
#' @param alpha2 Exclusion compatibility threshold. Default 0.05.
#' @param F_threshold Weak instrument filter. Default 10.
#' @param r2_threshold LD pruning threshold. Default 0.1.
#' @param alpha_mr MR significance threshold. Default 0.05.
#' @param n_boot Bootstrap resamples for SE. Default 1000.
#' @param V Optional precomputed m x m phenotypic covariance. NULL -> estimated
#'   inside compute_pairwise_conditional() from Y_residual each call.
#'
#' @return list(AB, BA, iv_AB, iv_BA) — same structure as before.
#' @export
bidirectional_mr <- function(marginal_effects,
                             traits,
                             Y_residual,
                             Y = NULL,
                             X_all = NULL,
                             ld_matrix = NULL,
                             alpha1,
                             alpha2 = 0.05,
                             F_threshold = 10,
                             r2_threshold = 0.1,
                             alpha_mr = 0.05,
                             n_boot = 1000L,
                             V = NULL,
                             screen_mode = c("pairwise", "full")) {
  screen_mode <- match.arg(screen_mode)
  stopifnot(length(traits) == 2L)
  trA <- traits[1L]
  trB <- traits[2L]

  # Precompute V once (optional) so both directions share the same covariance
  # estimate instead of re-estimating it inside each pairwise projection.
  if (is.null(V)) V <- stats::cov(as.matrix(Y_residual))

  # Helper: run MR for one direction
  run_direction <- function(exposure, outcome, exp_col, out_col) {
    # ---- IV selection: pairwise conditional screening (outcome | exposure) ----
    ivs <- select_instruments_pairwise(
      marginal_effects = marginal_effects,
      Y_residual       = Y_residual,
      X_loci           = X_all,
      exposure         = exposure,
      outcome          = outcome,
      alpha1           = alpha1,
      alpha2           = alpha2,
      F_threshold      = F_threshold,
      r2_threshold     = r2_threshold,
      ld_matrix        = ld_matrix,
      V                = V,
      screen_mode      = screen_mode
    )

    if (length(ivs) == 0L) {
      return(list(
        gamma = NA_real_, se = NA_real_, z = NA_real_,
        pval = NA_real_, sig = FALSE, n_iv = 0L,
        method = NA_character_, ivs = character(0L)
      ))
    }

    # ---- Extract IV effects on exposure and outcome (unchanged) ----------------
    get_eff <- function(trait) {
      rows <- marginal_effects[
        marginal_effects$TRAIT == trait &
          marginal_effects$SNPID %in% ivs,
      ]
      setNames(suppressWarnings(as.numeric(rows$A)), rows$SNPID)[ivs]
    }

    theta_exp <- get_eff(exposure)
    theta_out <- get_eff(outcome)

    ok <- !is.na(theta_exp) & !is.na(theta_out)
    theta_exp <- theta_exp[ok]
    theta_out <- theta_out[ok]
    ivs_used <- names(theta_exp)

    if (length(ivs_used) < 3L) {
      warning(sprintf(
        "Direction %s -> %s: only %d IV(s) have effects on both traits. Skipping MR.",
        exposure, outcome, length(ivs_used)
      ))
      return(list(
        gamma = NA_real_, se = NA_real_, z = NA_real_,
        pval = NA_real_, sig = FALSE, n_iv = length(ivs_used),
        method = NA_character_, ivs = ivs_used
      ))
    }

    # ---- Individual-level data for bootstrap SE (unchanged) --------------------
    X_iv_boot <- NULL
    Y_pair <- NULL
    if (!is.null(X_all) && !is.null(Y)) {
      cols <- intersect(ivs_used, colnames(X_all))
      if (length(cols) == length(ivs_used)) {
        X_iv_boot <- X_all[, cols, drop = FALSE]
        Y_pair <- Y[, c(exp_col, out_col), drop = FALSE]
      }
    }

    gls <- estimate_causal_gls(
      theta_exp = theta_exp,
      theta_out = theta_out,
      ld_matrix = ld_matrix,
      Y         = Y_pair,
      X_iv      = X_iv_boot,
      exp_col   = 1L,
      out_col   = 2L,
      n_boot    = n_boot
    )

    c(gls, list(
      sig = !is.na(gls$pval) && gls$pval < alpha_mr,
      ivs = ivs_used
    ))
  }

  # Resolve column indices for each trait in Y / Y_residual
  trait_col <- function(tr) {
    ref <- if (!is.null(Y)) Y else Y_residual
    if (!is.null(colnames(ref))) {
      idx <- which(colnames(ref) == tr)
      if (length(idx) == 1L) {
        return(idx)
      }
    }
    if (tr == trA) 1L else 2L
  }

  mr_AB <- run_direction(trA, trB, trait_col(trA), trait_col(trB))
  mr_BA <- run_direction(trB, trA, trait_col(trB), trait_col(trA))

  list(
    AB    = mr_AB,
    BA    = mr_BA,
    iv_AB = mr_AB$ivs,
    iv_BA = mr_BA$ivs
  )
}

# ------------------------------------------------------------------------------
# Step 4: Mediation consistency test (Section 2.3.2)
# ------------------------------------------------------------------------------

#' Reverse conditional consistency test for MR-supported causal pathway
#'
#' For a causal pathway A -> B supported by MR, verifies that controlling for
#' the upstream trait A systematically attenuates IV effects on the downstream
#' trait B. Tests whether delta_l = theta_{Bl}^{marg} - theta_{B|A,l}^{cond}
#' is systematically greater than zero using a one-sided Wilcoxon signed-rank test.
#'
#' Interpretation:
#' \itemize{
#'   \item delta > 0 for all IVs -> complete mediation (Class 3)
#'   \item delta > 0 but theta_{B|A}^{cond} still significant -> partial mediation (Class 4)
#'   \item delta not systematically positive -> inconsistent with A -> B
#' }
#'
#' @param theta_B_marg Named numeric vector. IV marginal effects on outcome B.
#'   Names are SNPID strings.
#' @param theta_B_cond Named numeric vector. IV conditional effects on B given A.
#'   Names are SNPID strings. Alignment is done by name intersection internally.
#'
#' @return Named list:
#' \describe{
#'   \item{pval}{One-sided Wilcoxon signed-rank p-value}
#'   \item{consistent}{TRUE if pval < 0.05}
#'   \item{delta_median}{Median of delta_l = theta_marg - theta_cond (positive = attenuation)}
#'   \item{delta}{Named numeric vector of per-IV delta values}
#' }
#'
#' @export
test_mediation_consistency <- function(theta_B_marg, theta_B_cond) {
  shared <- intersect(names(theta_B_marg), names(theta_B_cond))

  if (length(shared) < 3L) {
    warning(sprintf(
      "Only %d shared IV(s) between marginal and conditional effects. ",
      "Consistency test unreliable (minimum 3 required).",
      length(shared)
    ))
    return(list(
      pval         = NA_real_,
      consistent   = FALSE,
      delta_median = NA_real_,
      delta        = numeric(0L)
    ))
  }

  # delta_l = marginal effect - conditional effect
  # Positive delta means controlling for A attenuates the IV effect on B,
  # consistent with mediation through A.
  delta <- theta_B_marg[shared] - theta_B_cond[shared]

  wt <- tryCatch(
    wilcox.test(delta, mu = 0, alternative = "greater"),
    error = function(e) list(p.value = NA_real_)
  )

  list(
    pval         = wt$p.value,
    consistent   = !is.na(wt$p.value) && wt$p.value < 0.05,
    delta_median = median(delta),
    delta        = delta
  )
}


# ------------------------------------------------------------------------------
# (1) MVMR 估计器 —— 加进 mr.R
# ------------------------------------------------------------------------------

#' Multivariable Mendelian Randomization (多暴露直接因果效应)
#'
#' 对结局 outcome，把多个 exposures 同时作为暴露，用多元 GLS 估计每个暴露对
#' 结局的【直接】因果效应（控制其他暴露）。工具变量取"对至少一个暴露边际显著"
#' 的并集；LD 修剪后 R≈I，GLS 退化为多元 IVW。
#'
#' 与 bidirectional_mr 的区别：后者两两估计、靠成对条件投影挡混杂工具；MVMR
#' 把可观测混杂当协暴露显式纳入，从而吸收其后门效应。代价：必须预先知道协暴露
#' 身份（指错——如把下游后代当协暴露——会 over-adjustment 致偏）。
#'
#' @param marginal_effects data.frame(TRAIT, SNPID, A, SE, P_Value)
#' @param exposures 字符向量，长度>=2。第 1 个为焦点暴露，其余为协暴露。
#' @param outcome 结局性状名
#' @param alpha1 相关性阈值（用 l1$threshold）
#' @param F_threshold 弱工具阈值，默认 10
#' @param r2_threshold LD 修剪阈值，默认 0.1
#' @param ld_matrix 可选 LD 矩阵；NULL 则不修剪、R=I
#' @param focal 报告哪个暴露的直接效应；默认 exposures[1]
#' @return list(gamma 焦点暴露直接效应, se, z, pval, n_iv,
#'   gamma_all 所有暴露的直接效应向量, exposures, focal)
#' @export
mvmr_estimate <- function(marginal_effects,
                          exposures,
                          outcome,
                          alpha1,
                          F_threshold = 10,
                          r2_threshold = 0.1,
                          ld_matrix = NULL,
                          focal = NULL) {
  stopifnot(length(exposures) >= 2L)
  if (is.null(focal)) focal <- exposures[1L]

  pull <- function(trait) {
    rows <- marginal_effects[marginal_effects$TRAIT == trait, ]
    list(
      beta = setNames(suppressWarnings(as.numeric(rows$A)), rows$SNPID),
      se = setNames(suppressWarnings(as.numeric(rows$SE)), rows$SNPID),
      p = setNames(suppressWarnings(as.numeric(rows$P_Value)), rows$SNPID)
    )
  }

  # 工具：对任一暴露边际显著的并集
  inst <- character(0L)
  for (e in exposures) {
    pe <- pull(e)
    inst <- union(inst, names(pe$p)[!is.na(pe$p) & pe$p < alpha1])
  }
  if (length(inst) < length(exposures) + 1L) {
    warning(sprintf(
      "MVMR(%s -> %s)：工具数 %d 不足。",
      paste(exposures, collapse = ","), outcome, length(inst)
    ))
    return(list(
      gamma = NA_real_, se = NA_real_, z = NA_real_, pval = NA_real_,
      n_iv = length(inst), gamma_all = NULL, exposures = exposures, focal = focal
    ))
  }

  # 弱工具过滤：保留对【至少一个暴露】F>阈值的工具
  strong <- rep(FALSE, length(inst))
  names(strong) <- inst
  for (e in exposures) {
    pe <- pull(e)
    f <- (pe$beta[inst] / pe$se[inst])^2
    strong <- strong | (!is.na(f) & f > F_threshold)
  }
  inst <- inst[strong]

  # LD 修剪
  if (!is.null(ld_matrix) && length(inst) > 1L) {
    shared <- intersect(inst, rownames(ld_matrix))
    if (length(shared) > 1L) inst <- ld_prune_ivs(shared, ld_matrix[shared, shared, drop = FALSE], r2_threshold)
  }
  if (length(inst) < length(exposures) + 1L) {
    warning(sprintf(
      "MVMR(%s -> %s)：过滤后工具数 %d 不足。",
      paste(exposures, collapse = ","), outcome, length(inst)
    ))
    return(list(
      gamma = NA_real_, se = NA_real_, z = NA_real_, pval = NA_real_,
      n_iv = length(inst), gamma_all = NULL, exposures = exposures, focal = focal
    ))
  }

  # 设计矩阵 Theta (k x n_exp) 与结局向量 theta_out (k)
  Theta <- vapply(exposures, function(e) pull(e)$beta[inst], numeric(length(inst)))
  if (is.null(dim(Theta))) Theta <- matrix(Theta, ncol = length(exposures))
  colnames(Theta) <- exposures
  theta_out <- pull(outcome)$beta[inst]

  ok <- stats::complete.cases(Theta) & !is.na(theta_out)
  Theta <- Theta[ok, , drop = FALSE]
  theta_out <- theta_out[ok]
  k <- nrow(Theta)
  if (k < length(exposures) + 1L) {
    return(list(
      gamma = NA_real_, se = NA_real_, z = NA_real_, pval = NA_real_,
      n_iv = k, gamma_all = NULL, exposures = exposures, focal = focal
    ))
  }

  # 多元 GLS（R=I 退化为多元 IVW）
  if (is.null(ld_matrix)) {
    XtX <- crossprod(Theta) # Theta' Theta
    Xty <- crossprod(Theta, theta_out)
  } else {
    Ri <- tryCatch(solve(ld_matrix[rownames(Theta), rownames(Theta), drop = FALSE]),
      error = function(e) diag(k)
    )
    XtX <- t(Theta) %*% Ri %*% Theta
    Xty <- t(Theta) %*% Ri %*% theta_out
  }
  XtX_inv <- solve(XtX)
  gamma <- as.numeric(XtX_inv %*% Xty)
  names(gamma) <- exposures

  # 残差法标准误
  resid <- as.numeric(theta_out - Theta %*% gamma)
  sigma2 <- sum(resid^2) / max(k - length(exposures), 1L)
  cov_g <- sigma2 * XtX_inv
  se_all <- sqrt(diag(cov_g))
  names(se_all) <- exposures

  g_focal <- gamma[focal]
  se_focal <- se_all[focal]
  z <- g_focal / se_focal

  list(
    gamma = as.numeric(g_focal),
    se = as.numeric(se_focal),
    z = as.numeric(z),
    pval = 2 * stats::pnorm(-abs(z)),
    n_iv = k,
    gamma_all = gamma,
    exposures = exposures,
    focal = focal
  )
}


# ------------------------------------------------------------------------------
# (2) 驱动：confound_obs × MVMR（先只跑这一格）
# ------------------------------------------------------------------------------

#' confound_obs 上批量跑 MVMR(A,C -> B)，报告焦点暴露 A 的直接效应
#'
#' 真值 tau_AB = 0。预期 MVMR 把 gamma_A 修回 ~0、typeI ~名义水平。
#' 与成对两两MR（残余偏倚 ~0.10、typeI ~0.97）对比，凸显 MVMR 在可观测混杂下的价值。
#'
#' @return 1 行 data.frame
#' @export
run_mvmr_confound_obs <- function(n = 2000L, n_rep = 200L,
                                  n_conf = 15L, conf_pve = 0.04,
                                  tau_CA = 0.5, tau_CB = 0.5,
                                  alpha1_level = 0.05, verbose = TRUE) {
  gA <- numeric(n_rep)
  sig <- logical(n_rep)
  niv <- numeric(n_rep)

  for (s in seq_len(n_rep)) {
    if (verbose && s %% 20 == 0) message(sprintf("  rep %d / %d", s, n_rep))

    dat <- generate_dataset_IV(
      n = n, topology = "confound_obs", seed = s,
      n_conf = n_conf, conf_pve = conf_pve,
      tau_CA = tau_CA, tau_CB = tau_CB
    )

    l1 <- layer1_marginal_scan(dat$X, dat$Y, alpha = alpha1_level, verbose = FALSE)
    me <- with(l1$scan_result, data.frame(
      TRAIT = trait, SNPID = snp_id, A = beta, SE = se, P_Value = p_value,
      stringsAsFactors = FALSE
    ))

    mv <- tryCatch(
      mvmr_estimate(
        marginal_effects = me,
        exposures        = c("TraitA", "TraitC"), # 焦点 A + 协暴露 C（混杂）
        outcome          = "TraitB",
        alpha1           = l1$threshold,
        focal            = "TraitA"
      ),
      error = function(e) {
        warning(conditionMessage(e))
        NULL
      }
    )

    if (is.null(mv) || is.na(mv$gamma)) {
      gA[s] <- NA
      next
    }
    gA[s] <- mv$gamma
    sig[s] <- !is.na(mv$pval) && mv$pval < 0.05
    niv[s] <- mv$n_iv
  }

  ok <- !is.na(gA)
  data.frame(
    estimator = "MVMR(A,C->B)",
    topology = "confound_obs",
    tau_AB_true = 0,
    gamma_A_mean = mean(gA[ok]),
    gamma_A_sd = sd(gA[ok]),
    typeI_rate = mean(sig[ok]),
    mean_n_iv = mean(niv[ok]),
    n_valid = sum(ok),
    stringsAsFactors = FALSE
  )
}

# ---- 跑法 ----
# source 四个依赖后（mvmr_estimate 已并入 mr.R）：
#   res_mvmr <- run_mvmr_confound_obs(n = 2000, n_rep = 200, n_conf = 15, conf_pve = 0.04)
#   print(res_mvmr)
# 预期：gamma_A_mean ≈ 0（对比成对两两MR ~0.10）、typeI_rate ≈ 名义 5%（对比 ~0.97）。


# 放在 mr.R 里即可（内部函数，不需要export）
.compute_causal_adjusted <- function(Y, tau_matrix, traits) {
  # 输入：
  #   Y          n×m 原始表型矩阵
  #   tau_matrix m×m 当前tau估计值（Tau[j,i]=性状i对性状j的效应）
  #   traits     m维性状名向量
  #
  # 输出：
  #   Y_adj n×m 条件表型矩阵，第i列 = Y_i - Σ_j tau[i,j] × Y_j

  m <- length(traits)
  Y_adj <- Y

  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      if (i == j) next
      tau_ij <- tau_matrix[i, j]
      if (abs(tau_ij) > 1e-10) {
        # 从性状i里去掉"性状j通过因果链传递来的效应"
        Y_adj[, i] <- Y_adj[, i] - tau_ij * Y[, j]
      }
    }
  }
  Y_adj
}
