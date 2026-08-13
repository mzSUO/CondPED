#' Simulate CondPED data with population structural truth (revised v1.0)
#'
#' Generates individual-level multi-trait data under the unified
#' multi-signal phenotype model
#' \deqn{Y = WB + X_Q B_Q + U + E}
#' with a focal locus of one or more possibly linked causal signals,
#' plus the complete population structural truth (candidate traits,
#' full representation-loss map, minimum representative sets and
#' irreducible modules per signal).
#'
#' @details
#' All focal-region markers are excluded from the background GRM, and
#' \deqn{\Sigma_Q = B_Q^\top \Sigma_X B_Q,\qquad
#' \Sigma_G^{bg} = \Sigma_G^{tot} - \Sigma_Q}
#' Only generations with positive semi-definite \eqn{\Sigma_G^{bg}} are
#' accepted; if no generation qualifies within `max_attempts`, the
#' function returns status `invalid_covariance` (no silent clipping).
#'
#' Per-signal PVE scaling reuses the single-signal scaler (primary
#' signal: `locus_pve`; secondary signals: `secondary_signal_pve`);
#' these are marginal calibration targets only and are never summed
#' into a regional PVE. After row-wise scaling, the full
#' \eqn{\Sigma_Q} is recomputed and the truth records both the
#' realized per-signal PVE and the regional covariance contribution.
#'
#' Representability scenarios (R1/R2/R3) use the closed-form
#' target-rho construction (Methods): with R = A \ S*,
#' \eqn{\beta_R = \Gamma_{R|S^*}^\top \beta_{S^*} + \delta d},
#' \eqn{d^\top \Omega^{-1} d = 1},
#' \eqn{\delta = \sqrt{\rho^\star A_0 / (1 - \rho^\star)}}.
#' The full rho map is then exactly enumerated with the same Stage-1
#' loss path, and Rep/Irr are extracted from that map only.
#'
#' @param n Number of individuals.
#' @param m Number of traits.
#' @param p Number of markers.
#' @param maf_range Range for base minor allele frequencies.
#' @param h2 Trait heritabilities (length 1 or `m`).
#' @param R_G,R_E Optional `m x m` correlation matrices.
#' @param experiment One of `"signal_resolution"`,
#'   `"trait_representation"`, `"end_to_end"`.
#' @param scenario Scenario within the experiment (see the frozen
#'   registry; e.g. `"null"`, `"single_multi_trait"`,
#'   `"two_linked_trait_specific"`, `"highly_representable"`).
#' @param n_signals Optional consistency check on the scenario's
#'   signal count.
#' @param locus_pve Marginal PVE calibration target of the primary
#'   signal.
#' @param secondary_signal_pve Marginal PVE calibration target of each
#'   secondary signal.
#' @param local_ld Target LD level inside the focal region:
#'   `"none"` (0), `"low"` (0.1), `"moderate"` (0.3) or `"high"`
#'   (0.6).
#' @param target_r2 Optional explicit target pairwise r^2 inside the
#'   focal region (overrides `local_ld`).
#' @param local_region_size Number of markers in the focal region.
#' @param tolerance Primary representation-loss tolerance tau.
#' @param target_loss Optional closed-form target rho* for the
#'   representability scenarios (default 0.02).
#' @param target_representing_set Optional target representing set
#'   S* for the representability scenarios.
#' @param correlation Background trait correlation structure.
#' @param structured,n_groups,fst Population-structure settings for
#'   the background GRM.
#' @param truth_margin Minimum distance of structural-scenario losses
#'   from `tolerance`.
#' @param max_attempts Maximum generation attempts.
#' @param seed Optional random seed.
#' @param return_latent Logical; keep `U` and `E`.
#'
#' @return A list with components `Y`, `W`, `G`, `K_bg`,
#'   `focal_region` (marker ids of the focal region), `causal_index`,
#'   `truth`, `generator`, `status` and `diagnostics`. The truth
#'   fields follow the unified v1.0 schema (`loci`, `signals`,
#'   `causal_markers`, `signal_count`, `local_ld`,
#'   `realized_signal_pve`, `beta`, `candidate_traits`,
#'   `effect_breadth`, `effect_direction`, `subset_table`,
#'   `representation_map`, `minimum_representative_sets`,
#'   `irreducible_modules`, `tolerance_path`, `Sigma_G_total`,
#'   `Sigma_G_bg`, `Sigma_Q`, `Sigma_E`, `Sigma_P_total`,
#'   `Sigma_P_ref`).
#' @export
simulate_condped_data <- function(
  n = 1000L,
  m = 4L,
  p = 1000L,
  maf_range = c(0.05, 0.50),
  h2 = rep(0.50, m),
  R_G = NULL,
  R_E = NULL,
  experiment = c(
    "signal_resolution",
    "trait_representation",
    "end_to_end"
  ),
  scenario,
  n_signals = NULL,
  locus_pve = 0.01,
  secondary_signal_pve = 0.005,
  local_ld = c("none", "low", "moderate", "high"),
  target_r2 = NULL,
  local_region_size = 50L,
  tolerance = 0.10,
  target_loss = NULL,
  target_representing_set = NULL,
  correlation = c("block", "independent"),
  structured = TRUE,
  n_groups = 8L,
  fst = 0.05,
  truth_margin = 0.02,
  max_attempts = 200L,
  seed = NULL,
  return_latent = FALSE
) {
  experiment <- match.arg(experiment)
  if (missing(scenario) || !is.character(scenario) ||
      length(scenario) != 1L) {
    .stop_invalid_input("scenario must be a single scenario name.")
  }
  spec <- .scenario_spec(experiment, scenario)
  correlation <- match.arg(correlation)
  .check_count(n, "n", min = 2L)
  .check_count(m, "m", min = 1L)
  .check_count(p, "p", min = 3L)
  .check_count(n_groups, "n_groups", min = 2L)
  .check_count(max_attempts, "max_attempts", min = 1L)
  .check_count(local_region_size, "local_region_size", min = 2L)
  if (local_region_size > p) {
    .stop_invalid_input("local_region_size must not exceed p.")
  }
  if (!is.numeric(maf_range) || length(maf_range) != 2L ||
      any(!is.finite(maf_range)) || maf_range[1L] <= 0 ||
      maf_range[1L] >= maf_range[2L] || maf_range[2L] > 0.5) {
    .stop_invalid_input(
      "maf_range must satisfy 0 < maf_range[1] < maf_range[2] <= 0.5."
    )
  }
  if (!is.numeric(h2) || !(length(h2) %in% c(1L, m)) ||
      any(!is.finite(h2)) || any(h2 <= 0) || any(h2 >= 1)) {
    .stop_invalid_input(
      "h2 must have length 1 or m, with all entries in (0, 1)."
    )
  }
  h2 <- rep(h2, length.out = m)
  .check_prob01_vec <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 0 || x >= 1) {
      .stop_invalid_input("%s must be a single number in [0, 1).", name)
    }
  }
  .check_prob01_vec(locus_pve, "locus_pve")
  .check_prob01_vec(secondary_signal_pve, "secondary_signal_pve")
  .check_prob01_vec(tolerance, "tolerance")
  if (!is.null(target_loss)) {
    .check_prob01_vec(target_loss, "target_loss")
  }
  if (!is.null(target_r2)) {
    if (!is.numeric(target_r2) || length(target_r2) != 1L ||
        !is.finite(target_r2) || target_r2 < 0 || target_r2 > 1) {
      .stop_invalid_input("target_r2 must be a single number in [0, 1].")
    }
  }
  if (!is.numeric(truth_margin) || length(truth_margin) != 1L ||
      !is.finite(truth_margin) || truth_margin < 0) {
    .stop_invalid_input("truth_margin must be a single non-negative number.")
  }
  if (structured && (!is.numeric(fst) || length(fst) != 1L ||
                     !is.finite(fst) || fst <= 0 || fst >= 1)) {
    .stop_invalid_input("fst must be a single number in (0, 1).")
  }
  if (!is.logical(return_latent) || length(return_latent) != 1L ||
      is.na(return_latent)) {
    .stop_invalid_input("return_latent must be TRUE or FALSE.")
  }

  # ---- scenario spec / signal count -------------------------------------------
  q <- spec$q
  if (!is.null(n_signals)) {
    .check_count(n_signals, "n_signals", min = 0L)
    if (n_signals != q) {
      .stop_invalid_input(
        "scenario \"%s\" has %d signal(s) by definition; n_signals = %d conflicts.",
        scenario, q, n_signals
      )
    }
  }
  if (m < spec$min_m) {
    .stop_invalid_input(
      "scenario \"%s\" requires m >= %d.", scenario, spec$min_m
    )
  }

  # target LD inside the focal region: an explicitly supplied local_ld
  # or target_r2 always wins; otherwise the scenario default applies
  ld_arg <- if (missing(local_ld)) NULL else match.arg(local_ld)
  r2_target <- if (!is.null(target_r2)) {
    target_r2
  } else {
    ld_choice <- if (!is.null(ld_arg)) {
      ld_arg
    } else if (!is.null(spec$local_ld)) {
      spec$local_ld
    } else {
      "none"
    }
    unname(c(none = 0, low = 0.1, moderate = 0.3, high = 0.6)[ld_choice])
  }

  trait_names <- paste0("Trait", seq_len(m))
  ind_names <- paste0("Ind", seq_len(n))
  marker_names <- paste0("M", seq_len(p))

  # ---- covariance structures ---------------------------------------------------
  if (is.null(R_G)) {
    R_G <- if (correlation == "block") .block_cor(m) else diag(m)
  }
  if (is.null(R_E)) {
    R_E <- if (correlation == "block") .block_cor(m) else diag(m)
  }
  .check_correlation_matrix(R_G, m, "R_G")
  .check_correlation_matrix(R_E, m, "R_E")
  Sigma_G_total <- .cov_from_h2(h2, R_G)
  Sigma_E <- .cov_from_h2(1 - h2, R_E)
  Sigma_P_total <- Sigma_G_total + Sigma_E
  dimnames(Sigma_G_total) <- dimnames(Sigma_E) <-
    dimnames(Sigma_P_total) <- list(trait_names, trait_names)
  Sigma_P_ref <- Sigma_P_total   # frozen reference for all truth

  # resolve target representing set
  S_star <- NULL
  if (!is.null(target_representing_set)) {
    S_star <- .normalize_trait_set(target_representing_set, trait_names,
                                   "target_representing_set")
  } else if (!is.null(spec$default_S_star)) {
    S_star <- spec$default_S_star(m, trait_names)
  }
  rho_star <- if (!is.null(target_loss)) target_loss else 0.02

  if (!is.null(seed)) set.seed(seed)

  # ---- genotypes: background + focal LD region ----------------------------------
  maf <- stats::runif(p, min = maf_range[1L], max = maf_range[2L])
  region_start <- sample.int(p - local_region_size + 1L, 1L)
  region_cols <- region_start:(region_start + local_region_size - 1L)
  bg_cols <- setdiff(seq_len(p), region_cols)

  G <- matrix(0, nrow = n, ncol = p,
              dimnames = list(ind_names, marker_names))
  if (length(bg_cols) > 0L) {
    G[, bg_cols] <- .sim_genotypes(n, maf[bg_cols],
                                   structured = structured,
                                   n_groups = n_groups, fst = fst)
  }
  region <- .sim_ld_region(n, local_region_size, maf[region_cols],
                           r2_target)
  G[, region_cols] <- region$G

  # causal markers: evenly spaced inside the focal region
  causal_cols <- region_cols[round(seq(1, local_region_size,
                                       length.out = max(q, 1L)))]
  causal_index <- if (q > 0L) causal_cols else integer()
  position <- seq_len(p) * 1000
  chromosome <- rep("chr1", p)
  locus_id <- sprintf("chr1:%d-%d", min(position[region_cols]),
                      max(position[region_cols]))

  # ---- acceptance loop over population truth --------------------------------------
  psd_tol <- 1e-10
  tolerances <- unique(c(tolerance, 0.05, 0.10, 0.20))
  attempt <- 0L
  psd_failures <- 0L
  accepted <- FALSE
  B_Q <- matrix(0, nrow = max(q, 1L), ncol = m,
                dimnames = list(NULL, trait_names))
  truth_tabs <- NULL
  v_causal <- 2 * maf[causal_index] * (1 - maf[causal_index])
  Sigma_X <- matrix(0, max(q, 1L), max(q, 1L))
  if (q > 0L) {
    Sigma_X[] <- region$latent_r2
    diag(Sigma_X) <- 1
    Sigma_X <- Sigma_X * outer(sqrt(v_causal), sqrt(v_causal))
  }

  repeat {
    attempt <- attempt + 1L
    prop <- .propose_B_Q(spec, m, trait_names, Sigma_P_ref,
                         S_star, rho_star)
    B_Q <- prop$B_Q
    delta <- prop$delta
    # per-signal marginal PVE scaling (row-wise)
    if (q > 0L) {
      for (i in seq_len(q)) {
        tgt <- if (i == 1L) locus_pve else secondary_signal_pve
        nz <- which(B_Q[i, ] != 0)
        if (length(nz) > 0L && tgt > 0) {
          c_i <- .scale_signal_pve(B_Q[i, ], v_causal[i], tgt)
          B_Q[i, ] <- c_i * B_Q[i, ]
        }
      }
      # cross-architecture scale matching for R scenes
      if (spec$scale_match) {
        q_cur <- v_causal[1L] * as.numeric(
          crossprod(B_Q[1L, ], .safe_inverse(Sigma_P_ref)$inverse %*%
                      B_Q[1L, ]))
        q_target <- m * locus_pve
        if (q_cur > 0) B_Q <- B_Q * sqrt(q_target / q_cur)
      }
    }
    Sigma_Q <- if (q > 0L) {
      crossprod(B_Q, Sigma_X %*% B_Q)
    } else {
      matrix(0, m, m)
    }
    Sigma_G_bg <- Sigma_G_total - Sigma_Q
    min_eigen_bg <- .min_eigen_sym(Sigma_G_bg)
    if (min_eigen_bg < -psd_tol) {
      psd_failures <- psd_failures + 1L
      if (attempt >= max_attempts) break
      next
    }
    truth_tabs <- .truth_tables_per_signal(
      B_Q, Sigma_P_ref, tolerance, tolerances,
      locus_id, chromosome[causal_index], position[causal_index],
      marker_names[causal_index], trait_names
    )
    accepted <- .accept_scenario(spec, truth_tabs, tolerance,
                                 truth_margin)
    if (accepted || attempt >= max_attempts) break
  }

  # ---- data generation -------------------------------------------------------------
  Z_bg <- sweep(G[, bg_cols, drop = FALSE], 2L, 2 * maf[bg_cols], `-`)
  Z_bg <- sweep(Z_bg, 2L, sqrt(2 * maf[bg_cols] * (1 - maf[bg_cols])), `/`)
  K_bg <- .make_grm(Z_bg)
  dimnames(K_bg) <- list(ind_names, ind_names)
  eig_K <- eigen(K_bg, symmetric = TRUE)
  lambda_K <- eig_K$values
  K_rank <- sum(lambda_K > sqrt(.Machine$double.eps) * max(lambda_K))

  S_G <- .symmetric_psd_sqrt(Sigma_G_bg)
  S_E <- .symmetric_psd_sqrt(Sigma_E)
  Z_u <- matrix(stats::rnorm(n * m), nrow = n, ncol = m)
  U <- sweep(eig_K$vectors, 2L, sqrt(pmax(lambda_K, 0)), `*`) %*%
    Z_u %*% S_G
  E <- matrix(stats::rnorm(n * m), nrow = n, ncol = m) %*% S_E
  dimnames(U) <- dimnames(E) <- list(ind_names, trait_names)

  W <- matrix(1, nrow = n, ncol = 1,
              dimnames = list(ind_names, "Intercept"))
  Y <- matrix(0, n, m, dimnames = list(ind_names, trait_names))
  if (q > 0L) {
    Y <- G[, causal_index, drop = FALSE] %*% B_Q
  }
  Y <- Y + U + E

  # ---- truth assembly -----------------------------------------------------------------
  sigma2_P <- diag(Sigma_P_total)
  realized_pve <- list()
  beta_rows <- list()
  if (q > 0L) {
    for (i in seq_len(q)) {
      sid <- paste0(locus_id, "::S", i)
      realized_pve[[i]] <- data.frame(
        locus_id = locus_id, signal_id = sid, signal_order = i,
        target_pve = if (i == 1L) locus_pve else secondary_signal_pve,
        realized_marginal_signal_pve = mean(
          v_causal[i] * B_Q[i, ]^2 / sigma2_P),
        stringsAsFactors = FALSE
      )
      beta_rows[[i]] <- data.frame(
        locus_id = locus_id, signal_id = sid,
        representative_snp = marker_names[causal_index[i]],
        causal_marker = marker_names[causal_index[i]],
        trait = trait_names, beta = B_Q[i, ],
        stringsAsFactors = FALSE
      )
    }
  }
  pve_df <- if (length(realized_pve) > 0L) {
    do.call(rbind, realized_pve)
  } else {
    data.frame(
      locus_id = character(), signal_id = character(),
      signal_order = integer(), target_pve = numeric(),
      realized_marginal_signal_pve = numeric(), stringsAsFactors = FALSE
    )
  }
  beta_df <- if (length(beta_rows) > 0L) {
    do.call(rbind, beta_rows)
  } else {
    data.frame(
      locus_id = character(), signal_id = character(),
      representative_snp = character(), causal_marker = character(),
      trait = character(), beta = numeric(), stringsAsFactors = FALSE
    )
  }

  cand_traits <- stats::setNames(vector("list", q), character(q))
  breadth <- integer(q)
  direction <- character(q)
  if (q > 0L) {
    for (i in seq_len(q)) {
      sid <- paste0(locus_id, "::S", i)
      ct <- trait_names[B_Q[i, ] != 0]
      cand_traits[[i]] <- ct
      names(cand_traits)[i] <- sid
      breadth[i] <- length(ct)
      direction[i] <- .direction_of(B_Q[i, ][B_Q[i, ] != 0])
    }
  }

  truth <- list(
    loci = data.frame(
      locus_id = locus_id, chromosome = "chr1",
      start = min(position[region_cols]),
      end = max(position[region_cols]),
      n_region_markers = local_region_size, n_signals = q,
      stringsAsFactors = FALSE
    ),
    signals = if (q > 0L) {
      data.frame(
        locus_id = locus_id,
        signal_id = paste0(locus_id, "::S", seq_len(q)),
        signal_order = seq_len(q),
        representative_snp = marker_names[causal_index],
        chromosome = "chr1", position = position[causal_index],
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        locus_id = character(), signal_id = character(),
        signal_order = integer(), representative_snp = character(),
        chromosome = character(), position = numeric(),
        stringsAsFactors = FALSE
      )
    },
    causal_markers = marker_names[causal_index],
    signal_count = stats::setNames(q, locus_id),
    local_ld = data.frame(
      locus_id = locus_id, target_r2 = r2_target,
      realized_mean_r2 = region$mean_r2, realized_max_r2 = region$max_r2,
      stringsAsFactors = FALSE
    ),
    realized_signal_pve = pve_df,
    beta = beta_df,
    candidate_traits = cand_traits,
    effect_breadth = stats::setNames(breadth, names(cand_traits)),
    effect_direction = stats::setNames(direction, names(cand_traits)),
    subset_table = truth_tabs$subset_table,
    representation_map = truth_tabs$representation_map,
    minimum_representative_sets = truth_tabs$minimum_representative_sets,
    irreducible_modules = truth_tabs$irreducible_modules,
    tolerance_path = truth_tabs$tolerance_path,
    B_Q = B_Q,
    Sigma_X = Sigma_X,
    Sigma_Q = Sigma_Q,
    Sigma_G_total = Sigma_G_total,
    Sigma_G_bg = Sigma_G_bg,
    Sigma_E = Sigma_E,
    Sigma_P_total = Sigma_P_total,
    Sigma_P_ref = Sigma_P_ref
  )

  status <- .new_status(
    ok = accepted,
    code = if (accepted) "ok" else "invalid_covariance",
    message = if (accepted) "" else sprintf(
      paste0(
        "No generation satisfied the population acceptance rule for ",
        "scenario \"%s\" within %d attempts (PSD failures: %d). The ",
        "returned object contains the last (unaccepted) generation."
      ),
      scenario, max_attempts, psd_failures
    )
  )

  list(
    Y = Y,
    W = W,
    G = G,
    K_bg = K_bg,
    focal_region = marker_names[region_cols],
    causal_index = causal_index,
    truth = truth,
    latent = if (isTRUE(return_latent)) list(U = U, E = E) else NULL,
    generator = list(
      seed = seed,
      attempts = attempt,
      settings = list(
        n = n, m = m, p = p, maf_range = maf_range, h2 = h2,
        experiment = experiment, scenario = scenario,
        n_signals = q, locus_pve = locus_pve,
        secondary_signal_pve = secondary_signal_pve,
        local_ld = local_ld, target_r2 = r2_target,
        local_region_size = local_region_size,
        tolerance = tolerance, target_loss = rho_star,
        target_representing_set = S_star,
        correlation = correlation, structured = structured,
        n_groups = n_groups, fst = fst, truth_margin = truth_margin,
        max_attempts = max_attempts, return_latent = return_latent
      )
    ),
    status = status,
    diagnostics = list(
      attempts = attempt,
      psd_failures = psd_failures,
      accepted = accepted,
      mean_diag_K = mean(diag(K_bg)),
      trace_K_over_n = sum(diag(K_bg)) / n,
      K_eigen_sd = stats::sd(lambda_K),
      K_rank = K_rank,
      min_eigen_Sigma_G_bg = min_eigen_bg,
      truth_margin_achieved = truth_tabs$margin_achieved
    )
  )
}

# ---- internal helpers ---------------------------------------------------------

#' Check a positive integer-like scalar argument
#' @keywords internal
.check_count <- function(x, name, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != round(x) || x < min) {
    .stop_invalid_input("%s must be a single integer >= %d.", name, min)
  }
  invisible(TRUE)
}

#' Exchangeable correlation matrix
#' @keywords internal
.exchangeable_cor <- function(m, r) {
  matrix(r, m, m) + diag(1 - r, m)
}

#' Block correlation matrix (within-block 0.6, between-block 0.2)
#' @keywords internal
.block_cor <- function(m, within = 0.6, between = 0.2) {
  R <- matrix(between, m, m)
  block_id <- (seq_len(m) + 1L) %/% 2L
  same <- outer(block_id, block_id, `==`)
  R[same] <- within
  diag(R) <- 1
  R
}

#' Validate a correlation matrix argument
#' @keywords internal
.check_correlation_matrix <- function(R, m, name) {
  if (!is.matrix(R) || !is.numeric(R) || nrow(R) != m || ncol(R) != m) {
    .stop_invalid_input("%s must be a numeric %d x %d matrix.", name, m, m)
  }
  if (anyNA(R) || max(abs(R - t(R))) > 1e-10) {
    .stop_invalid_input(
      "%s must be symmetric and free of missing values.", name
    )
  }
  if (max(abs(diag(R) - 1)) > 1e-10) {
    .stop_invalid_input("%s must have unit diagonal.", name)
  }
  if (.min_eigen_sym(R) <= 0) {
    .stop_invalid_input("%s must be positive definite.", name)
  }
  invisible(TRUE)
}

#' Build a covariance matrix from marginal variances and a correlation matrix
#' @keywords internal
.cov_from_h2 <- function(var, R) {
  R * outer(sqrt(var), sqrt(var))
}

#' Smallest eigenvalue of a symmetric matrix
#' @keywords internal
.min_eigen_sym <- function(A) {
  min(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
}

#' Symmetric square root of a PSD matrix
#' @keywords internal
.symmetric_psd_sqrt <- function(Sigma) {
  eig <- eigen(Sigma, symmetric = TRUE)
  sweep(eig$vectors, 2L, sqrt(pmax(eig$values, 0)), `*`) %*% t(eig$vectors)
}

#' True conditional contrasts from a phenotypic covariance matrix
#'
#' Legacy helper still used by the `fit_mt_null()` compatibility
#' fields (`compute_gamma`).
#' @keywords internal
.contrasts_from_cov <- function(Sigma_P, trait_names = colnames(Sigma_P)) {
  m <- ncol(Sigma_P)
  gamma <- matrix(0, nrow = m, ncol = m - 1L,
                  dimnames = list(trait_names, NULL))
  C <- diag(m)
  dimnames(C) <- list(trait_names, trait_names)
  conditional_variance <- numeric(m)
  for (i in seq_len(m)) {
    inv <- .safe_inverse(Sigma_P[-i, -i, drop = FALSE])
    if (inv$status == "failed") {
      stop(sprintf("Failed to invert the phenotypic covariance submatrix for trait %d.", i),
           call. = FALSE)
    }
    g <- drop(inv$inverse %*% Sigma_P[-i, i])
    gamma[i, ] <- g
    C[-i, i] <- -g
    conditional_variance[i] <- Sigma_P[i, i] - sum(Sigma_P[i, -i] * g)
  }
  list(gamma = gamma, C = C, conditional_variance = conditional_variance)
}

#' Generate genotype dosages with optional group structure
#' @keywords internal
.sim_genotypes <- function(n, maf, structured, n_groups, fst) {
  p <- length(maf)
  if (!structured) {
    return(matrix(
      stats::rbinom(n * p, size = 2, prob = rep(maf, each = n)),
      nrow = n, ncol = p
    ))
  }
  sizes <- rep(n %/% n_groups, n_groups)
  remainder <- n %% n_groups
  if (remainder > 0L) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
  group <- rep(seq_len(n_groups), times = sizes)
  shape1 <- maf * (1 - fst) / fst
  shape2 <- (1 - maf) * (1 - fst) / fst
  p_group <- matrix(
    stats::rbeta(p * n_groups, rep(shape1, times = n_groups),
                 rep(shape2, times = n_groups)),
    nrow = n_groups, byrow = TRUE
  )
  matrix(
    stats::rbinom(n * p, size = 2, prob = p_group[group, ]),
    nrow = n, ncol = p
  )
}

#' Generate a focal-region genotype block with target local LD
#'
#' Latent-factor Gaussian copula: marker j's latent variable is
#' a * F + sqrt(1 - a^2) * e_j with F shared across the region; the
#' factor loading `a` is calibrated once (deterministically, within
#' the RNG stream) so the realized mean pairwise dosage r^2 matches
#' `target_r2` up to discretization noise.
#'
#' @param n Individuals.
#' @param size Markers in the region.
#' @param maf_vec Per-marker MAFs.
#' @param target_r2 Target mean pairwise r^2 (0 gives independent
#'   dosages).
#' @return List with `G` (n x size dosages), `latent_r2` (model
#'   correlation used for Sigma_X), `mean_r2`, `max_r2` (realized).
#' @keywords internal
.sim_ld_region <- function(n, size, maf_vec, target_r2) {
  mean_off_r2 <- function(G) {
    cr2 <- stats::cor(G)^2
    mean(cr2[upper.tri(cr2)])
  }
  if (target_r2 <= 0) {
    G <- matrix(
      stats::rbinom(n * size, 2, rep(maf_vec, each = n)),
      n, size
    )
    return(list(G = G, latent_r2 = 0,
                mean_r2 = 0, max_r2 = 0))
  }
  # calibrate the latent loading on a FIXED latent noise draw, so the
  # response is a smooth deterministic function of `a` (bisection)
  F_lat <- stats::rnorm(n)
  E_lat <- matrix(stats::rnorm(n * size), n, size)
  build <- function(a) {
    Z <- matrix(a * F_lat, n, size) + sqrt(1 - a^2) * E_lat
    G <- matrix(0L, n, size)
    for (j in seq_len(size)) {
      f <- maf_vec[j]
      t1 <- stats::qnorm((1 - f)^2)
      t2 <- stats::qnorm(1 - f^2)
      G[, j] <- as.integer(Z[, j] > t1) + as.integer(Z[, j] > t2)
    }
    G
  }
  r2_of <- function(a) mean_off_r2(build(a))
  lo <- 0; hi <- 0.9999
  if (r2_of(hi) <= target_r2) {
    a <- hi                      # cannot reach the target; keep the max
  } else {
    for (iter in 1:40) {
      mid <- (lo + hi) / 2
      if (r2_of(mid) < target_r2) lo <- mid else hi <- mid
    }
    a <- (lo + hi) / 2
  }
  G <- build(a)
  cr2 <- stats::cor(G)^2
  off <- cr2[upper.tri(cr2)]
  list(
    G = G,
    latent_r2 = a^2,
    mean_r2 = mean(off),
    max_r2 = max(off)
  )
}

#' Per-signal marginal PVE scaler (single-signal frozen formula)
#'
#' `c = sqrt(target / (v * mean(beta_row[nz]^2)))` with
#' `nz = beta_row != 0`; identical to the legacy single-signal scaler.
#' @keywords internal
.scale_signal_pve <- function(beta_row, v, target_pve) {
  nz <- beta_row != 0
  if (!any(nz) || target_pve <= 0) return(1)
  sqrt(target_pve / (v * mean(beta_row[nz]^2)))
}

#' Scenario registry for the three formal simulation families
#'
#' @param experiment Experiment name.
#' @param scenario Scenario name.
#' @return A list describing the scenario.
#' @keywords internal
.scenario_spec <- function(experiment, scenario) {
  pat <- function(...) c(...)          # shorthand
  registry <- list(
    signal_resolution = list(
      null = list(q = 0L, kind = "null", min_m = 1L),
      single_multi_trait = list(
        q = 1L, kind = "fixed_dirs", min_m = 2L,
        dirs = function(m) list(c(1, 0.8, rep(0, m - 2L)))
      ),
      two_linked_trait_specific = list(
        q = 2L, kind = "fixed_dirs", min_m = 2L, local_ld = "moderate",
        dirs = function(m) list(c(1, rep(0, m - 1L)),
                                c(0, 1, rep(0, m - 2L)))
      ),
      two_heterogeneous = list(
        q = 2L, kind = "fixed_dirs", min_m = 3L,
        dirs = function(m) list(c(1, 0.8, rep(0, m - 2L)),
                                c(0, 0, 1, rep(0, m - 3L)))
      )
    ),
    trait_representation = list(
      trait_specific = list(
        q = 1L, kind = "fixed_dirs", min_m = 1L,
        dirs = function(m) list(c(1, rep(0, m - 1L)))
      ),
      two_trait_concordant = list(
        q = 1L, kind = "fixed_dirs", min_m = 2L,
        dirs = function(m) list(c(1, 0.8, rep(0, m - 2L)))
      ),
      two_trait_antagonistic = list(
        q = 1L, kind = "fixed_dirs", min_m = 2L,
        dirs = function(m) list(c(1, -0.8, rep(0, m - 2L)))
      ),
      broad_concordant = list(
        q = 1L, kind = "fixed_dirs", min_m = 2L,
        dirs = function(m) list(1 - 0.2 * seq(0, length.out = m))
      ),
      highly_representable = list(
        q = 1L, kind = "target_rho", min_m = 2L, accept = "R1",
        scale_match = TRUE,
        default_S_star = function(m, tn) tn[1L]
      ),
      partially_representable = list(
        q = 1L, kind = "target_rho", min_m = 3L, accept = "R2",
        scale_match = TRUE,
        default_S_star = function(m, tn) tn[1:2]
      ),
      strongly_nonredundant = list(
        q = 1L, kind = "nonredundant", min_m = 2L, accept = "R3",
        scale_match = TRUE
      )
    ),
    end_to_end = list(
      single_highly_representable = list(
        q = 1L, kind = "target_rho", min_m = 2L, accept = "R1",
        scale_match = TRUE,
        default_S_star = function(m, tn) tn[1L]
      ),
      single_nonredundant = list(
        q = 1L, kind = "nonredundant", min_m = 2L, accept = "R3",
        scale_match = TRUE
      ),
      linked_pseudo_multitrait = list(
        q = 2L, kind = "fixed_dirs", min_m = 2L, local_ld = "high",
        dirs = function(m) list(c(1, rep(0, m - 1L)),
                                c(0, 1, rep(0, m - 2L)))
      ),
      mixed_multisignal = list(
        q = 2L, kind = "fixed_dirs", min_m = 3L, local_ld = "moderate",
        dirs = function(m) list(c(1, 0.8, rep(0, m - 2L)),
                                c(0, 0, 1, rep(0, m - 3L)))
      )
    )
  )
  fam <- registry[[experiment]]
  if (is.null(fam)) {
    .stop_invalid_input("unknown experiment \"%s\".", experiment)
  }
  spec <- fam[[scenario]]
  if (is.null(spec)) {
    .stop_invalid_input(
      "unknown scenario \"%s\" for experiment \"%s\". Available: %s.",
      scenario, experiment, paste(names(fam), collapse = ", ")
    )
  }
  spec$experiment <- experiment
  spec$scenario <- scenario
  if (is.null(spec$accept)) spec$accept <- "none"
  if (is.null(spec$scale_match)) spec$scale_match <- FALSE
  spec
}

#' Propose the q x m effect matrix for one attempt
#' @keywords internal
.propose_B_Q <- function(spec, m, trait_names, Sigma_P, S_star,
                         rho_star) {
  q <- spec$q
  B_Q <- matrix(0, nrow = max(q, 1L), ncol = m,
                dimnames = list(NULL, trait_names))
  delta <- NA_real_
  if (q == 0L || spec$kind == "null") {
    return(list(B_Q = B_Q[0, , drop = FALSE], delta = delta))
  }
  if (spec$kind == "fixed_dirs") {
    dirs <- spec$dirs(m)
    for (i in seq_len(q)) B_Q[i, ] <- dirs[[i]]
    return(list(B_Q = B_Q, delta = delta))
  }
  if (spec$kind == "nonredundant") {
    beta <- sign(stats::rnorm(m)) * (1 + 0.3 * abs(stats::rnorm(m)))
    B_Q[1L, ] <- beta
    return(list(B_Q = B_Q, delta = delta))
  }
  # target_rho closed form
  A <- trait_names
  beta_S <- stats::rnorm(length(S_star))
  w <- stats::rnorm(m - length(S_star))
  out <- .generate_target_rho_beta(S_star, beta_S, Sigma_P[A, A],
                                   rho_star, w)
  B_Q[1L, ] <- out$beta
  delta <- out$delta
  list(B_Q = B_Q, delta = delta)
}

#' Closed-form target-rho effect construction (Methods)
#'
#' With S* the target representing set and R = A \ S*:
#' \deqn{\beta_R = \Gamma^\top \beta_S + \delta d,\qquad
#'       d^\top \Omega^{-1} d = 1,\qquad
#'       \delta = \sqrt{\rho^\star A_0 / (1 - \rho^\star)}.}
#' The sign convention is always "+ delta d".
#'
#' @param S_star Target representing set (trait names).
#' @param beta_S Effects on S_star.
#' @param Sigma_AA Population covariance of the candidate traits.
#' @param rho_star Target loss in [0, 1).
#' @param w Random direction seed (length |R|).
#' @return List with `beta` (full length-m named vector) and `delta`.
#' @keywords internal
.generate_target_rho_beta <- function(S_star, beta_S, Sigma_AA,
                                      rho_star, w) {
  A <- colnames(Sigma_AA)
  S <- match(S_star, A)
  R <- setdiff(seq_along(A), S)
  if (length(R) == 0L) {
    .stop_invalid_input(
      "the target-rho construction requires a non-empty complement."
    )
  }
  if (!is.numeric(rho_star) || rho_star < 0 || rho_star >= 1) {
    .stop_invalid_input("rho_star must be in [0, 1).")
  }
  inv_SS <- .safe_inverse(Sigma_AA[S, S, drop = FALSE])
  Gamma <- inv_SS$inverse %*% Sigma_AA[S, R, drop = FALSE]
  Omega <- Sigma_AA[R, R, drop = FALSE] -
    crossprod(Sigma_AA[S, R, drop = FALSE],
              inv_SS$inverse %*% Sigma_AA[S, R, drop = FALSE])
  Omega <- (Omega + t(Omega)) / 2
  d <- drop(Omega %*% w)
  d <- d / sqrt(sum(w * d))
  A0 <- sum(beta_S * (inv_SS$inverse %*% beta_S))
  delta <- if (rho_star == 0) {
    0
  } else {
    sqrt(rho_star / (1 - rho_star) * A0)
  }
  beta <- numeric(length(A))
  beta[S] <- beta_S
  beta[R] <- drop(crossprod(Gamma, beta_S)) + delta * d
  names(beta) <- A
  list(beta = beta, delta = delta)
}

#' Exact truth representation map for one signal
#'
#' Enumerates all subsets of the signal's candidate set with the
#' unique Stage-1 loss path and extracts Rep/Irr from the map only.
#'
#' @param beta Named effect vector of the signal (length m).
#' @param Sigma_P Reference covariance.
#' @param tolerance Primary tolerance.
#' @param tolerances All tolerances for extraction.
#' @param ids List with `locus_id`, `signal_id`, `representative_snp`.
#' @return List with `subset_table`, `representation_map`,
#'   `minimum_representative_sets`, `irreducible_modules`,
#'   `tolerance_path`, `margin_achieved`.
#' @keywords internal
.compute_truth_representation_map <- function(beta, Sigma_P, tolerance,
                                              tolerances, ids) {
  trait_names <- names(beta)
  A <- trait_names[beta != 0]
  empty <- .empty_truth_tables(ids)
  k <- length(A)
  if (k == 0L) return(empty)

  Sigma_AA <- Sigma_P[A, A, drop = FALSE]
  beta_A <- beta[A]
  enum <- .enumerate_subsets(A, "all")
  rows <- vector("list", length(enum$sets))
  for (si in seq_along(enum$sets)) {
    out <- .compute_subset_loss(beta_A, Sigma_AA, enum$sets[[si]],
                                trait_names = A)
    rows[[si]] <- data.frame(
      locus_id = ids$locus_id, signal_id = ids$signal_id,
      representative_snp = ids$representative_snp,
      marker_id = ids$representative_snp,
      set_id = si,
      representing_set = I(list(out$representing_set)),
      representing_key = out$representing_key,
      complement_set = I(list(out$complement_set)),
      complement_key = out$complement_key,
      set_size = length(out$representing_set),
      conditional_effect = I(list(out$eta)),
      full_qform = out$full_qform,
      subset_qform = out$subset_qform,
      residual_qform = out$residual_qform,
      representation_loss = out$representation_loss,
      feasible_primary = !is.na(out$representation_loss) &&
        out$representation_loss <= tolerance + 1e-10,
      status = out$status$code,
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, rows)
  rep_map <- tab[, c("locus_id", "signal_id", "representing_set",
                     "representing_key", "complement_set",
                     "complement_key", "set_size",
                     "representation_loss", "feasible_primary")]

  reps_list <- mods_list <- path_list <- list()
  for (tol in tolerances) {
    r <- .extract_minimum_representative_sets(tab, tol)
    mo <- .extract_irreducible_modules(tab, tol)
    if (nrow(r) > 0L) {
      reps_list[[length(reps_list) + 1L]] <- cbind(
        locus_id = ids$locus_id, signal_id = ids$signal_id,
        r[, setdiff(names(r), "status")], stringsAsFactors = FALSE)
    }
    if (nrow(mo) > 0L) {
      mods_list[[length(mods_list) + 1L]] <- cbind(
        locus_id = ids$locus_id, signal_id = ids$signal_id,
        mo[, setdiff(names(mo), "status")], stringsAsFactors = FALSE)
    }
    path_list[[length(path_list) + 1L]] <- data.frame(
      locus_id = ids$locus_id, signal_id = ids$signal_id,
      tolerance = tol,
      minimum_set_size = if (nrow(r) > 0L) r$set_size[1L] else NA_integer_,
      n_minimum_sets = nrow(r), n_irreducible_modules = nrow(mo),
      stringsAsFactors = FALSE
    )
  }
  bind <- function(lst, proto) {
    if (length(lst) == 0L) return(proto)
    out <- do.call(rbind, lst)
    rownames(out) <- NULL
    out
  }
  losses <- tab$representation_loss
  list(
    subset_table = tab,
    representation_map = rep_map,
    minimum_representative_sets = bind(reps_list,
                                       empty$minimum_representative_sets),
    irreducible_modules = bind(mods_list, empty$irreducible_modules),
    tolerance_path = bind(path_list, empty$tolerance_path),
    margin_achieved = if (all(is.finite(losses))) {
      min(abs(losses - tolerance))
    } else {
      NA_real_
    }
  )
}

#' Empty truth tables prototype
#' @keywords internal
.empty_truth_tables <- function(ids) {
  list(
    subset_table = data.frame(
      locus_id = character(), signal_id = character(),
      representative_snp = character(), marker_id = character(),
      set_id = integer(),
      representing_set = I(list()), representing_key = character(),
      complement_set = I(list()), complement_key = character(),
      set_size = integer(), conditional_effect = I(list()),
      full_qform = numeric(), subset_qform = numeric(),
      residual_qform = numeric(), representation_loss = numeric(),
      feasible_primary = logical(), status = character(),
      stringsAsFactors = FALSE
    ),
    representation_map = data.frame(
      locus_id = character(), signal_id = character(),
      representing_set = I(list()), representing_key = character(),
      complement_set = I(list()), complement_key = character(),
      set_size = integer(), representation_loss = numeric(),
      feasible_primary = logical(), stringsAsFactors = FALSE
    ),
    minimum_representative_sets = data.frame(
      locus_id = character(), signal_id = character(),
      tolerance = numeric(), solution_id = integer(),
      trait_set = I(list()), trait_key = character(),
      set_size = integer(), representation_loss = numeric(),
      n_tied_solutions = integer(), stringsAsFactors = FALSE
    ),
    irreducible_modules = data.frame(
      locus_id = character(), signal_id = character(),
      tolerance = numeric(), module_id = integer(),
      trait_set = I(list()), trait_key = character(),
      module_size = integer(), complement_set = I(list()),
      complement_key = character(), complement_loss = numeric(),
      stringsAsFactors = FALSE
    ),
    tolerance_path = data.frame(
      locus_id = character(), signal_id = character(),
      tolerance = numeric(), minimum_set_size = integer(),
      n_minimum_sets = integer(), n_irreducible_modules = integer(),
      stringsAsFactors = FALSE
    ),
    margin_achieved = NA_real_
  )
}

#' Truth tables for all signals of a replicate
#' @keywords internal
.truth_tables_per_signal <- function(B_Q, Sigma_P, tolerance, tolerances,
                                     locus_id, chr, pos, causal_ids,
                                     trait_names) {
  q <- nrow(B_Q)
  if (q == 0L) {
    return(.empty_truth_tables(list(locus_id = character(),
                                    signal_id = character(),
                                    representative_snp = character())))
  }
  all <- lapply(seq_len(q), function(i) {
    ids <- list(locus_id = locus_id,
                signal_id = paste0(locus_id, "::S", i),
                representative_snp = causal_ids[i])
    .compute_truth_representation_map(
      B_Q[i, ], Sigma_P, tolerance, tolerances, ids
    )
  })
  merge_field <- function(field) {
    out <- do.call(rbind, lapply(all, `[[`, field))
    rownames(out) <- NULL
    out
  }
  losses <- merge_field("subset_table")$representation_loss
  list(
    subset_table = merge_field("subset_table"),
    representation_map = merge_field("representation_map"),
    minimum_representative_sets = merge_field("minimum_representative_sets"),
    irreducible_modules = merge_field("irreducible_modules"),
    tolerance_path = merge_field("tolerance_path"),
    margin_achieved = if (all(is.finite(losses))) {
      min(abs(losses - tolerance))
    } else {
      NA_real_
    }
  )
}

#' Validate a representation architecture against its truth map
#'
#' @param architecture `"R1"`, `"R2"` or `"R3"`.
#' @param tabs One signal's truth tables.
#' @param tolerance Primary tolerance.
#' @param margin Required distance from the tolerance boundary.
#' @return List with `architecture`, `accepted`,
#'   `min_rep_cardinality`, `min_proper_subset_rho`,
#'   `target_subset_rho`, `diagnostics`.
#' @keywords internal
.validate_representation_architecture <- function(architecture, tabs,
                                                  tolerance, margin) {
  tab <- tabs$subset_table
  k <- max(tab$set_size)
  losses <- tab$representation_loss
  reps <- tabs$minimum_representative_sets
  reps <- reps[reps$tolerance == tolerance, ]
  min_card <- if (nrow(reps) > 0L) min(reps$set_size) else NA_integer_
  proper <- losses[tab$set_size < k]
  min_proper <- if (length(proper) > 0L) min(proper) else NA_real_
  singles <- losses[tab$set_size == 1L]
  pairs <- losses[tab$set_size == 2L]
  margin_ok <- margin <= 0 ||
    all(abs(losses - tolerance) >= margin)
  accepted <- switch(architecture,
    R1 = identical(min_card, 1L) && margin_ok,
    R2 = identical(min_card, 2L) &&
      all(singles > tolerance) && any(pairs <= tolerance) && margin_ok,
    R3 = all(proper > tolerance) && margin_ok,
    FALSE
  )
  list(
    architecture = architecture,
    accepted = isTRUE(accepted),
    min_rep_cardinality = min_card,
    min_proper_subset_rho = min_proper,
    target_subset_rho = NA_real_,
    diagnostics = list(
      tolerance = tolerance, margin = margin,
      min_singleton_rho = if (length(singles)) min(singles) else NA_real_,
      min_pair_rho = if (length(pairs)) min(pairs) else NA_real_
    )
  )
}

#' Scenario acceptance rule (population truth only)
#' @keywords internal
.accept_scenario <- function(spec, tabs, tolerance, margin) {
  if (spec$accept == "none") return(TRUE)
  if (nrow(tabs$subset_table) == 0L) return(FALSE)
  v <- .validate_representation_architecture(spec$accept, tabs,
                                             tolerance, margin)
  v$accepted
}

#' Frozen descriptive direction rule for an effect vector
#' @keywords internal
.direction_of <- function(effects) {
  B <- length(effects)
  if (B == 0L) return("not_applicable")
  if (B == 1L) return("single_trait")
  sgn <- sign(effects)
  if (all(sgn == sgn[1L])) return("concordant")
  if (B == 2L) return("antagonistic")
  "mixed"
}
