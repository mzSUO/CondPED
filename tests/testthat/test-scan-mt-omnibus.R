# Tests for scan_mt_omnibus() (S3)

# ---- helpers -----------------------------------------------------------------

# Explicit dense construction of V, X0, P, Z for a single marker.
dense_scan_components <- function(null_fit, x, Y) {
  n <- nrow(Y)
  m <- ncol(Y)
  # Recover the ORIGINAL design from its rotation: W = U %*% W_tilde.
  W <- null_fit$rotation$U %*% null_fit$rotation$W_tilde
  SG <- null_fit$Sigma_G
  SE <- null_fit$Sigma_E
  K <- null_fit$rotation$U %*% diag(null_fit$rotation$lambda) %*% t(null_fit$rotation$U)
  V <- kronecker(SG, K) + kronecker(SE, diag(n))
  X0 <- kronecker(diag(m), W)
  Z <- kronecker(diag(m), matrix(x, ncol = 1))
  y <- as.vector(Y)
  Vc <- chol(V)
  Vinv <- chol2inv(Vc)
  XtVinvX <- crossprod(X0, Vinv %*% X0)
  G <- solve(XtVinvX)
  P <- Vinv - Vinv %*% X0 %*% G %*% crossprod(X0, Vinv)
  U <- crossprod(Z, P %*% y)
  J <- crossprod(Z, P %*% Z)
  list(U = drop(U), J = J)
}

make_small_null <- function(n = 12, m = 2, seed = 1) {
  set.seed(seed)
  Z <- matrix(rnorm(n * 200), n, 200)
  K <- tcrossprod(Z) / 200
  K <- K / mean(diag(K))
  SG <- matrix(0.3, m, m) + diag(0.7, m)
  SE <- matrix(-0.1, m, m) + diag(1.1, m)
  SE <- (SE + t(SE)) / 2
  if (min(eigen(SE, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    SE <- matrix(0.1, m, m) + diag(0.9, m)
  }
  W <- cbind(1, rnorm(n))
  B <- matrix(rnorm(2 * m), 2, m)
  u <- t(chol(K)) %*% matrix(rnorm(n * m), n, m) %*% chol(SG)
  Y <- W %*% B + u + matrix(rnorm(n * m), n, m) %*% chol(SE)
  fit <- fit_mt_null(Y, W = W[, 2, drop = FALSE], K = K, n_starts = 2L,
                     control = list(maxit = 200))
  fit
}

# ---- m = 1 degenerate --------------------------------------------------------

test_that("m = 1 score matches the direct single-trait GLS score", {
  set.seed(11)
  n <- 60
  Z <- matrix(rnorm(n * 200), n, 200)
  K <- tcrossprod(Z) / 200
  K <- K / mean(diag(K))
  W <- cbind(1, rnorm(n))
  Y <- matrix(W %*% c(0.5, -0.8) +
                drop(t(chol(K)) %*% rnorm(n)) * sqrt(0.6) +
                rnorm(n, sd = sqrt(0.4)), ncol = 1)
  fit <- fit_mt_null(Y, W = W[, 2, drop = FALSE], K = K, n_starts = 2L,
                     control = list(maxit = 200))

  set.seed(12)
  x <- stats::rbinom(n, 2, 0.3)
  scan <- scan_mt_omnibus(fit, matrix(x, ncol = 1), marker_ids = "M1",
                          return_effects = TRUE)

  # Direct single-trait GLS
  V <- fit$Sigma_G[1, 1] * K + fit$Sigma_E[1, 1] * diag(n)
  Vinv <- solve(V)
  XtVinvX <- crossprod(W, Vinv %*% W)
  P <- Vinv - Vinv %*% W %*% solve(XtVinvX) %*% crossprod(W, Vinv)
  U_d <- sum(x * (P %*% Y))
  J_d <- sum(x * (P %*% x))
  Q_d <- U_d^2 / J_d
  p_d <- stats::pchisq(Q_d, df = 1, lower.tail = FALSE)
  beta_d <- U_d / J_d
  se_d <- sqrt(1 / J_d)

  expect_equal(scan$omnibus$Q[1], Q_d, tolerance = 1e-7)
  expect_equal(scan$omnibus$df[1], 1L)
  expect_equal(scan$omnibus$p_value[1], p_d, tolerance = 1e-7)
  eff <- scan$effects$effects_long
  expect_equal(eff$beta[1], beta_d, tolerance = 1e-7)
  expect_equal(eff$se[1], se_d, tolerance = 1e-7)
})

# ---- dense对照 ---------------------------------------------------------------

test_that("rotated block computation matches explicit dense V/P/Z", {
  fit <- make_small_null(n = 12, m = 2, seed = 21)
  n <- nrow(fit$rotation$Y_tilde)
  Y <- fit$rotation$U %*% fit$rotation$Y_tilde
  set.seed(22)
  x <- stats::rbinom(n, 2, 0.3)

  rot <- fit$rotation
  x_tilde <- crossprod(rot$U, x)
  AM_arr <- CondPED:::.precompute_AM(rot)
  Ar <- CondPED:::.precompute_Ar(rot)
  A_arr <- rot$Vinv
  block <- CondPED:::.gls_block_components(x_tilde, AM_arr, Ar, A_arr,
                                           rot$XtVinvX_inv,
                                           sqrt(.Machine$double.eps))

  dense <- dense_scan_components(fit, x, Y)
  expect_equal(block$U, dense$U, tolerance = 1e-7)
  expect_equal(unname(block$J), unname(dense$J), tolerance = 1e-7)

  scan <- scan_mt_omnibus(fit, matrix(x, ncol = 1), marker_ids = "M1")
  Q_d <- sum(dense$U * solve(dense$J, dense$U))
  expect_equal(scan$omnibus$Q[1], Q_d, tolerance = 1e-7)
})

# ---- effect identities -------------------------------------------------------

test_that("effect estimates satisfy beta = J^+ U and covariance = J^+", {
  fit <- make_small_null(n = 30, m = 2, seed = 31)
  set.seed(32)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 5, 2, 0.3), ncol = 5)
  scan <- scan_mt_omnibus(fit, G, return_effects = TRUE)
  est <- estimate_mt_effects(fit, G, loci = 1:5)

  expect_equal(scan$effects$beta, est$beta, tolerance = 1e-10)

  for (i in seq_len(5)) {
    ef <- est$effects_long[est$effects_long$marker_id == paste0("M", i), ]
    expect_equal(ef$se^2, pmax(diag(est$covariance[, , i]), 0),
                 tolerance = 1e-10)
    expect_equal(ef$z, ef$beta / ef$se, tolerance = 1e-10)
    expect_equal(ef$p_value,
                 stats::pchisq(ef$z^2, df = 1, lower.tail = FALSE),
                 tolerance = 1e-10)
  }
})

# ---- singular J --------------------------------------------------------------

test_that("rank-zero J returns NA and df = 0", {
  fit <- make_small_null(n = 20, m = 2, seed = 41)
  # An all-zero marker has x_tilde = 0, hence J = 0 exactly.
  x <- rep(0, nrow(fit$rotation$Y_tilde))
  est <- estimate_mt_effects(fit, matrix(x, ncol = 1), loci = 1L)
  expect_true(all(is.na(est$beta)))
  expect_true(all(is.na(est$effects_long$se)))

  # The same rank-zero block through the internal component builder
  rot <- fit$rotation
  AM_arr <- CondPED:::.precompute_AM(rot)
  Ar <- CondPED:::.precompute_Ar(rot)
  block <- CondPED:::.gls_block_components(
    numeric(nrow(rot$Y_tilde)), AM_arr, Ar, rot$Vinv, rot$XtVinvX_inv,
    sqrt(.Machine$double.eps))
  expect_equal(block$rank, 0L)
})

test_that("rank-one J reports df = rank(J)", {
  # Construct an A array whose first block is rank 1, then isolate it.
  # AM_arr and Ar are zero so the correction term vanishes.
  m <- 2L
  qm <- 2L
  A1 <- matrix(c(1, 1, 1, 1), m, m)  # rank 1
  A_arr <- array(0, dim = c(m, m, 4))
  A_arr[, , 1] <- A1
  AM_arr <- array(0, dim = c(m, qm, 4))
  Ar <- matrix(rnorm(m * 4), m, 4)
  x_tilde <- c(1, 0, 0, 0)
  block <- CondPED:::.gls_block_components(x_tilde, AM_arr, Ar, A_arr,
                                           diag(qm),
                                           sqrt(.Machine$double.eps))
  expect_equal(block$rank, 1L)
  expect_equal(block$status, "rank_deficient")
})

# ---- marker filtering --------------------------------------------------------

test_that("markers are filtered with reasons and retained in output", {
  fit <- make_small_null(n = 30, m = 2, seed = 51)
  n <- nrow(fit$rotation$Y_tilde)
  G <- matrix(0, n, 6)
  G[, 1] <- stats::rbinom(n, 2, 0.3)            # ok
  G[, 2] <- 0L                                   # monomorphic / all-zero
  G[, 3] <- rep(c(0, 2), length.out = n)         # maf = 0.5 exactly -> kept
  G[, 4] <- c(1L, rep(0L, n - 1L))               # low maf (1/60 < 0.05)
  G[, 5] <- stats::rbinom(n, 2, 0.3)             # varied dosages...
  G[1:5, 5] <- NA                                # ...then 5 NA but still ok
  G[, 6] <- NA                                   # all NA
  colnames(G) <- c("ok", "mono", "half", "low", "partial", "allna")

  scan <- scan_mt_omnibus(fit, G)
  expect_equal(scan$diagnostics$n_tested, 3L)    # ok + half + partial
  expect_equal(scan$diagnostics$n_filtered, 3L)  # mono, low, allna

  expect_equal(scan$omnibus$status[scan$omnibus$marker_id == "mono"],
               "filtered")
  expect_equal(scan$omnibus$filter_reason[scan$omnibus$marker_id == "mono"],
               "monomorphic")
  expect_equal(scan$omnibus$filter_reason[scan$omnibus$marker_id == "low"],
               "low_maf")
  expect_equal(scan$omnibus$filter_reason[scan$omnibus$marker_id == "allna"],
               "all_missing")
  # maf = 0.5 marker is kept
  expect_equal(scan$omnibus$status[scan$omnibus$marker_id == "half"], "ok")
  # partial NA keeps n_eff
  expect_equal(scan$omnibus$n_eff[scan$omnibus$marker_id == "partial"], n - 5L)
})

# ---- consistency scan return_effects vs estimate ------------------------------

test_that("scan(return_effects = TRUE) matches estimate_mt_effects", {
  fit <- make_small_null(n = 40, m = 3, seed = 61)
  set.seed(62)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 8, 2, 0.3), ncol = 8)
  scan <- scan_mt_omnibus(fit, G, return_effects = TRUE)
  est <- estimate_mt_effects(fit, G, loci = 1:8)
  expect_equal(scan$effects$beta, est$beta, tolerance = 1e-10)
  expect_equal(scan$effects$effects_long, est$effects_long, tolerance = 1e-10)
  expect_equal(scan$effects$covariance, est$covariance, tolerance = 1e-10)
})

# ---- chunk sizes --------------------------------------------------------------

test_that("chunk_size does not change results", {
  fit <- make_small_null(n = 25, m = 2, seed = 71)
  set.seed(72)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 7, 2, 0.3), ncol = 7)
  s1 <- scan_mt_omnibus(fit, G, chunk_size = 1L)
  s2 <- scan_mt_omnibus(fit, G, chunk_size = 2L)
  s_all <- scan_mt_omnibus(fit, G, chunk_size = 100L)
  expect_equal(s1$omnibus$Q, s2$omnibus$Q, tolerance = 1e-10)
  expect_equal(s1$omnibus$Q, s_all$omnibus$Q, tolerance = 1e-10)
})

# ---- null_fit validation ------------------------------------------------------

test_that("invalid null_fit objects raise condped_invalid_input", {
  fit <- make_small_null(n = 20, m = 2, seed = 81)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 3, 2, 0.3), ncol = 3)

  bad1 <- fit
  class(bad1) <- "list"
  expect_error(scan_mt_omnibus(bad1, G), class = "condped_invalid_input")

  bad2 <- fit
  bad2$status$ok <- FALSE
  expect_error(scan_mt_omnibus(bad2, G), class = "condped_invalid_input")

  bad3 <- fit
  bad3$rotation <- NULL
  expect_error(scan_mt_omnibus(bad3, G), class = "condped_invalid_input")

  expect_error(estimate_mt_effects(bad1, G, loci = 1),
               class = "condped_invalid_input")
})

# ---- marker_ids length mismatch ----------------------------------------------

test_that("marker_ids length mismatch is rejected", {
  fit <- make_small_null(n = 20, m = 2, seed = 91)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 3, 2, 0.3), ncol = 3)
  expect_error(scan_mt_omnibus(fit, G, marker_ids = c("a", "b")),
               class = "condped_invalid_input")
  expect_error(estimate_mt_effects(fit, G, loci = 1:2,
                                   marker_ids = c("a", "b", "c", "d")),
               class = "condped_invalid_input")
})

# ---- return flags --------------------------------------------------------------

test_that("return_score and return_effects control optional outputs", {
  fit <- make_small_null(n = 20, m = 2, seed = 101)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 3, 2, 0.3), ncol = 3)
  s0 <- scan_mt_omnibus(fit, G)
  s1 <- scan_mt_omnibus(fit, G, return_score = TRUE)
  s2 <- scan_mt_omnibus(fit, G, return_effects = TRUE)
  expect_null(s0$score)
  expect_null(s0$effects)
  expect_equal(length(s1$score), 3L)
  expect_type(s1$score[[1]], "double")
  expect_equal(length(s2$effects$beta), 6L)
})
