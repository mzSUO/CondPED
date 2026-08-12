# Tests for .build_conditional_projection() and .conditional_mt_scan()
# (Stage 6B-2), covering Anchors 1-6.

bcp <- CondPED:::.build_conditional_projection
cms <- CondPED:::.conditional_mt_scan

make_null <- function(n = 12L, m = 2L, seed = 1L) {
  set.seed(seed)
  Y <- matrix(stats::rnorm(n * m), n, m,
              dimnames = list(NULL, paste0("Trait", seq_len(m))))
  W <- cbind(1, stats::rnorm(n))
  colnames(W) <- c("Intercept", "W1")
  Zg <- matrix(stats::rnorm(n * 40), n, 40)
  K <- CondPED:::.make_grm(Zg)
  fit <- fit_mt_null(Y, W = W, K = K, n_starts = 2L,
                     control = list(maxit = 200L))
  list(fit = fit, K = K, W = W, Y = Y)
}

# package stacking order: column (a-1)*m + t (covariate a, trait t),
# row (t-1)*n + j
make_design <- function(B, n, m) {
  qq <- ncol(B)
  X <- matrix(0, n * m, qq * m)
  for (a in seq_len(qq)) {
    for (t in seq_len(m)) {
      X[(t - 1L) * n + seq_len(n), (a - 1L) * m + t] <- B[, a]
    }
  }
  X
}

dense_P <- function(fit, K, W, X_C = NULL) {
  n <- nrow(K)
  m <- ncol(fit$rotation$Y_tilde)
  V <- kronecker(fit$Sigma_G, K) + kronecker(fit$Sigma_E, diag(n))
  Vinv <- solve(V)
  Xstar <- make_design(W, n, m)
  if (!is.null(X_C)) {
    Xstar <- cbind(Xstar, make_design(X_C, n, m))
  }
  S <- crossprod(Xstar, Vinv %*% Xstar)
  Sinv <- CondPED:::.safe_inverse(S)$inverse
  P <- Vinv - Vinv %*% Xstar %*% Sinv %*% crossprod(Xstar, Vinv)
  list(P = P, S = S, Xstar = Xstar, Vinv = Vinv)
}

# ---- Anchor 1: empty conditioning set reduces to the original scan ----------

test_that("Anchor 1: empty X_C reproduces P, U, J, Q of the original scan", {
  env <- make_null()
  fit <- env$fit
  set.seed(2)
  G <- matrix(stats::rbinom(12 * 6, 2, 0.3), ncol = 6)
  proj0 <- bcp(fit, NULL)
  rot <- fit$rotation
  expect_equal(proj0$S_XX, rot$XtVinvX, tolerance = 1e-8)
  expect_equal(proj0$G_inv, rot$XtVinvX_inv, tolerance = 1e-8)
  expect_equal(proj0$beta_star, rot$beta_fixed, tolerance = 1e-8)
  expect_equal(proj0$residual_tilde, rot$residual_tilde, tolerance = 1e-8)
  expect_identical(proj0$A_arr, rot$Vinv)

  scan <- scan_mt_omnibus(fit, G, return_effects = TRUE)
  cscan <- cms(proj0, G)
  expect_equal(cscan$conditional$Q, scan$omnibus$Q, tolerance = 1e-10)
  expect_equal(cscan$conditional$p_value, scan$omnibus$p_value,
               tolerance = 1e-10)
  expect_equal(unname(cscan$effects$beta),
               unname(scan$effects$beta), tolerance = 1e-10)
})

# ---- Anchor 2: one conditioning SNP equals the explicit dense residual maker --

test_that("Anchor 2: block projection matches explicit dense P_C", {
  env <- make_null(n = 8L)
  fit <- env$fit
  set.seed(3)
  x_c <- stats::rbinom(8, 2, 0.3)
  proj <- bcp(fit, x_c)
  ref <- dense_P(fit, env$K, env$W, matrix(x_c, ncol = 1))
  expect_equal(proj$S_XX, ref$S, tolerance = 1e-8)
  yvec <- as.vector(env$Y)
  beta_dense <- CondPED:::.safe_inverse(ref$S)$inverse %*%
    crossprod(ref$Xstar, ref$Vinv %*% yvec)
  expect_equal(proj$beta_star, drop(beta_dense), tolerance = 1e-8)
  # raw residuals y - X* beta
  resid_dense <- yvec - ref$Xstar %*% beta_dense
  resid_block <- as.vector(proj$U %*% proj$residual_tilde)
  expect_equal(drop(resid_block), drop(resid_dense), tolerance = 1e-8)
  # the score input P_C y = V^-1 (y - X* beta), in rotated blocks
  Py <- matrix(ref$P %*% yvec, nrow = 8L)
  Py_rot <- crossprod(proj$U, Py)   # n x m, row j = coordinate block
  expect_equal(proj$Ar[, 1], Py_rot[1, ], tolerance = 1e-8)
  expect_equal(proj$Ar, t(Py_rot), tolerance = 1e-8)
})

# ---- Anchor 3: conditional score equals direct GLS ---------------------------

test_that("Anchor 3: conditional U/J/Q match direct dense GLS", {
  env <- make_null(n = 8L)
  fit <- env$fit
  set.seed(4)
  x_c <- matrix(stats::rbinom(8, 2, 0.3), ncol = 1)
  x_j <- stats::rbinom(8, 2, 0.3)
  proj <- bcp(fit, x_c)
  ref <- dense_P(fit, env$K, env$W, x_c)
  m <- 2L
  Z_j <- make_design(matrix(x_j, ncol = 1), 8L, m)
  U_dense <- crossprod(Z_j, ref$P %*% as.vector(env$Y))
  J_dense <- crossprod(Z_j, ref$P %*% Z_j)
  x_tilde <- crossprod(fit$rotation$U, x_j)
  block <- CondPED:::.gls_block_components(
    x_tilde, proj$AM_arr, proj$Ar, proj$A_arr, proj$G_inv,
    sqrt(.Machine$double.eps)
  )
  expect_equal(block$U, drop(U_dense), tolerance = 1e-8)
  expect_equal(block$J, J_dense, tolerance = 1e-8)
  Q_dense <- sum(U_dense * (CondPED:::.safe_inverse(J_dense)$inverse %*%
                              U_dense))
  cscan <- cms(proj, matrix(x_j, ncol = 1))
  expect_equal(cscan$conditional$Q[1], Q_dense, tolerance = 1e-8)
})

# ---- Anchor 4: conditioning removes LD-induced association ---------------------

test_that("Anchor 4: conditioning on the causal SNP removes the LD shadow", {
  set.seed(5)
  n <- 400
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- x1 + stats::rbinom(n, 1, 0.02)          # high LD with x1
  x2 <- pmin(x2, 2)
  Y <- matrix(x1 * 0.8 + stats::rnorm(n), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  K <- diag(n)
  fit <- fit_mt_null(Y, K = K, n_starts = 2L,
                     control = list(maxit = 200L))
  G <- cbind(x1, x2)
  colnames(G) <- c("causal", "shadow")
  marg <- scan_mt_omnibus(fit, G)
  p_marginal <- marg$omnibus$p_value[marg$omnibus$marker_id == "shadow"]
  expect_lt(p_marginal, 1e-4)          # x2 marginally associated via LD

  proj <- bcp(fit, x1)
  cond <- cms(proj, G)
  p_cond <- cond$conditional$p_value[
    cond$conditional$marker_id == "shadow"]
  Q_marg <- marg$omnibus$Q[marg$omnibus$marker_id == "shadow"]
  Q_cond <- cond$conditional$Q[cond$conditional$marker_id == "shadow"]
  expect_gt(p_cond, 0.05)
  expect_lt(Q_cond, Q_marg * 0.2)

  # cross-check against plain OLS partialling-out (m = 1, K = I):
  # the conditional score test tracks the t^2 of x2 in lm(Y ~ x1 + x2)
  # (same order of magnitude; the null-model V vs residual sigma^2
  # make them close but not identical)
  tt <- summary(stats::lm(Y ~ x1 + x2))$coefficients["x2", "t value"]^2
  expect_equal(Q_cond / tt, 1, tolerance = 0.5)
})

# ---- Anchor 5: V stays fixed ---------------------------------------------------

test_that("Anchor 5: Sigma_G, Sigma_E and the V blocks never change", {
  env <- make_null()
  fit <- env$fit
  proj0 <- bcp(fit, NULL)
  proj1 <- bcp(fit, stats::rbinom(12, 2, 0.3))
  proj2 <- bcp(fit, cbind(stats::rbinom(12, 2, 0.3),
                          stats::rbinom(12, 2, 0.4)))
  expect_identical(proj0$A_arr, fit$rotation$Vinv)
  expect_identical(proj1$A_arr, fit$rotation$Vinv)
  expect_identical(proj2$A_arr, fit$rotation$Vinv)
  expect_identical(proj1$Sigma_G, fit$Sigma_G)
  expect_identical(proj2$Sigma_G, fit$Sigma_G)
  expect_identical(proj1$Sigma_E, fit$Sigma_E)
  expect_identical(proj2$Sigma_E, fit$Sigma_E)
  # only the projection moves
  expect_false(isTRUE(all.equal(proj0$S_XX, proj1$S_XX)))
})

# ---- Anchor 6: collinear conditioning input --------------------------------------

test_that("Anchor 6: rank-deficient conditioning uses pinv + status, no crash", {
  env <- make_null()
  fit <- env$fit
  x1 <- stats::rbinom(12, 2, 0.3)
  X_C <- cbind(x1, x1)                      # exact duplicate columns
  proj <- bcp(fit, X_C)
  expect_identical(proj$status$code, "rank_deficient")
  expect_true(proj$used_pseudoinverse)
  expect_lt(proj$rank_S_XX, ncol(proj$S_XX))
  set.seed(6)
  G <- matrix(stats::rbinom(12 * 4, 2, 0.3), ncol = 4)
  cscan <- cms(proj, G)
  expect_identical(nrow(cscan$conditional), 4L)
  expect_true(all(cscan$conditional$status %in%
                    c("ok", "rank_deficient", "numerical_error")))
})

test_that("input validation", {
  env <- make_null()
  expect_error(bcp(env$fit, matrix(1, 5, 1)),
               class = "condped_invalid_input")
  expect_error(bcp(list()), class = "condped_invalid_input")
  proj <- bcp(env$fit, NULL)
  expect_error(cms(list(), matrix(1, 12, 2)),
               class = "condped_invalid_input")
})
