# Tests for Stage 6B-4: final joint signal-specific effect estimation.

fjs <- CondPED:::.fit_joint_signal_effects

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

# ---- A. one signal: joint == marginal GLS ---------------------------------------

test_that("A: one-signal joint beta equals the marginal GLS beta", {
  env <- make_null()
  set.seed(2)
  G <- matrix(stats::rbinom(12 * 3, 2, 0.3), ncol = 3,
              dimnames = list(NULL, paste0("s", 1:3)))
  joint <- fjs(env$fit, G, "s2")
  marg <- estimate_mt_effects(env$fit, G, targets = "s2")
  expect_equal(unname(joint$beta[1, ]),
               unname(marg$beta[1, ]), tolerance = 1e-10)
  expect_equal(unname(joint$covariance[, , 1]), marg$covariance[, , 1],
               tolerance = 1e-10)
})

# ---- B. two signals: joint GLS equals explicit dense GLS -------------------------

test_that("B: two-signal joint GLS matches explicit dense computation", {
  env <- make_null(n = 10L)
  fit <- env$fit
  n <- 10L; m <- 2L
  set.seed(3)
  X_S <- matrix(stats::rbinom(n * 2, 2, 0.3), n, 2,
                dimnames = list(NULL, c("a", "b")))
  joint <- fjs(fit, X_S, c("a", "b"))

  V <- kronecker(fit$Sigma_G, env$K) + kronecker(fit$Sigma_E, diag(n))
  Vinv <- solve(V)
  X0 <- make_design(env$W, n, m)
  S0 <- crossprod(X0, Vinv %*% X0)
  P <- Vinv - Vinv %*% X0 %*% CondPED:::.safe_inverse(S0)$inverse %*%
    crossprod(X0, Vinv)
  Z_S <- make_design(X_S, n, m)
  U_dense <- crossprod(Z_S, P %*% as.vector(env$Y))
  J_dense <- crossprod(Z_S, P %*% Z_S)
  Jinv_dense <- CondPED:::.safe_inverse(J_dense)$inverse
  beta_dense <- Jinv_dense %*% U_dense

  expect_equal(joint$J, J_dense, tolerance = 1e-8)
  expect_equal(joint$U, drop(U_dense), tolerance = 1e-8)
  expect_equal(joint$J_inv, Jinv_dense, tolerance = 1e-8)
  expect_equal(as.vector(t(joint$beta)), drop(beta_dense),
               tolerance = 1e-8)
})

# ---- C. conditioned estimate == full joint block ----------------------------------

test_that("C: target conditioned on the others equals the joint block", {
  env <- make_null(n = 14L)
  set.seed(4)
  G <- matrix(stats::rbinom(14 * 4, 2, 0.3), ncol = 4,
              dimnames = list(NULL, c("s1", "s2", "s3", "s4")))
  S <- c("s1", "s2", "s3")
  joint <- fjs(env$fit, G, S)
  for (i in seq_along(S)) {
    est <- estimate_mt_effects(
      env$fit, G, targets = S[i],
      conditioning_sets = list(setdiff(S, S[i]))
    )
    expect_equal(unname(est$beta[1, ]), unname(joint$beta[i, ]),
                 tolerance = 1e-10)
    expect_equal(est$covariance[, , 1], joint$covariance[, , i],
                 tolerance = 1e-10)
    expect_equal(unname(est$se[1, ]), unname(joint$se[i, ]),
                 tolerance = 1e-10)
  }
})

# ---- D. LD: joint beta closer to the truth than marginal --------------------------

test_that("D: joint effects beat marginals under LD", {
  set.seed(5)
  n <- 800
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- x1 + stats::rbinom(n, 1, 0.15)         # strong but separable LD
  x2 <- pmin(x2, 2)
  b_true <- c(1.0, 0.5)
  Y <- matrix(x1 * b_true[1] + x2 * b_true[2] +
                stats::rnorm(n, sd = 0.5), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  fit <- fit_mt_null(Y, K = diag(n), n_starts = 2L,
                     control = list(maxit = 300L))
  G <- cbind(s1 = x1, s2 = x2)
  marg <- estimate_mt_effects(fit, G, targets = c("s1", "s2"))
  joint <- fjs(fit, G, c("s1", "s2"))
  err_marg <- abs(unname(marg$beta[1, ]) - b_true)
  err_joint <- abs(unname(joint$beta[, 1]) - b_true)
  expect_true(all(err_joint < err_marg))
  expect_true(all(err_joint < 0.15))
})

# ---- E. input order invariance ----------------------------------------------------

test_that("E: signal order does not matter after aligning by id", {
  env <- make_null()
  set.seed(6)
  G <- matrix(stats::rbinom(12 * 3, 2, 0.3), ncol = 3,
              dimnames = list(NULL, c("s1", "s2", "s3")))
  j12 <- fjs(env$fit, G, c("s1", "s2"))
  j21 <- fjs(env$fit, G, c("s2", "s1"))
  expect_equal(j21$beta["s1", ], j12$beta["s1", ], tolerance = 1e-12)
  expect_equal(j21$beta["s2", ], j12$beta["s2", ], tolerance = 1e-12)
  expect_equal(j21$covariance[, , "s2"], j12$covariance[, , "s2"],
               tolerance = 1e-12)
  # joint covariance block (s1, s2) symmetric under reordering
  m <- 2L
  expect_equal(j21$joint_covariance[1:m, (m + 1):(2 * m)],
               t(j12$joint_covariance[1:m, (m + 1):(2 * m)]),
               tolerance = 1e-12)
})

# ---- F. rank-deficient joint design -------------------------------------------------

test_that("F: duplicate signals trigger pinv + status, no crash", {
  env <- make_null()
  set.seed(7)
  x <- stats::rbinom(12, 2, 0.3)
  G <- cbind(d1 = x, d2 = x, ok = stats::rbinom(12, 2, 0.3))
  joint <- fjs(env$fit, G, c("d1", "d2", "ok"))
  expect_identical(joint$status$code, "rank_deficient")
  expect_true(joint$used_pseudoinverse)
  expect_lt(joint$rank, 6L)
  expect_true(all(is.finite(joint$beta["ok", ])))
})

# ---- G. global covariance stays fixed ------------------------------------------------

test_that("G: the joint path never re-estimates the null covariance", {
  env <- make_null()
  fit <- env$fit
  G <- matrix(stats::rbinom(12 * 2, 2, 0.3), ncol = 2,
              dimnames = list(NULL, c("s1", "s2")))
  sg <- fit$Sigma_G; se_ <- fit$Sigma_E; vinv <- fit$rotation$Vinv
  fjs(fit, G, c("s1", "s2"))
  expect_identical(fit$Sigma_G, sg)
  expect_identical(fit$Sigma_E, se_)
  expect_identical(fit$rotation$Vinv, vinv)
})

# ---- public resolve_locus_signals() -----------------------------------------------------

test_that("resolve_locus_signals: selection then joint refit, clean outputs", {
  set.seed(8)
  n <- 500
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- stats::rbinom(n, 2, 0.3)
  x3 <- pmin(x1 + stats::rbinom(n, 1, 0.03), 2)   # LD proxy of x1
  Y <- matrix(x1 + 0.9 * x2 + stats::rnorm(n, sd = 0.8), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, S2 = x2, proxy = x3)
  pos <- c(S1 = 1000, S2 = 2000, proxy = 3000)
  chr <- rep("chr1", 3)
  names(chr) <- names(pos) <- colnames(G)
  fit <- fit_mt_null(Y, K = diag(n), n_starts = 2L,
                     control = list(maxit = 300L))
  scan <- scan_mt_omnibus(fit, G)
  lo <- define_associated_loci(
    scan, G, chr, pos,
    selected_markers = c("S1", "S2"),
    method = "physical", window_bp = 5000
  )
  out <- resolve_locus_signals(fit, G, lo)
  expect_identical(out$diagnostics$n_loci, 1L)
  expect_identical(sort(out$signals$representative_snp),
                   sort(c("S1", "S2")))
  # beta/se come from the joint model of {S1, S2}
  joint <- fjs(fit, G, out$signals$representative_snp)
  expect_equal(out$beta$beta[out$beta$representative_snp == "S1"],
               unname(joint$beta["S1", ]), tolerance = 1e-12)
  expect_equal(out$beta$se[out$beta$representative_snp == "S2"],
               unname(joint$se["S2", ]), tolerance = 1e-12)
  expect_equal(dim(out$covariance[[1]]), c(2L, 2L))
  expect_identical(out$status$code, "ok")
})

test_that("resolve_locus_signals validates inputs", {
  env <- make_null()
  G <- matrix(stats::rbinom(12 * 2, 2, 0.3), ncol = 2,
              dimnames = list(NULL, c("a", "b")))
  lo <- list(
    loci = data.frame(locus_id = "chr1:1-2", lead_snp = "a",
                      lead_p = 0.001, stringsAsFactors = FALSE),
    membership = data.frame(locus_id = "chr1:1-2",
                            marker_id = c("a", "b"),
                            position = c(1, 2),
                            stringsAsFactors = FALSE)
  )
  expect_error(resolve_locus_signals(env$fit, G, lo,
                                     signal_adjust = "fixed"),
               class = "condped_invalid_input")
  expect_error(resolve_locus_signals(env$fit, G, list()),
               class = "condped_invalid_input")
})
