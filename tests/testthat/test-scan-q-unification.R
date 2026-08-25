# Regression: scan_mt_omnibus() must produce stable Q/df/p/rank/status
# after its Q logic was unified into the shared `.q_from_block()` path
# (one formula, one implementation).
#
# Q/p constants for markers 1-4 and 6 were captured from the pre-refactor
# inline implementation on this fixed fixture; they are checked with a
# numerical tolerance because cross-BLAS/cross-platform runs differ at
# the ~1e-10 level (no bitwise equality across BLAS is required).
#
# Marker 5 is an all-ones dosage, exactly collinear with the intercept in
# W, so its residualized genotype and its information matrix J are
# theoretically zero. The earlier frozen values (rank_J = 2, Q = 0, p = 1,
# status = "ok") were an artifact of ranking a pure rounding-noise J
# (entries ~1e-14 against an information scale of ~86) relative to its
# own largest eigenvalue, which is BLAS-dependent. `.gls_block_components()`
# now applies a deterministic scale-aware zero floor, so marker 5 follows
# the frozen rank-deficient policy on every platform:
# rank_J = 0, df = 0, Q/p = NA, status = "rank_deficient".

test_that("scan output is unchanged after unifying the Q path", {
  set.seed(20260813)
  n <- 24; m <- 2
  Y <- matrix(rnorm(n * m), n, m,
              dimnames = list(NULL, paste0("Trait", 1:m)))
  W <- cbind(1, rnorm(n)); colnames(W) <- c("Intercept", "W1")
  K <- CondPED:::.make_grm(matrix(rnorm(n * 40), n, 40))
  fit <- fit_mt_null(Y, W = W, K = K, n_starts = 2L,
                     control = list(maxit = 200))
  G <- matrix(rbinom(n * 6, 2, 0.3), ncol = 6)
  colnames(G) <- paste0("snp", 1:6)
  G[3, 2] <- NA
  G[, 5] <- 1
  G[1:12, 6] <- 0; G[13:24, 6] <- 2
  d <- scan_mt_omnibus(fit, G)$omnibus

  expect_equal(d$Q,
               c(3.878602161532, 2.195136159882, 0.208152413939,
                 1.785045852548, NA_real_, 2.506233655448),
               tolerance = 1e-8)
  expect_identical(d$df, c(2L, 2L, 2L, 2L, 0L, 2L))
  expect_equal(d$p_value,
               c(0.143804422339, 0.333681584701, 0.901156620390,
                 0.409621004410, NA_real_, 0.285613200963),
               tolerance = 1e-8)
  expect_identical(d$rank_J, c(2L, 2L, 2L, 2L, 0L, 2L))
  expect_identical(d$status, c(rep("ok", 4L), "rank_deficient", "ok"))
})
