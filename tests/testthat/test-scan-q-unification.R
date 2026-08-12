# Regression: scan_mt_omnibus() must produce byte-identical Q/df/p/
# rank/status after its Q logic was unified into the shared
# `.q_from_block()` path (one formula, one implementation).
# The expected constants below were captured from the pre-refactor
# inline implementation on this fixed fixture (including a zero-variance
# all-ones dosage marker and a partially missing marker).

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
                 1.785045852548, 0.000000000000, 2.506233655448),
               tolerance = 1e-10)
  expect_identical(d$df, c(2L, 2L, 2L, 2L, 2L, 2L))
  expect_equal(d$p_value,
               c(0.143804422339, 0.333681584701, 0.901156620390,
                 0.409621004410, 1.000000000000, 0.285613200963),
               tolerance = 1e-10)
  expect_identical(d$rank_J, c(2L, 2L, 2L, 2L, 2L, 2L))
  expect_identical(d$status, rep("ok", 6L))
})
