# Regression: the scale-aware zero floor in `.gls_block_components()`
# (Stage 7.3.5) must declare rank 0 only for structurally/numerically
# zero information (marker collinear with the fixed-effect design),
# never for weak-but-estimable polymorphic markers.

test_that("zero floor separates structural zero from weak polymorphic markers", {
  set.seed(20260815)
  n <- 24; m <- 2
  Y <- matrix(rnorm(n * m), n, m,
              dimnames = list(NULL, paste0("Trait", 1:m)))
  W <- cbind(1, rnorm(n)); colnames(W) <- c("Intercept", "W1")
  K <- CondPED:::.make_grm(matrix(rnorm(n * 40), n, 40))
  fit <- fit_mt_null(Y, W = W, K = K, n_starts = 2L,
                     control = list(maxit = 200))

  z <- rnorm(n)
  G <- cbind(
    ones = rep(1, n),                    # intercept-collinear: structural zero
    numzero = 1 + 1e-12 * z,             # residual information ~1e-24: numerical zero
    weak = 1 + 0.05 * z,                 # weak but estimable polymorphic-ish dosage
    ordinary = rbinom(n, 2, 0.3)         # ordinary polymorphic marker
  )
  d <- scan_mt_omnibus(fit, G)$omnibus

  # Structural / numerical zero information: rank 0, frozen rank-deficient
  # policy (Q/p = NA, df = 0).
  expect_identical(d$rank_J[d$marker_id == "ones"], 0L)
  expect_identical(d$status[d$marker_id == "ones"], "rank_deficient")
  expect_identical(d$df[d$marker_id == "ones"], 0L)
  expect_true(is.na(d$Q[d$marker_id == "ones"]))
  expect_true(is.na(d$p_value[d$marker_id == "ones"]))

  expect_identical(d$rank_J[d$marker_id == "numzero"], 0L)
  expect_identical(d$status[d$marker_id == "numzero"], "rank_deficient")

  # Weak-but-estimable and ordinary markers: full rank m, status ok.
  expect_identical(d$rank_J[d$marker_id == "weak"], 2L)
  expect_identical(d$status[d$marker_id == "weak"], "ok")
  expect_true(is.finite(d$Q[d$marker_id == "weak"]))

  expect_identical(d$rank_J[d$marker_id == "ordinary"], 2L)
  expect_identical(d$status[d$marker_id == "ordinary"], "ok")
})
