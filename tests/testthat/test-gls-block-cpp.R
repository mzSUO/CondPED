# Golden-standard regression tests for the Rcpp GLS block kernel.
# The R reference `.gls_block_components_r()` and the C++ kernel behind
# `.gls_block_components()` / `.gls_blocks_chunk()` must agree to 1e-10 on
# U/J/J_inv/Q/df/p/beta and exactly on rank/condition/status, including the
# rank-deficient zero floor and all scan-level filter edge cases.

make_scan_fixture <- function(seed = 20260826) {
  set.seed(seed)
  n <- 60; m <- 3
  Y <- matrix(rnorm(n * m), n, m, dimnames = list(NULL, paste0("T", 1:m)))
  W <- cbind(1, rnorm(n)); colnames(W) <- c("Intercept", "W1")
  K <- CondPED:::.make_grm(matrix(rnorm(n * 30), n, 30))
  fit <- fit_mt_null(Y, W = W, K = K, control = list(maxit = 200))
  G <- matrix(rbinom(n * 8, 2, 0.3), n, 8)
  colnames(G) <- paste0("snp", 1:8)
  G[2, 3] <- NA                     # partially missing (mean-imputed in chunk)
  G[, 5] <- 1                       # intercept-collinear: zero floor -> rank 0
  G[, 6] <- 2                       # monomorphic (maf = 0 -> filtered)
  G[1, 7] <- 0; G[2:n, 7] <- 2      # low MAF (filtered by maf_min)
  G[, 8] <- NA_integer_             # all missing (filtered)
  list(fit = fit, G = G, n = n, m = m)
}

test_that("C++ kernel matches the R reference per marker (1e-10)", {
  fx <- make_scan_fixture()
  fit <- fx$fit; m <- fx$m
  rot <- fit$rotation
  A_arr <- rot$Vinv
  AM_arr <- CondPED:::.precompute_AM(rot)
  Ar <- CondPED:::.precompute_Ar(rot)
  G_inv <- rot$XtVinvX_inv
  tol <- sqrt(.Machine$double.eps)

  G_imp <- fx$G
  cm <- colMeans(fx$G, na.rm = TRUE)
  for (s in seq_len(ncol(G_imp))) {
    if (anyNA(G_imp[, s]) && is.finite(cm[s])) {
      G_imp[is.na(G_imp[, s]), s] <- cm[s]
    }
  }
  Xt <- crossprod(rot$U, G_imp)

  chunk <- CondPED:::.gls_blocks_chunk(Xt, AM_arr, Ar, A_arr, G_inv, tol)
  expect_identical(length(chunk), ncol(Xt))
  ## only markers 1:5 pass the scan filter and therefore read their block
  ## in production (snp6/7/8 are filtered before block use; the all-NA
  ## marker yields NaN blocks in C++ and a hard error in the R reference,
  ## which is unreachable in the scan flow)
  for (s in 1:5) {
    br <- CondPED:::.gls_block_components_r(Xt[, s], AM_arr, Ar, A_arr,
                                            G_inv, tol)
    bc1 <- CondPED:::.gls_block_components(Xt[, s], AM_arr, Ar, A_arr,
                                           G_inv, tol)
    bcc <- chunk[[s]]
    for (bc in list(bc1, bcc)) {
      expect_equal(bc$U, br$U, tolerance = 1e-10)
      expect_equal(bc$J, br$J, tolerance = 1e-10)
      expect_equal(bc$J_inv, br$J_inv, tolerance = 1e-10)
      expect_equal(bc$beta, br$beta, tolerance = 1e-10)
      expect_identical(bc$rank, br$rank)
      expect_identical(bc$status, br$status)
      ## condition is a derived ratio of LAPACK eigenvalues; last-bit
      ## differences between R eigen() and arma eig_sym are expected
      expect_equal(bc$condition, br$condition, tolerance = 1e-10)
      ## shared Q convention must agree too
      qr_ <- CondPED:::.q_from_block(br, m)
      qc_ <- CondPED:::.q_from_block(bc, m)
      expect_equal(qc_$Q, qr_$Q, tolerance = 1e-10)
      expect_equal(qc_$p_value, qr_$p_value, tolerance = 1e-10)
      expect_identical(qc_$status, qr_$status)
    }
  }
  ## snp5 is intercept-collinear: the zero floor must give rank 0
  expect_identical(chunk[[5]]$rank, 0L)
  expect_identical(chunk[[5]]$status, "rank_deficient")
  expect_identical(chunk[[5]]$condition, Inf)
})

test_that("scan_mt_omnibus on the C++ path matches R-reference blocks", {
  fx <- make_scan_fixture()
  d <- scan_mt_omnibus(fx$fit, fx$G)$omnibus
  expect_identical(d$status,
                   c("ok", "ok", "ok", "ok", "rank_deficient", "filtered",
                     "filtered", "filtered"))
  expect_identical(d$filter_reason,
                   c("", "", "", "", "", "monomorphic", "low_maf",
                     "all_missing"))

  ## recompute tested markers with the R reference and compare Q/df/p
  rot <- fx$fit$rotation
  A_arr <- rot$Vinv
  AM_arr <- CondPED:::.precompute_AM(rot)
  Ar <- CondPED:::.precompute_Ar(rot)
  G_inv <- rot$XtVinvX_inv
  tol <- sqrt(.Machine$double.eps)
  m <- fx$m
  for (s in 1:5) {
    x <- fx$G[, s]
    if (anyNA(x)) x[is.na(x)] <- mean(x, na.rm = TRUE)
    xt <- drop(crossprod(rot$U, x))
    br <- CondPED:::.gls_block_components_r(xt, AM_arr, Ar, A_arr, G_inv, tol)
    qs <- CondPED:::.q_from_block(br, m)
    expect_equal(d$Q[s], qs$Q, tolerance = 1e-10)
    expect_equal(d$p_value[s], qs$p_value, tolerance = 1e-10)
    expect_identical(d$df[s], as.integer(br$rank))
    expect_identical(d$rank_J[s], as.integer(br$rank))
    expect_identical(d$status[s], qs$status)
  }
  ## filtered markers keep the frozen filtered schema
  for (s in 6:8) {
    expect_true(is.na(d$Q[s]))
    expect_identical(d$df[s], 0L)
    expect_identical(d$rank_J[s], 0L)
  }
})

test_that("chunk kernel handles return_score / return_effects paths", {
  fx <- make_scan_fixture()
  d0 <- scan_mt_omnibus(fx$fit, fx$G)
  d1 <- scan_mt_omnibus(fx$fit, fx$G, return_score = TRUE,
                        return_effects = TRUE)
  expect_identical(d0$omnibus, d1$omnibus)
  expect_identical(length(d1$score), ncol(fx$G))
  ## effects of the rank-deficient marker are NA, others finite
  eff <- d1$effects$effects_long
  expect_true(all(is.na(eff$beta[eff$marker_id == "snp5"])))
  expect_true(all(is.finite(eff$beta[eff$marker_id == "snp1"])))
  expect_identical(d1$diagnostics$n_rank_deficient, 1L)
  expect_identical(d1$diagnostics$n_filtered, 3L)
})

test_that("C++ kernel preserves m x m matrix shapes when m = 1", {
  set.seed(11)
  n <- 60
  Z <- matrix(rnorm(n * 200), n, 200)
  K <- tcrossprod(Z) / 200; K <- K / mean(diag(K))
  W <- cbind(1, rnorm(n))
  Y <- matrix(W %*% c(0.5, -0.8) +
                drop(t(chol(K)) %*% rnorm(n)) * sqrt(0.6) +
                rnorm(n, sd = sqrt(0.4)), ncol = 1)
  fit <- fit_mt_null(Y, W = W[, 2, drop = FALSE], K = K, n_starts = 2L,
                     control = list(maxit = 200))
  rot <- fit$rotation
  A_arr <- rot$Vinv
  AM_arr <- CondPED:::.precompute_AM(rot)
  Ar <- CondPED:::.precompute_Ar(rot)
  G_inv <- rot$XtVinvX_inv
  tol <- sqrt(.Machine$double.eps)
  set.seed(12)
  x <- rbinom(n, 2, 0.3)
  xt <- drop(crossprod(rot$U, matrix(x, ncol = 1)))
  br <- CondPED:::.gls_block_components_r(xt, AM_arr, Ar, A_arr, G_inv, tol)
  bc <- CondPED:::.gls_block_components(xt, AM_arr, Ar, A_arr, G_inv, tol)
  expect_true(is.matrix(bc$J) && all(dim(bc$J) == c(1L, 1L)))
  expect_true(is.matrix(bc$J_inv) && all(dim(bc$J_inv) == c(1L, 1L)))
  expect_equal(bc$J, br$J, tolerance = 1e-10)
  expect_equal(bc$J_inv, br$J_inv, tolerance = 1e-10)
  expect_identical(bc$rank, br$rank)
  expect_identical(bc$status, br$status)
  ## end-to-end m=1 scan path (regression: dim-dropping broke effects)
  scan <- scan_mt_omnibus(fit, matrix(x, ncol = 1), marker_ids = "M1",
                          return_effects = TRUE)
  expect_identical(scan$omnibus$status, "ok")
  expect_true(is.finite(scan$effects$effects_long$se))
})
