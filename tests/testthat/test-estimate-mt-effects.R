# Tests for estimate_mt_effects() (S3)

make_small_null <- function(n = 12, m = 2, seed = 1) {
  set.seed(seed)
  Z <- matrix(rnorm(n * 200), n, 200)
  K <- tcrossprod(Z) / 200
  K <- K / mean(diag(K))
  SG <- matrix(0.3, m, m) + diag(0.7, m)
  SE <- matrix(0.1, m, m) + diag(0.9, m)
  W <- cbind(1, rnorm(n))
  B <- matrix(rnorm(2 * m), 2, m)
  u <- t(chol(K)) %*% matrix(rnorm(n * m), n, m) %*% chol(SG)
  Y <- W %*% B + u + matrix(rnorm(n * m), n, m) %*% chol(SE)
  fit_mt_null(Y, W = W[, 2, drop = FALSE], K = K, n_starts = 2L,
              control = list(maxit = 200))
}

test_that("beta = J^+ U and covariance = J^+", {
  fit <- make_small_null(n = 25, m = 3, seed = 11)
  set.seed(12)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 4, 2, 0.3), ncol = 4)
  est <- estimate_mt_effects(fit, G, targets = 1:4)

  for (i in seq_len(4)) {
    x_tilde <- crossprod(fit$rotation$U, G[, i])
    AM_arr <- CondPED:::.precompute_AM(fit$rotation)
    Ar <- CondPED:::.precompute_Ar(fit$rotation)
    A_arr <- fit$rotation$Vinv
    block <- CondPED:::.gls_block_components(
      x_tilde, AM_arr, Ar, A_arr, fit$rotation$XtVinvX_inv,
      sqrt(.Machine$double.eps))
    expect_equal(unname(est$beta[i, ]), block$beta, tolerance = 1e-10)
    expect_equal(est$covariance[, , i], block$J_inv, tolerance = 1e-10)
  }
})

test_that("loci can be specified by character marker_id", {
  fit <- make_small_null(n = 20, m = 2, seed = 21)
  set.seed(22)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 5, 2, 0.3), ncol = 5)
  colnames(G) <- paste0("snp", 1:5)
  est_num <- estimate_mt_effects(fit, G, targets = c(1, 3, 5))
  est_chr <- estimate_mt_effects(fit, G, targets = c("snp1", "snp3", "snp5"))
  expect_equal(est_num$beta, est_chr$beta, tolerance = 1e-10)
})

test_that("return_covariance = FALSE suppresses covariance", {
  fit <- make_small_null(n = 20, m = 2, seed = 31)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 2, 2, 0.3), ncol = 2)
  est <- estimate_mt_effects(fit, G, targets = 1, return_covariance = FALSE)
  expect_null(est$covariance)
})

test_that("rank-deficient loci produce NA standard errors", {
  fit <- make_small_null(n = 18, m = 2, seed = 41)
  # All-zero genotype -> J = 0 (rank deficient)
  x <- rep(0, nrow(fit$rotation$Y_tilde))
  est <- estimate_mt_effects(fit, matrix(x, ncol = 1), targets = 1L)
  expect_true(all(is.na(est$beta)))
  expect_true(all(is.na(est$effects_long$se)))
  expect_equal(est$diagnostics$n_rank_deficient, 1L)
})

test_that("effect dimensions match contract", {
  fit <- make_small_null(n = 20, m = 4, seed = 51)
  set.seed(52)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 3, 2, 0.3), ncol = 3)
  est <- estimate_mt_effects(fit, G, targets = 1:3)
  expect_equal(dim(est$beta), c(3L, 4L))
  expect_equal(dim(est$covariance), c(4L, 4L, 3L))
  expect_equal(nrow(est$effects_long), 12L)
  expect_named(est$effects_long, c("marker_id", "trait", "beta", "se", "z",
                                   "p_value"))
})

test_that("v1.0 interface fields: se, genotype_variance, locus_table, trait_names", {
  fit <- make_small_null(n = 20, m = 4, seed = 51)
  set.seed(52)
  G <- matrix(stats::rbinom(nrow(fit$rotation$Y_tilde) * 3, 2, 0.3), ncol = 3)
  colnames(G) <- paste0("snp", 1:3)
  est <- estimate_mt_effects(fit, G, targets = 1:3)
  # se equals sqrt(diag(covariance)) and matches effects_long
  expect_equal(unname(est$se[2, ]), sqrt(diag(est$covariance[, , 2])),
               tolerance = 1e-12)
  expect_equal(as.vector(t(est$se)), est$effects_long$se)
  # dimnames lock the ID/order contract
  expect_identical(rownames(est$beta), paste0("snp", 1:3))
  expect_identical(colnames(est$beta), est$trait_names)
  expect_identical(est$trait_names, fit$trait_names)
  # genotype variance from the raw dosages
  expect_equal(unname(est$genotype_variance), unname(apply(G, 2, stats::var)))
  expect_identical(names(est$genotype_variance), paste0("snp", 1:3))
  # locus_table columns and values
  expect_identical(names(est$locus_table),
                   c("marker_id", "maf", "genotype_variance", "rank_J",
                     "condition_J", "status"))
  expect_identical(est$locus_table$marker_id, paste0("snp", 1:3))
  expect_equal(est$locus_table$genotype_variance,
               unname(est$genotype_variance))
  expect_true(all(est$locus_table$status == "ok"))
})
