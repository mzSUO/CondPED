# Fast deterministic unit tests for simulate_condped_data().
# Heavy statistical calibration (empirical variances, eigenvalue dispersion
# gates) lives in inst/validation/validate-simulate-condped-data.R.

test_that("return structure and dimensions match the contract", {
  sim <- simulate_condped_data(n = 120, m = 4, p = 300, n_qtl = 2, seed = 11)
  expect_equal(dim(sim$Y), c(120L, 4L))
  expect_equal(dim(sim$W), c(120L, 1L))
  expect_equal(dim(sim$G), c(120L, 300L))
  expect_equal(dim(sim$K_bg), c(120L, 120L))
  expect_length(sim$qtl_index, 2L)
  expect_true(all(sim$qtl_index %in% seq_len(300)))

  truth <- sim$truth
  expect_equal(dim(truth$B_Q), c(2L, 4L))
  expect_equal(dim(truth$beta), c(2L, 4L))
  expect_identical(truth$B_Q, truth$beta)
  expect_length(truth$A, 2L)
  expect_length(truth$D, 2L)
  expect_equal(dim(truth$eta), c(2L, 4L))
  expect_equal(dim(truth$locus_pve), c(2L, 4L))
  expect_equal(dim(truth$conditional_pve), c(2L, 4L))
  for (nm in c("Sigma_Q", "Sigma_G_bg", "Sigma_G_total", "Sigma_E",
               "Sigma_P_total")) {
    expect_equal(dim(truth[[nm]]), c(4L, 4L), label = nm)
  }
  expect_equal(dim(truth$gamma), c(4L, 3L))
  expect_equal(dim(truth$C), c(4L, 4L))

  expect_equal(dim(sim$latent$U), c(120L, 4L))
  expect_equal(dim(sim$latent$E), c(120L, 4L))
  expect_null(simulate_condped_data(
    n = 60, m = 2, p = 100, return_latent = FALSE, seed = 1
  )$latent)

  expect_true(sim$status$ok)
  expect_identical(sim$status$code, "ok")
  expect_true(all(c("mean_diag_K", "K_eigen_sd", "K_rank",
                    "min_eigen_Sigma_G_bg") %in% names(sim$diagnostics)))
  expect_true(all(c("seed", "attempts", "settings") %in% names(sim$generator)))
})

test_that("K_bg is scaled to mean diagonal 1", {
  sim <- simulate_condped_data(n = 150, m = 4, p = 400, seed = 21)
  expect_lt(abs(sim$diagnostics$mean_diag_K - 1), 1e-10)
  expect_lt(abs(mean(diag(sim$K_bg)) - 1), 1e-10)
})

test_that("true contrasts satisfy the model-scale orthogonality identity", {
  sim <- simulate_condped_data(n = 100, m = 4, p = 200, seed = 31)
  Sigma_P <- sim$truth$Sigma_P_total
  gamma <- sim$truth$gamma
  C <- sim$truth$C
  for (i in seq_len(4)) {
    # ||Sigma_{-i,-i} gamma_{i,-i} - Sigma_{-i,i}||_inf < 1e-10
    resid <- Sigma_P[-i, -i] %*% gamma[i, ] - Sigma_P[-i, i]
    expect_lt(max(abs(resid)), 1e-10)
    # contrast vector: c_i = e_i - E_{-i} gamma_{i,-i}
    expect_identical(unname(C[i, i]), 1)
    expect_equal(unname(C[-i, i]), -gamma[i, ], tolerance = 1e-12)
  }
})

test_that("eta identity holds for every architecture", {
  for (arch in c("null", "single_trait", "shared_same", "shared_opposite",
                 "dense", "covariance_aligned", "conditional_deviation",
                 "projection_induced")) {
    sim <- simulate_condped_data(
      n = 80, m = 4, p = 200, architecture = arch,
      n_qtl = 2, delta = 0.5, seed = 41
    )
    expect_equal(
      sim$truth$eta, sim$truth$beta %*% sim$truth$C,
      tolerance = 1e-12, label = paste("eta = beta %*% C for", arch)
    )
    # D is always a subset of A
    for (l in seq_len(2)) {
      expect_true(all(sim$truth$D[[l]] %in% sim$truth$A[[l]]),
                  label = paste("D subset of A for", arch))
    }
  }
})

test_that("covariance_aligned forces the target-trait eta to zero", {
  sim <- simulate_condped_data(
    n = 100, m = 4, p = 200, architecture = "covariance_aligned",
    target_trait = 1L, n_qtl = 2, seed = 51
  )
  expect_lt(max(abs(sim$truth$eta[, 1])), 1e-10)
  # the effect itself is generally non-zero on the target trait
  expect_true(any(sim$truth$beta[, 1] != 0))
})

test_that("conditional_deviation produces non-zero target-trait deviation", {
  sim <- simulate_condped_data(
    n = 100, m = 4, p = 200, architecture = "conditional_deviation",
    target_trait = 1L, delta = 0.5, n_qtl = 2, seed = 55
  )
  expect_true(all(sim$truth$eta[, 1] > 0))
  # and the deviation places the target trait in D
  for (l in seq_len(2)) {
    expect_true(1L %in% sim$truth$D[[l]])
  }
})

test_that("projection_induced has beta = 0 but eta != 0 on the target trait", {
  sim <- simulate_condped_data(
    n = 100, m = 4, p = 200, architecture = "projection_induced",
    target_trait = 1L, n_qtl = 2, seed = 61
  )
  expect_identical(as.numeric(sim$truth$beta[, 1]), c(0, 0))
  expect_true(all(abs(sim$truth$eta[, 1]) > 1e-3))
  # projection-induced conditional signal must not enter A
  for (l in seq_len(2)) {
    expect_false(1L %in% sim$truth$A[[l]])
  }
})

test_that("covariance ground truth is consistent and PSD", {
  sim <- simulate_condped_data(n = 100, m = 4, p = 200,
                               architecture = "dense", n_qtl = 2, seed = 71)
  truth <- sim$truth
  expect_equal(truth$Sigma_G_bg + truth$Sigma_Q, truth$Sigma_G_total,
               tolerance = 1e-12)
  expect_equal(truth$Sigma_G_total + truth$Sigma_E, truth$Sigma_P_total,
               tolerance = 1e-12)
  # total phenotypic variance is standardised to 1 per trait
  expect_equal(unname(diag(truth$Sigma_P_total)), rep(1, 4),
               tolerance = 1e-12)
  expect_gte(min(eigen(truth$Sigma_G_bg, symmetric = TRUE,
                       only.values = TRUE)$values), -1e-10)
  expect_gte(sim$diagnostics$min_eigen_Sigma_G_bg, -1e-10)
})

test_that("PSD enforcement rescales and keeps every truth component in sync", {
  # Deliberately impossible target: 3 dense loci with average PVE 0.4 each
  # cannot fit inside a total genetic variance of 0.5.
  sim <- simulate_condped_data(
    n = 100, m = 4, p = 200, architecture = "dense", n_qtl = 3,
    locus_pve = 0.4, max_attempts = 2L, seed = 81
  )
  expect_true(sim$generator$scaled)
  expect_length(sim$status$warnings, 1L)
  truth <- sim$truth
  expect_gte(sim$diagnostics$min_eigen_Sigma_G_bg, -1e-8)
  expect_equal(truth$Sigma_G_bg + truth$Sigma_Q, truth$Sigma_G_total,
               tolerance = 1e-8)
  expect_equal(truth$eta, truth$beta %*% truth$C, tolerance = 1e-10)
  for (l in seq_len(3)) {
    expect_identical(truth$A[[l]], which(truth$beta[l, ] != 0))
    expect_true(all(truth$D[[l]] %in% truth$A[[l]]))
  }

  # the projection identity survives scaling as well
  sim_ca <- simulate_condped_data(
    n = 100, m = 4, p = 200, architecture = "covariance_aligned",
    n_qtl = 3, locus_pve = 0.4, max_attempts = 2L, seed = 81
  )
  expect_true(sim_ca$generator$scaled)
  expect_lt(max(abs(sim_ca$truth$eta[, 1])), 1e-10)
})

test_that("the same seed reproduces the data set exactly", {
  a <- simulate_condped_data(n = 100, m = 4, p = 300, n_qtl = 2, seed = 91)
  b <- simulate_condped_data(n = 100, m = 4, p = 300, n_qtl = 2, seed = 91)
  expect_identical(a, b)
  c <- simulate_condped_data(n = 100, m = 4, p = 300, n_qtl = 2, seed = 92)
  expect_false(identical(a$Y, c$Y))
})

test_that("group structure yields eigenvalue dispersion in K_bg", {
  sim_s <- simulate_condped_data(n = 200, m = 2, p = 2000, structured = TRUE,
                                 n_groups = 8L, fst = 0.05, seed = 101)
  sim_u <- simulate_condped_data(n = 200, m = 2, p = 2000, structured = FALSE,
                                 seed = 101)
  # deterministic in-environment: mild Balding-Nichols differentiation
  # (fst = 0.05, 8 groups) lifts sd(eigen K) by a stable factor of ~1.8
  expect_gt(sim_s$diagnostics$K_eigen_sd,
            1.5 * sim_u$diagnostics$K_eigen_sd)
})

test_that("null architecture carries no focal signal", {
  sim <- simulate_condped_data(n = 80, m = 4, p = 200,
                               architecture = "null", n_qtl = 2, seed = 111)
  expect_identical(as.numeric(sim$truth$beta), rep(0, 8))
  expect_identical(as.numeric(sim$truth$Sigma_Q), rep(0, 16))
  for (l in seq_len(2)) {
    expect_length(sim$truth$A[[l]], 0L)
    expect_length(sim$truth$D[[l]], 0L)
  }
  expect_true(sim$status$ok)
})

test_that("input validation rejects invalid arguments", {
  expect_error(simulate_condped_data(n = 1), "n must be")
  expect_error(simulate_condped_data(m = 1), "m must be")
  expect_error(simulate_condped_data(maf_range = c(0.6, 0.7)), "maf_range")
  expect_error(simulate_condped_data(target_trait = 99L), "target_trait")
  expect_error(simulate_condped_data(h2 = 1.5), "h2")
  expect_error(simulate_condped_data(p = 10, n_qtl = 20L), "n_qtl")
  expect_error(simulate_condped_data(structured = TRUE, fst = 0), "fst")
  R_bad <- diag(4); R_bad[1, 2] <- R_bad[2, 1] <- 0.5
  expect_error(
    simulate_condped_data(R_G = R_bad * 2), "R_G"
  )
})
