# Tests for fit_mt_null() (S2): rotated REML correctness, log-Cholesky
# parameterisation, m = 1 degenerate case, weak-identifiability warning.

# ---- helpers ---------------------------------------------------------------

# Dense REML/ML reference: builds the full nm x nm matrices (small n only).
dense_null_loglik <- function(theta, m, K, Y, W, reml = TRUE,
                              chol_floor = 1e-8) {
  n <- nrow(Y)
  q <- ncol(W)
  ncp <- m * (m + 1) / 2
  SG <- CondPED:::.logchol_unpack(theta[seq_len(ncp)], m, chol_floor)
  SE <- CondPED:::.logchol_unpack(theta[ncp + seq_len(ncp)], m, chol_floor)
  V <- kronecker(SG, K) + kronecker(SE, diag(n))
  X0 <- kronecker(diag(m), W)
  y <- as.vector(Y)
  Vc <- chol(V)
  Vinv <- chol2inv(Vc)
  logdetV <- 2 * sum(log(diag(Vc)))
  XtVinvX <- crossprod(X0, Vinv %*% X0)
  XtVinvy <- crossprod(X0, Vinv %*% y)
  beta <- solve(XtVinvX, XtVinvy)
  r <- y - X0 %*% beta
  quad <- sum(r * (Vinv %*% r))
  if (reml) {
    logdetX <- as.numeric(determinant(XtVinvX, logarithm = TRUE)$modulus)
    -0.5 * ((n * m - q * m) * log(2 * pi) + logdetV + logdetX + quad)
  } else {
    -0.5 * (n * m * log(2 * pi) + logdetV + quad)
  }
}

make_null_problem <- function(n = 50, m = 2, q = 2, seed = 1) {
  set.seed(seed)
  Z <- matrix(rnorm(n * 300), n, 300)
  K <- tcrossprod(Z) / 300
  K <- K / mean(diag(K))
  SG <- matrix(c(1, 0.3, 0.3, 0.8), m, m)[seq_len(m), seq_len(m), drop = FALSE]
  SE <- matrix(c(0.9, -0.2, -0.2, 1.1), 2, 2)[seq_len(m), seq_len(m), drop = FALSE]
  W <- cbind(1, rnorm(n))[, seq_len(q), drop = FALSE]
  B <- matrix(rnorm(q * m), q, m)
  u <- t(chol(K)) %*% matrix(rnorm(n * m), n, m) %*% chol(SG)
  Y <- W %*% B + u + matrix(rnorm(n * m), n, m) %*% chol(SE)
  list(K = K, Y = Y, W = W, SG = SG, SE = SE, m = m, q = q)
}

rotated_inputs <- function(prob) {
  eig <- eigen(prob$K, symmetric = TRUE)
  list(
    lambda = eig$values,
    Y_tilde = crossprod(eig$vectors, prob$Y),
    W_tilde = crossprod(eig$vectors, prob$W)
  )
}

# ---- rotated vs dense objective --------------------------------------------

test_that("rotated REML matches the dense reference to 1e-7", {
  prob <- make_null_problem(n = 12, seed = 11)
  rot <- rotated_inputs(prob)
  thetas <- list(
    c(CondPED:::.logchol_pack(prob$SG), CondPED:::.logchol_pack(prob$SE)),
    rnorm(6, sd = 0.5),
    rnorm(6, sd = 2)
  )
  for (theta in thetas) {
    for (reml in c(TRUE, FALSE)) {
      val_rot <- CondPED:::.reml_objective_rotated(
        theta, prob$m, rot$lambda, rot$Y_tilde, rot$W_tilde, reml = reml)
      val_dense <- dense_null_loglik(theta, prob$m, prob$K, prob$Y, prob$W,
                                     reml = reml)
      expect_lt(abs(val_rot - val_dense), 1e-7)
    }
  }
})

test_that("analytic gradient matches finite differences", {
  prob <- make_null_problem(n = 12, seed = 12)
  rot <- rotated_inputs(prob)
  theta <- rnorm(6, sd = 0.4)
  vg <- CondPED:::.reml_value_grad(theta, prob$m, rot$lambda, rot$Y_tilde,
                                   rot$W_tilde)
  fd <- vapply(seq_along(theta), function(p) {
    h <- 1e-6 * (1 + abs(theta[p]))
    th2 <- theta; th2[p] <- th2[p] + h
    (CondPED:::.reml_objective_rotated(th2, prob$m, rot$lambda, rot$Y_tilde,
                                       rot$W_tilde) - vg$value) / h
  }, numeric(1))
  expect_lt(max(abs(vg$gradient - fd) / (1 + abs(fd))), 1e-4)
})

test_that("ML and REML objectives differ by the fixed-effect correction", {
  prob <- make_null_problem(n = 12, seed = 13)
  rot <- rotated_inputs(prob)
  theta <- c(CondPED:::.logchol_pack(prob$SG), CondPED:::.logchol_pack(prob$SE))
  ll_r <- CondPED:::.reml_objective_rotated(theta, prob$m, rot$lambda,
                                            rot$Y_tilde, rot$W_tilde,
                                            reml = TRUE)
  ll_m <- CondPED:::.reml_objective_rotated(theta, prob$m, rot$lambda,
                                            rot$Y_tilde, rot$W_tilde,
                                            reml = FALSE)
  expect_false(isTRUE(all.equal(ll_r, ll_m)))
  # REML - ML = -1/2 (log|X'V^{-1}X| - q m log 2pi)
  dense_r <- dense_null_loglik(theta, prob$m, prob$K, prob$Y, prob$W, TRUE)
  dense_m <- dense_null_loglik(theta, prob$m, prob$K, prob$Y, prob$W, FALSE)
  expect_equal(ll_r - ll_m, dense_r - dense_m, tolerance = 1e-7)
})

# ---- m = 1 degenerate case ---------------------------------------------------

test_that("m = 1 degenerates to the univariate linear mixed model", {
  set.seed(21)
  n <- 40
  Z <- matrix(rnorm(n * 300), n, 300)
  K <- tcrossprod(Z) / 300
  K <- K / mean(diag(K))
  sg2 <- 1.2
  se2 <- 0.7
  W <- cbind(1, rnorm(n))
  Y <- matrix(W %*% c(0.5, -1) +
                drop(t(chol(K)) %*% rnorm(n)) * sqrt(sg2) +
                rnorm(n, sd = sqrt(se2)), ncol = 1)
  eig <- eigen(K, symmetric = TRUE)
  rot <- list(lambda = eig$values,
              Y_tilde = crossprod(eig$vectors, Y),
              W_tilde = crossprod(eig$vectors, W))
  theta <- c(log(sqrt(sg2)), log(sqrt(se2)))
  val_rot <- CondPED:::.reml_objective_rotated(theta, 1L, rot$lambda,
                                               rot$Y_tilde, rot$W_tilde)
  val_dense <- dense_null_loglik(theta, 1L, K, Y, W, reml = TRUE)
  expect_equal(val_rot, val_dense, tolerance = 1e-7)

  fit <- fit_mt_null(Y, W = W[, -1, drop = FALSE], K = K, n_starts = 2L,
                     control = list(maxit = 300))
  expect_s3_class(fit, "condped_mt_null")
  expect_equal(dim(fit$gamma), c(1L, 0L))
  expect_equal(fit$Sigma_P_ref, fit$Sigma_G + fit$Sigma_E)
})

# ---- full fit ----------------------------------------------------------------

test_that("fit_mt_null returns positive-definite covariances and valid output", {
  sim <- simulate_condped_data(n = 150, m = 3, p = 300,
                               architecture = "null", seed = 31)
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 2L,
                     control = list(maxit = 300))
  expect_s3_class(fit, "condped_mt_null")
  expect_true(all(diag(fit$Sigma_G) > 0))
  expect_true(all(diag(fit$Sigma_E) > 0))
  for (S in list(fit$Sigma_G, fit$Sigma_E, fit$Sigma_P_ref)) {
    expect_true(all(eigen(S, symmetric = TRUE, only.values = TRUE)$values > 0))
  }
  expect_equal(fit$Sigma_P_ref, fit$Sigma_G + fit$Sigma_E)
  expect_equal(dim(fit$fixed_effects), c(1L, 3L))
  expect_equal(dim(fit$gamma), c(3L, 2L))
  expect_equal(dim(fit$contrasts), c(3L, 3L))
  expect_length(fit$conditional_variance, 3L)
  expect_true(is.finite(fit$logLik))
  expect_named(fit$convergence,
               c("code", "message", "gradient_norm", "iterations",
                 "selected_start", "starts"))
  expect_true(all(c("U", "lambda", "Y_tilde", "W_tilde", "Vinv", "logdet_Vj",
                    "XtVinvX", "XtVinvX_inv", "beta_fixed", "residual_tilde")
                  %in% names(fit$rotation)))
  expect_false(fit$diagnostics$identifiable_warning)
})

test_that("the optimum beats the starting point in log-likelihood", {
  sim <- simulate_condped_data(n = 150, m = 2, p = 300,
                               architecture = "null", seed = 32)
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 2L,
                     control = list(maxit = 300))
  theta_start <- c(CondPED:::.logchol_pack(diag(0.5 * apply(sim$Y, 2, stats::var))),
                   CondPED:::.logchol_pack(diag(0.5 * apply(sim$Y, 2, stats::var))))
  rot <- list(lambda = fit$rotation$lambda,
              Y_tilde = fit$rotation$Y_tilde,
              W_tilde = fit$rotation$W_tilde)
  ll_start <- CondPED:::.reml_objective_rotated(theta_start, 2L, rot$lambda,
                                                rot$Y_tilde, rot$W_tilde)
  expect_gte(fit$logLik, ll_start - 1e-6)
})

test_that("rotation matrices reproduce the model identity", {
  sim <- simulate_condped_data(n = 100, m = 2, p = 200,
                               architecture = "null", seed = 33)
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 1L,
                     control = list(maxit = 100))
  U <- fit$rotation$U
  lambda <- fit$rotation$lambda
  # U is orthogonal and reconstructs K
  expect_equal(crossprod(U), diag(nrow(U)), tolerance = 1e-8)
  expect_equal(unname(U %*% (lambda * t(U))), unname(sim$K_bg),
               tolerance = 1e-6)
  # rotated data match U'Y and U'W
  expect_equal(fit$rotation$Y_tilde, crossprod(U, sim$Y),
               tolerance = 1e-10)
  # Vinv blocks invert lambda_j Sigma_G + Sigma_E
  j <- 5L
  V_j <- lambda[j] * fit$Sigma_G + fit$Sigma_E
  expect_equal(unname(fit$rotation$Vinv[, , j] %*% V_j), diag(2),
               tolerance = 1e-8)
})

# ---- weak identifiability -----------------------------------------------------

test_that("K near identity triggers identifiable_warning", {
  set.seed(41)
  n <- 50
  Y <- matrix(rnorm(n * 2), n, 2)
  fit_i <- fit_mt_null(Y, K = diag(n), n_starts = 1L,
                       control = list(maxit = 100))
  expect_true(fit_i$diagnostics$identifiable_warning)
  expect_true(any(grepl("weakly identifiable", fit_i$status$warnings)))

  sim <- simulate_condped_data(n = 100, m = 2, p = 200,
                               architecture = "null", seed = 42)
  fit_k <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 1L,
                       control = list(maxit = 100))
  expect_false(fit_k$diagnostics$identifiable_warning)
})

# ---- reproducibility and options ----------------------------------------------

test_that("set.seed gives reproducible fits and return_rotation = FALSE works", {
  sim <- simulate_condped_data(n = 100, m = 2, p = 200,
                               architecture = "null", seed = 51)
  set.seed(1)
  fit1 <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 3L,
                      control = list(maxit = 150))
  set.seed(1)
  fit2 <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 3L,
                      control = list(maxit = 150))
  expect_identical(fit1$logLik, fit2$logLik)
  expect_identical(fit1$Sigma_G, fit2$Sigma_G)

  fit3 <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 1L,
                      control = list(maxit = 100), return_rotation = FALSE)
  expect_null(fit3$rotation)

  fit4 <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 1L,
                      control = list(maxit = 100), compute_gamma = FALSE)
  expect_null(fit4$gamma)
  expect_null(fit4$contrasts)
})

test_that("log-Cholesky pack/unpack round-trips", {
  S <- matrix(c(2, 0.4, 0.4, 1), 2, 2)
  theta <- CondPED:::.logchol_pack(S)
  expect_length(theta, 3L)
  expect_equal(CondPED:::.logchol_unpack(theta, 2L), S, tolerance = 1e-10)
  # floored diagonal keeps the matrix positive definite
  S_pd <- CondPED:::.logchol_unpack(c(-50, -50, 0), 2L, chol_floor = 1e-4)
  expect_gt(min(eigen(S_pd, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("gamma satisfies the orthogonality identity on the fitted Sigma_P", {
  sim <- simulate_condped_data(n = 150, m = 3, p = 300,
                               architecture = "null", seed = 61)
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 2L,
                     control = list(maxit = 300))
  SP <- fit$Sigma_P_ref
  for (i in seq_len(3L)) {
    lhs <- SP[-i, -i] %*% fit$gamma[i, ]
    expect_equal(drop(lhs), SP[-i, i], tolerance = 1e-8)
  }
  # contrasts and conditional variance are consistent with gamma
  for (i in seq_len(3L)) {
    c_i <- fit$contrasts[, i]
    expect_equal(unname(c_i[i]), 1)
    expect_equal(unname(c_i[-i]), -fit$gamma[i, ], tolerance = 1e-10)
    expect_equal(fit$conditional_variance[i],
                 SP[i, i] - sum(fit$gamma[i, ] * SP[-i, i]),
                 tolerance = 1e-8)
  }
})

test_that("residual_tilde is consistent with the GLS fixed effects", {
  sim <- simulate_condped_data(n = 100, m = 2, p = 200,
                               architecture = "null", seed = 62)
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 2L,
                     control = list(maxit = 200))
  rot <- fit$rotation
  # residual_tilde = Y_tilde - W_tilde %*% fixed_effects
  expect_equal(rot$residual_tilde,
               rot$Y_tilde - rot$W_tilde %*% fit$fixed_effects,
               tolerance = 1e-10)
  # beta_fixed is the vectorised version of fixed_effects
  expect_equal(matrix(rot$beta_fixed, nrow = 2L), unname(t(fit$fixed_effects)))
  # XtVinvX_inv inverts XtVinvX
  expect_equal(unname(rot$XtVinvX_inv %*% rot$XtVinvX),
               diag(nrow(rot$XtVinvX)), tolerance = 1e-8)
  # logdet_Vj matches a direct determinant computation at one coordinate
  j <- 3L
  V_j <- rot$lambda[j] * fit$Sigma_G + fit$Sigma_E
  expect_equal(rot$logdet_Vj[j],
               as.numeric(determinant(V_j, logarithm = TRUE)$modulus),
               tolerance = 1e-8)
})

test_that("every start is recorded and the selected one is the best legal", {
  sim <- simulate_condped_data(n = 100, m = 2, p = 200,
                               architecture = "null", seed = 63)
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 4L,
                     control = list(maxit = 200))
  st <- fit$convergence$starts
  expect_equal(nrow(st), 4L)
  expect_named(st, c("start", "logLik", "code", "converged", "legal",
                     "selected"))
  expect_equal(sum(st$selected), 1L)
  sel <- which(st$selected)
  expect_equal(sel, fit$convergence$selected_start)
  legal_ll <- st$logLik[st$legal]
  expect_gte(fit$logLik, max(legal_ll) - 1e-8)
  # distinct starting points actually differ
  expect_gt(length(unique(round(st$logLik, 8))), 1L)
})

# ---- invalid input handling ---------------------------------------------------

test_that("invalid inputs raise condped_invalid_input errors", {
  sim <- simulate_condped_data(n = 60, m = 2, p = 120,
                               architecture = "null", seed = 71)
  Y <- sim$Y
  K <- sim$K_bg

  # Y / K row mismatch
  expect_error(fit_mt_null(Y[-1, ], K = K), class = "condped_invalid_input")
  # non-square K
  expect_error(fit_mt_null(Y, K = K[-1, ]), class = "condped_invalid_input")
  # non-symmetric K (beyond numerical noise)
  K_asym <- K
  K_asym[1, 2] <- K_asym[1, 2] + 0.5
  expect_error(fit_mt_null(Y, K = K_asym), class = "condped_invalid_input")
  # NA / Inf in Y, W, K
  Y_na <- Y; Y_na[1, 1] <- NA
  expect_error(fit_mt_null(Y_na, K = K), class = "condped_invalid_input")
  Y_inf <- Y; Y_inf[1, 1] <- Inf
  expect_error(fit_mt_null(Y_inf, K = K), class = "condped_invalid_input")
  K_inf <- K; K_inf[2, 2] <- Inf
  expect_error(fit_mt_null(Y, K = K_inf), class = "condped_invalid_input")
  W_inf <- matrix(Inf, nrow(Y), 1)
  expect_error(fit_mt_null(Y, W = W_inf, K = K),
               class = "condped_invalid_input")
  # W row mismatch
  expect_error(fit_mt_null(Y, W = matrix(1, nrow(Y) - 1L, 1), K = K),
               class = "condped_invalid_input")
  # n_starts not a positive integer
  expect_error(fit_mt_null(Y, K = K, n_starts = 0),
               class = "condped_invalid_input")
  expect_error(fit_mt_null(Y, K = K, n_starts = -2),
               class = "condped_invalid_input")
  expect_error(fit_mt_null(Y, K = K, n_starts = 2.5),
               class = "condped_invalid_input")
  # h2_start outside (0, 1)
  expect_error(fit_mt_null(Y, K = K, h2_start = 0),
               class = "condped_invalid_input")
  expect_error(fit_mt_null(Y, K = K, h2_start = 1),
               class = "condped_invalid_input")
  expect_error(fit_mt_null(Y, K = K, h2_start = NA),
               class = "condped_invalid_input")
  # non-PSD K
  K_npd <- K
  K_npd[1, ] <- 10
  K_npd[, 1] <- 10
  diag(K_npd)[1] <- 1
  expect_error(fit_mt_null(Y, K = K_npd), class = "condped_invalid_input")
})

test_that("v1.0 interface fields: Sigma_P_ref, trait_names, individual_ids", {
  sim <- simulate_condped_data(n = 60, m = 3, p = 100,
                               architecture = "null", seed = 77)
  fit <- fit_mt_null(sim$Y, K = sim$K_bg,
                     control = list(maxit = 100))
  expect_equal(fit$Sigma_P_ref, fit$Sigma_G + fit$Sigma_E)
  expect_null(fit[["Sigma_P"]])   # old field name retired (exact match)
  expect_identical(fit$trait_names, colnames(sim$Y))
  expect_identical(fit$individual_ids, rownames(sim$Y))
  expect_true(is.finite(fit$diagnostics$min_eigen_Sigma_P_ref))
  expect_true(is.finite(fit$diagnostics$condition_Sigma_P_ref))
})
