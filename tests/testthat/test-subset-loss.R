# Tests for .compute_subset_loss() (v1.0 Stage 1).
#
# Covers the Runbook 3.5 hand-calculation anchors A-C and the frozen
# mathematical invariants (interface contract sections 3.4-3.6, 5.3,
# 9.1). No subset enumeration, no representative-set or module
# extraction in this stage.

csl <- CondPED:::.compute_subset_loss

# ---- Anchor A: independent traits -----------------------------------------

test_that("anchor A: independent traits give the hand-computed loss", {
  Sigma <- diag(3)
  dimnames(Sigma) <- list(c("A", "B", "C"), c("A", "B", "C"))
  beta <- c(A = 1, B = 2, C = 0)
  out <- csl(beta, Sigma, "A")
  expect_equal(out$full_qform, 5)
  expect_equal(out$subset_qform, 1)
  expect_equal(out$residual_qform, 4)
  expect_equal(out$representation_loss, 0.8)
  expect_equal(out$decomposition_error, 0)
  expect_identical(out$status$code, "ok")
  expect_true(out$status$ok)
  expect_identical(out$representing_set, "A")
  expect_identical(out$complement_set, c("B", "C"))
  expect_identical(out$representing_key, "A")
  expect_identical(out$complement_key, "B|C")
})

# ---- Anchor B: exact representation through correlation -------------------

test_that("anchor B: correlated traits allow exact representation", {
  r <- 0.6
  Sigma <- matrix(c(1, r, r, 1), 2, 2,
                  dimnames = list(c("A", "B"), c("A", "B")))
  beta <- c(A = 1, B = r)
  out <- csl(beta, Sigma, "A")
  expect_equal(drop(out$Gamma), r)
  expect_equal(unname(out$eta), 0, tolerance = 1e-12)
  expect_equal(out$residual_qform, 0, tolerance = 1e-12)
  expect_equal(out$representation_loss, 0, tolerance = 1e-12)
})

# ---- Anchor C: scale invariance -------------------------------------------

test_that("anchor C: loss is invariant to beta scaling and trait rescaling", {
  Sigma <- matrix(c(1, 0.4, 0.1,
                    0.4, 2, 0.3,
                    0.1, 0.3, 1.5), 3, 3,
                  dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  beta <- c(A = 0.8, B = -0.5, C = 1.2)
  ref <- csl(beta, Sigma, "A")$representation_loss

  # beta scaled by 3
  expect_equal(csl(3 * beta, Sigma, "A")$representation_loss, ref,
               tolerance = 1e-10)

  # trait units rescaled by a nonzero diagonal D: Y -> D Y implies
  # Sigma -> D Sigma D and beta -> D beta
  d <- c(2, 0.5, 3)
  Sigma2 <- diag(d) %*% Sigma %*% diag(d)
  beta2 <- beta * d
  expect_equal(csl(beta2, Sigma2, "A")$representation_loss, ref,
               tolerance = 1e-10)
})

# ---- boundary sets ----------------------------------------------------------

test_that("empty set: loss 1, subset 0, residual full, no Gamma", {
  Sigma <- diag(c(1, 2))
  dimnames(Sigma) <- list(c("A", "B"), c("A", "B"))
  beta <- c(A = 1, B = 1)
  for (empty in list(NULL, character(), integer())) {
    out <- csl(beta, Sigma, empty)
    expect_identical(out$representation_loss, 1)
    expect_identical(out$subset_qform, 0)
    expect_identical(out$residual_qform, out$full_qform)
    expect_null(out$Gamma)
    expect_identical(out$representing_key, "<empty>")
    expect_identical(out$complement_key, "A|B")
    expect_identical(out$rank_Sigma_SS, 0L)
  }
})

test_that("full set: loss 0, residual 0, subset full", {
  Sigma <- diag(c(1, 2))
  dimnames(Sigma) <- list(c("A", "B"), c("A", "B"))
  beta <- c(A = 1, B = 1)
  out <- csl(beta, Sigma, c("B", "A"))   # input order must not matter
  expect_identical(out$representation_loss, 0)
  expect_identical(out$residual_qform, 0)
  expect_identical(out$subset_qform, out$full_qform)
  expect_identical(out$representing_key, "A|B")
  expect_identical(out$complement_key, "<empty>")
  expect_identical(out$rank_Omega, 0L)
})

# ---- mathematical identities ------------------------------------------------

test_that("normal equations, Schur complement and decomposition hold", {
  Sigma <- matrix(c(2.0, 0.5, 0.2,
                    0.5, 1.5, 0.4,
                    0.2, 0.4, 1.0), 3, 3,
                  dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  beta <- c(A = 0.7, B = -1.1, C = 0.4)
  S <- c("A", "B"); R <- "C"
  out <- csl(beta, Sigma, S)

  # normal equations: Sigma_SS %*% Gamma == Sigma_SR (base solve as
  # independent reference is allowed in tests)
  expect_equal(unname(out$Gamma),
               unname(solve(Sigma[S, S], Sigma[S, R, drop = FALSE])),
               tolerance = 1e-10)
  # eta identity
  expect_equal(unname(out$eta),
               unname(beta[R] - drop(crossprod(
                 solve(Sigma[S, S], Sigma[S, R, drop = FALSE]), beta[S]))),
               tolerance = 1e-10)
  # Schur complement
  expect_equal(unname(out$Omega),
               unname(Sigma[R, R, drop = FALSE] -
                        Sigma[R, S, drop = FALSE] %*%
                        solve(Sigma[S, S], Sigma[S, R, drop = FALSE])),
               tolerance = 1e-10)
  # Mahalanobis decomposition: full = subset + residual
  expect_equal(out$full_qform, out$subset_qform + out$residual_qform,
               tolerance = 1e-10)
  expect_lt(out$decomposition_error, 1e-10)
  # loss consistent with the ratio
  expect_equal(out$representation_loss,
               out$residual_qform / out$full_qform, tolerance = 1e-12)
})

test_that("loss lies in [0, 1] and is monotone on nested pairs", {
  Sigma <- matrix(c(2.0, 0.9, 0.3,
                    0.9, 1.5, 0.5,
                    0.3, 0.5, 1.0), 3, 3,
                  dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  beta <- c(A = 1.2, B = 0.4, C = -0.7)
  rho <- function(S) csl(beta, Sigma, S)$representation_loss
  for (S in list(NULL, "A", "B", "C", c("A", "B"), c("A", "C"),
                 c("B", "C"), c("A", "B", "C"))) {
    l <- rho(S)
    expect_gte(l, 0)
    expect_lte(l, 1)
  }
  # nested pairs: S subset of T implies rho(T) <= rho(S)
  expect_lte(rho(c("A", "B")), rho("A") + 1e-10)
  expect_lte(rho(c("A", "B")), rho("B") + 1e-10)
  expect_lte(rho("A"), rho(NULL) + 1e-10)
  expect_lte(rho(c("A", "B", "C")), rho(c("B", "C")) + 1e-10)
})

test_that("trait permutation does not change the result", {
  Sigma <- matrix(c(2.0, 0.5, 0.2,
                    0.5, 1.5, 0.4,
                    0.2, 0.4, 1.0), 3, 3,
                  dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  beta <- c(A = 0.7, B = -1.1, C = 0.4)
  perm <- c(3, 1, 2)
  out1 <- csl(beta, Sigma, "B")
  out2 <- csl(beta[perm], Sigma[perm, perm], "B")
  expect_equal(out2$representation_loss, out1$representation_loss,
               tolerance = 1e-12)
  expect_equal(out2$full_qform, out1$full_qform, tolerance = 1e-12)
  # set keys follow each call's own trait_names order, never the order
  # traits were listed in representing_set
  expect_identical(out1$complement_key, "A|C")
  expect_identical(out2$complement_key, "C|A")  # permuted reference order
  out3 <- csl(beta, Sigma, c("B"))             # same call, key stable
  expect_identical(csl(beta, Sigma, "B")$complement_key,
                   csl(beta, Sigma, c("B"))$complement_key)
})

# ---- numerical states --------------------------------------------------------

test_that("near-zero full quadratic form is not_applicable", {
  Sigma <- diag(2)
  dimnames(Sigma) <- list(c("A", "B"), c("A", "B"))
  out <- csl(c(A = 0, B = 0), Sigma, "A")
  expect_identical(out$status$code, "not_applicable")
  expect_true(out$status$ok)
  expect_true(is.na(out$representation_loss))
  tiny <- csl(c(A = 1e-8, B = 0), Sigma, "A")
  expect_identical(tiny$status$code, "not_applicable")  # full = 1e-16 <= 1e-12
})

test_that("rank-deficient matrices: diagnostics kept, loss not usable", {
  # Sigma_SS singular: trait A has zero variance
  Sigma <- diag(c(0, 1, 1))
  dimnames(Sigma) <- list(c("A", "B", "C"), c("A", "B", "C"))
  beta <- c(A = 1, B = 2, C = 3)
  out <- csl(beta, Sigma, "A")
  expect_identical(out$status$code, "rank_deficient")
  expect_false(out$status$ok)                 # cannot enter optimisation
  expect_true(is.na(out$representation_loss)) # no formal loss offered
  expect_true(out$used_pseudoinverse)
  expect_identical(out$rank_Sigma_SS, 0L)
  # pinv diagnostics remain available: Gamma = 0, eta = beta_R,
  # residual = 13, raw loss = 1
  expect_equal(out$residual_qform, 13)
  expect_equal(out$raw_representation_loss, 1)

  # the same rule applies to the boundary branches
  out_empty <- csl(beta, Sigma, character())
  expect_identical(out_empty$status$code, "rank_deficient")
  expect_false(out_empty$status$ok)
  expect_true(is.na(out_empty$representation_loss))
  expect_equal(out_empty$raw_representation_loss, 1)
  out_full <- csl(beta, Sigma, c("A", "B", "C"))
  expect_identical(out_full$status$code, "rank_deficient")
  expect_false(out_full$status$ok)
  expect_true(is.na(out_full$representation_loss))
  expect_equal(out_full$raw_representation_loss, 0)
})

test_that("gross out-of-range losses are marked unstable, never truncated", {
  # Indefinite Omega drives the residual quadratic form negative beyond
  # any floating-point tolerance; the raw value must be kept, not
  # clipped, and no formal loss may be offered.
  Sigma <- diag(c(2, -1, 3))
  dimnames(Sigma) <- list(c("A", "B", "C"), c("A", "B", "C"))
  beta <- c(A = 2, B = 0.2, C = 0.3)
  out <- csl(beta, Sigma, "A")
  # hand: full = 2 - 0.04 + 0.03 = 1.99; subset = 2; residual = -0.01
  expect_equal(out$full_qform, 1.99, tolerance = 1e-12)
  expect_equal(out$residual_qform, -0.01, tolerance = 1e-12)
  expect_equal(out$raw_representation_loss, -0.01 / 1.99, tolerance = 1e-12)
  expect_identical(out$status$code, "unstable")
  expect_false(out$status$ok)
  expect_true(is.na(out$representation_loss))
})

test_that("only floating-point-scale excursions are truncated", {
  # A mathematically exact boundary case computed with rounding error:
  # rho may come out as 1 + 1e-16; it is truncated to 1 with status ok.
  r <- 0.6
  Sigma <- matrix(c(1, r, r, 1), 2, 2,
                  dimnames = list(c("A", "B"), c("A", "B")))
  beta <- c(A = 1, B = r)
  out <- csl(beta, Sigma, character())   # exact 1 by construction
  expect_identical(out$representation_loss, 1)
  expect_identical(out$status$code, "ok")
  expect_true(out$status$ok)
})

test_that("high condition numbers produce warnings, not failures", {
  Sigma <- diag(c(1, 1e-11))
  dimnames(Sigma) <- list(c("A", "B"), c("A", "B"))
  out <- csl(c(A = 1, B = 1), Sigma, "A", inverse_tol = 1e-13)
  expect_identical(out$status$code, "ok")
  expect_match(paste(out$status$warnings, collapse = " "),
               "condition number")
})

# ---- set normalisation and input validation -----------------------------------

test_that("set normalisation: dedup, order, integer positions, keys", {
  ns <- CondPED:::.normalize_trait_set
  key <- CondPED:::.trait_set_key
  expect_identical(ns(c("B", "A", "B"), c("A", "B", "C")), c("A", "B"))
  expect_identical(ns(c(3, 1), c("A", "B", "C")), c("A", "C"))
  expect_identical(ns(NULL, c("A", "B")), character())
  expect_identical(key(character()), "<empty>")
  expect_identical(key(c("A", "C")), "A|C")
  expect_error(ns("Z", c("A", "B")), class = "condped_invalid_input")
  expect_error(ns(9, c("A", "B")), class = "condped_invalid_input")
})

test_that("invalid inputs raise condped_invalid_input", {
  Sigma <- diag(2)
  dimnames(Sigma) <- list(c("A", "B"), c("A", "B"))
  expect_error(csl(c(1, 2), Sigma, "A"), class = "condped_invalid_input")
  expect_error(csl(c(A = 1, B = NA), Sigma, "A"),
               class = "condped_invalid_input")
  expect_error(csl(c(A = 1, B = 2), Sigma[, 1, drop = FALSE], "A"),
               class = "condped_invalid_input")
  expect_error(csl(c(A = 1, B = 2), matrix(c(1, 0, 0.5, 1), 2), "A"),
               class = "condped_invalid_input")  # asymmetric, unnamed
  expect_error(csl(c(A = 1, B = 2), Sigma, "Z"),
               class = "condped_invalid_input")
  expect_error(csl(c(A = 1, B = 2), Sigma, "A", inverse_tol = -1),
               class = "condped_invalid_input")
  expect_error(csl(c(A = 1, B = 2), Sigma, "A", trait_names = "A"),
               class = "condped_invalid_input")
})

test_that("Sigma dimnames are reordered to match beta", {
  Sigma <- matrix(c(1, 0.5, 0.5, 2), 2, 2,
                  dimnames = list(c("B", "A"), c("B", "A")))
  beta <- c(A = 1, B = 0.5)
  out <- csl(beta, Sigma, "A")
  expect_identical(out$status$code, "ok")
  # independent reference with explicit ordering
  ref <- csl(beta, Sigma[c("A", "B"), c("A", "B")], "A")
  expect_equal(out$representation_loss, ref$representation_loss,
               tolerance = 1e-12)
})
