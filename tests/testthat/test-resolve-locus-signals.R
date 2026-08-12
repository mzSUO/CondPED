# Tests for the stepwise signal resolution core (Stage 6B-3),
# Anchors 1-10.

bcp <- CondPED:::.build_conditional_projection
resolve <- CondPED:::.resolve_one_locus_signals

mk_fit1 <- function(Y) {
  fit_mt_null(Y, K = diag(nrow(Y)), n_starts = 2L,
              control = list(maxit = 300L))
}

# Build a resolve() call around a hand-made genotype set.
run_resolve <- function(G, Y, lead, members = colnames(G),
                        screening_rule = "within_locus_bonferroni",
                        alpha_signal = 0.05, max_signals = 5L,
                        pos = NULL, chr = NULL, ...) {
  fit <- mk_fit1(Y)
  n <- nrow(G)
  marker_ids <- colnames(G)
  if (is.null(pos)) pos <- seq_along(marker_ids) * 1000
  if (is.null(chr)) chr <- rep("chr1", length(marker_ids))
  names(pos) <- names(chr) <- marker_ids
  resolve(
    locus_id = "chr1:test", lead_snp = lead, members = members,
    null_fit = fit, G = G, marker_ids = marker_ids,
    chromosome = chr, position = pos,
    screening_rule = screening_rule, alpha_signal = alpha_signal,
    max_signals = max_signals, ...
  )
}

# ---- Anchor 1: one-signal locus ------------------------------------------------

test_that("Anchor 1: single causal signal with LD proxies resolves to one", {
  set.seed(11)
  n <- 400
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- pmin(x1 + stats::rbinom(n, 1, 0.05), 2)   # LD proxy
  x3 <- stats::rbinom(n, 2, 0.3)                  # null
  Y <- matrix(x1 + stats::rnorm(n), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, proxy = x2, null = x3)
  out <- run_resolve(G, Y, lead = "S1")
  expect_identical(out$selected, "S1")
  expect_identical(out$diagnostics$stop_reason, "no_passing_marker")
  expect_identical(out$diagnostics$n_signals, 1L)
  expect_identical(out$signals$conditioning_key[1], "<empty>")
  expect_true(is.na(out$signals$conditional_p_raw[1]))
})

# ---- Anchor 2: two resolvable signals -------------------------------------------

test_that("Anchor 2: two independent causal signals both resolve", {
  set.seed(12)
  n <- 500
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- stats::rbinom(n, 2, 0.3)
  x3 <- stats::rbinom(n, 2, 0.3)
  Y <- matrix(x1 + 0.9 * x2 + stats::rnorm(n, sd = 0.8), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, S2 = x2, null = x3)
  out <- run_resolve(G, Y, lead = "S1")
  expect_identical(out$selected, c("S1", "S2"))
  expect_identical(out$signals$conditioning_key[2], "S1")
  expect_identical(out$diagnostics$n_signals, 2L)
})

# ---- Anchor 3: LD shadow does not become a signal --------------------------------

test_that("Anchor 3: pure LD shadow is not selected after conditioning", {
  set.seed(13)
  n <- 500
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- pmin(x1 + stats::rbinom(n, 1, 0.01), 2)   # near-perfect shadow
  x3 <- stats::rbinom(n, 2, 0.3)
  Y <- matrix(x1 + stats::rnorm(n, sd = 0.8), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, shadow = x2, null = x3)
  out <- run_resolve(G, Y, lead = "S1")
  expect_identical(out$selected, "S1")
  expect_false("shadow" %in% out$selected)
})

# ---- Anchor 4: iterative conditioning updates ---------------------------------------

test_that("Anchor 4: conditioning set grows and P_C is rebuilt each step", {
  set.seed(14)
  n <- 600
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- stats::rbinom(n, 2, 0.3)
  x3 <- stats::rbinom(n, 2, 0.3)
  Y <- matrix(x1 + x2 + x3 + stats::rnorm(n, sd = 0.9), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, S2 = x2, S3 = x3)
  out <- run_resolve(G, Y, lead = "S1", max_signals = 5L)
  expect_identical(out$diagnostics$n_signals, 3L)
  expect_identical(out$signals$n_conditioning, c(0L, 1L, 2L))
  # second and third signals are data-driven (S2/S3 in some order)
  expect_identical(out$signals$conditioning_before_selection[[2]], "S1")
  second <- out$selected[2]
  expect_true(second %in% c("S2", "S3"))
  expect_identical(sort(out$signals$conditioning_before_selection[[3]]),
                   sort(c("S1", second)))
  # the projection really was rebuilt: step-2 and step-3 scan tables
  # (over different candidate sets, from different P_C) differ
  expect_false(isTRUE(all.equal(out$history[[1]], out$history[[2]])))
  expect_identical(nrow(out$history[[1]]), 2L)   # 3 markers - 1 selected
  expect_identical(nrow(out$history[[2]]), 1L)   # 3 markers - 2 selected
})

# ---- Anchor 5: local Bonferroni denominator --------------------------------------------

test_that("Anchor 5: M_remaining is the pre-scan candidate count", {
  set.seed(15)
  n <- 600
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- stats::rbinom(n, 2, 0.3)
  filler <- replicate(3, stats::rbinom(n, 2, 0.3))
  Y <- matrix(x1 + x2 + stats::rnorm(n, sd = 0.8), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, S2 = x2, f1 = filler[, 1], f2 = filler[, 2],
             f3 = filler[, 3])
  out <- run_resolve(G, Y, lead = "S1")
  sig <- out$signals
  # step 2: 5 members - 1 selected -> M = 4
  expect_identical(sig$M_remaining[2], 4L)
  expect_identical(nrow(out$history[[1]]), 4L)
  expect_equal(sig$conditional_p_adjusted[2],
               min(1, 4 * sig$conditional_p_raw[2]))
  # step 3: 5 - 2 -> M = 3
  expect_identical(sig$M_remaining[3], 3L)
  expect_equal(sig$conditional_p_adjusted[3],
               min(1, 3 * sig$conditional_p_raw[3]),
               tolerance = 1e-12)
})

# ---- Anchor 6: max_signals caps selection -------------------------------------------------

test_that("Anchor 6: max_signals = 2 caps at two signals", {
  set.seed(16)
  n <- 600
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- stats::rbinom(n, 2, 0.3)
  x3 <- stats::rbinom(n, 2, 0.3)
  Y <- matrix(x1 + x2 + x3 + stats::rnorm(n, sd = 0.9), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, S2 = x2, S3 = x3)
  out <- run_resolve(G, Y, lead = "S1", max_signals = 2L)
  expect_length(out$selected, 2L)
  expect_identical(out$diagnostics$stop_reason, "max_signals")
})

# ---- Anchor 7: deterministic tie-break -------------------------------------------------------

test_that("Anchor 7: exact p ties break by chromosome/position/marker_id", {
  set.seed(17)
  n <- 500
  x1 <- stats::rbinom(n, 2, 0.3)
  xd <- stats::rbinom(n, 2, 0.3)
  Y <- matrix(x1 + 0.6 * xd + stats::rnorm(n, sd = 0.8), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, TIEpos200 = xd, TIEpos100 = xd)
  pos <- c(S1 = 500, TIEpos200 = 200, TIEpos100 = 100)
  out <- run_resolve(G, Y, lead = "S1", pos = pos, max_signals = 2L)
  # identical conditional statistics -> tie -> smaller position wins
  expect_equal(out$history[[1]]$p_value[
    out$history[[1]]$marker_id == "TIEpos200"],
    out$history[[1]]$p_value[
      out$history[[1]]$marker_id == "TIEpos100"])
  expect_identical(out$selected[2], "TIEpos100")
  # input column order must not matter
  out2 <- run_resolve(G[, c(2, 3, 1)], Y, lead = "S1",
                      pos = pos[c(2, 3, 1)], max_signals = 2L)
  expect_identical(out2$selected, out$selected)
})

# ---- Anchor 8: marker-order invariance -------------------------------------------------------

test_that("Anchor 8: shuffling marker order keeps the signal sequence", {
  set.seed(18)
  n <- 600
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- stats::rbinom(n, 2, 0.3)
  x3 <- stats::rbinom(n, 2, 0.3)
  Y <- matrix(x1 + x2 + stats::rnorm(n, sd = 0.8), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, S2 = x2, N1 = x3)
  pos <- c(S1 = 100, S2 = 200, N1 = 300)
  base <- run_resolve(G, Y, lead = "S1", pos = pos)
  perm <- c(3, 1, 2)
  shuf <- run_resolve(G[, perm], Y, lead = "S1", pos = pos[perm])
  expect_identical(shuf$selected, base$selected)
  expect_equal(shuf$signals$conditional_p_raw,
               base$signals$conditional_p_raw, tolerance = 1e-9)
})

# ---- Anchor 9: unusable candidates cannot win -------------------------------------------------

test_that("Anchor 9: constant/collinear candidates never win, no crash", {
  set.seed(19)
  n <- 500
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- stats::rbinom(n, 2, 0.3)
  const <- rep(1, n)                 # zero variance
  Y <- matrix(x1 + 0.9 * x2 + stats::rnorm(n, sd = 0.8), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, dup = x1, const = const, S2 = x2)
  expect_no_error(
    out <- run_resolve(G, Y, lead = "S1")
  )
  expect_false("dup" %in% out$selected)
  expect_false("const" %in% out$selected)
  expect_identical(out$selected[1:2], c("S1", "S2"))
})

# ---- Anchor 10: very high LD between true signals --------------------------------------------

test_that("Anchor 10: near-unresolvable high-LD signals stay stable", {
  set.seed(20)
  n <- 600
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- pmin(x1 + stats::rbinom(n, 1, 0.01), 2)   # r2 ~ 0.99
  Y <- matrix(x1 + x2 + stats::rnorm(n), ncol = 1,
              dimnames = list(NULL, "Trait1"))
  G <- cbind(S1 = x1, S2 = x2)
  out <- run_resolve(G, Y, lead = "S1", max_signals = 3L)
  expect_true(out$diagnostics$n_signals <= 2L)
  expect_identical(out$selected[1], "S1")
  expect_identical(anyDuplicated(out$selected), 0L)
  expect_true(out$status$ok)
})

# ---- validation ---------------------------------------------------------------------------

test_that("input validation", {
  set.seed(21)
  n <- 100
  G <- cbind(a = stats::rbinom(n, 2, 0.3), b = stats::rbinom(n, 2, 0.3))
  Y <- matrix(stats::rnorm(n), ncol = 1, dimnames = list(NULL, "Trait1"))
  expect_error(run_resolve(G, Y, lead = "zzz"),
               class = "condped_invalid_input")
  expect_error(run_resolve(G, Y, lead = "a", max_signals = 0L),
               class = "condped_invalid_input")
  expect_error(run_resolve(G, Y, lead = "a", alpha_signal = 2),
               class = "condped_invalid_input")
})
