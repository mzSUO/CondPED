# Tests for Stage 6C: ASSET comparator adapter + comparison pipelines,
# Anchors 36-43. ASSET is NOT installed in this environment: the adapter
# must degrade to status "asset_not_available" (backend = NULL) and all
# ASSET-side behaviour is exercised through a fake recording backend.

asset_cmp <- CondPED:::.run_asset_comparison
pipelines <- CondPED:::.run_comparison_pipelines

# Fake ASSET backend: records every call and returns the frozen raw
# shape (p + subset) using a deterministic |z| >= 2 inclusion rule.
mk_recording_backend <- function() {
  store <- new.env()
  store$calls <- list()
  backend <- function(beta, se, z, Sigma_Z, trait_names, two_sided) {
    store$calls[[length(store$calls) + 1L]] <- list(
      beta = beta, se = se, z = z, Sigma_Z = Sigma_Z,
      trait_names = trait_names, two_sided = two_sided
    )
    keep <- trait_names[abs(z[trait_names]) >= 2]
    if (length(keep) == 0L) keep <- trait_names[which.max(abs(z))]
    list(p = stats::pchisq(max(z^2), df = 1, lower.tail = FALSE),
         subset = keep)
  }
  list(backend = backend, store = store)
}

# Correlated Gaussian noise with correlation rho (no MASS dependency).
mk_corr_noise <- function(n, rho, m = 2L) {
  Sig <- matrix(rho, m, m)
  diag(Sig) <- 1
  matrix(stats::rnorm(n * m), n, m) %*% chol(Sig)
}

mk_locus_object <- function(locus_id, lead_snp, lead_p, members, pos,
                            r2_to_lead) {
  list(
    loci = data.frame(
      locus_id = locus_id, chromosome = "chr1",
      start = min(pos), end = max(pos),
      lead_snp = lead_snp, lead_p = lead_p,
      n_significant_markers = length(members),
      n_region_markers = length(members), status = "ok",
      stringsAsFactors = FALSE
    ),
    membership = data.frame(
      locus_id = locus_id, marker_id = members,
      significant_in_marginal_scan = TRUE,
      r2_to_lead = r2_to_lead, position = pos,
      stringsAsFactors = FALSE
    )
  )
}

# ---- Shared dataset A: two LINKED trait-specific causal variants -----------
# v1 -> T1, v2 -> T2, cor(v1, v2)^2 ~ 0.55-0.6: the marginal lead sees
# both traits (pseudo multi-trait), conditional resolution separates them.
set.seed(601)
n_a <- 500
x1 <- stats::rbinom(n_a, 2, 0.35)
x2 <- pmin(x1 + stats::rbinom(n_a, 1, 0.6), 2)
x3 <- stats::rbinom(n_a, 2, 0.3)
E_a <- mk_corr_noise(n_a, 0.4)
Y_a <- cbind(T1 = 1.0 * x1 + E_a[, 1], T2 = 1.0 * x2 + E_a[, 2])
G_a <- cbind(v1 = x1, v2 = x2, n1 = x3)
r2_linked <- stats::cor(x1, x2)^2
fit_a <- fit_mt_null(Y_a, K = diag(n_a), n_starts = 1L,
                     control = list(maxit = 200))
scan_a <- scan_mt_omnibus(fit_a, G_a)
locus_a <- mk_locus_object(
  "chr1:L1", "v1", scan_a$omnibus$p_value[scan_a$omnibus$marker_id == "v1"],
  members = c("v1", "v2", "n1"), pos = c(100, 200, 300),
  r2_to_lead = c(1, r2_linked, 0)
)
fb_a <- mk_recording_backend()
out_a <- pipelines(fit_a, G_a, locus_a, scan_a,
                   trait_names = c("T1", "T2"),
                   asset_backend = fb_a$backend)

# ---- Shared dataset B: one PLEIOTROPIC causal variant ----------------------
# c1 -> T1 and T2 (both survive within-signal Holm -> 2-trait candidate).
set.seed(602)
n_b <- 500
c1 <- stats::rbinom(n_b, 2, 0.3)
c2 <- pmin(c1 + stats::rbinom(n_b, 1, 0.5), 2)   # LD proxy
E_b <- mk_corr_noise(n_b, 0.4)
Y_b <- cbind(T1 = 1.0 * c1 + E_b[, 1], T2 = 0.9 * c1 + E_b[, 2])
G_b <- cbind(c1 = c1, c2 = c2)
fit_b <- fit_mt_null(Y_b, K = diag(n_b), n_starts = 1L,
                     control = list(maxit = 200))
scan_b <- scan_mt_omnibus(fit_b, G_b)
locus_b <- mk_locus_object(
  "chr1:L1", "c1", scan_b$omnibus$p_value[scan_b$omnibus$marker_id == "c1"],
  members = c("c1", "c2"), pos = c(100, 200),
  r2_to_lead = c(1, stats::cor(c1, c2)^2)
)
fb_b <- mk_recording_backend()
out_b <- pipelines(fit_b, G_b, locus_b, scan_b,
                   trait_names = c("T1", "T2"),
                   asset_backend = fb_b$backend)

# Hand-built correlated trait-z correlation matrix for direct adapter tests.
tn3 <- c("T1", "T2", "T3")
Sigma_Z3 <- matrix(c(1, .5, .2,
                     .5, 1, .3,
                     .2, .3, 1), 3, 3,
                   dimnames = list(tn3, tn3))
beta3 <- c(T1 = 1.0, T2 = 2.0, T3 = 0.5)
se3 <- c(T1 = 0.2, T2 = 0.3, T3 = 0.25)

# ---- Anchor 36: frozen output schema of the adapter ------------------------

test_that("Anchor 36: adapter output schema is exactly the frozen one", {
  fb <- mk_recording_backend()
  res <- asset_cmp(beta3, se = se3, Sigma_Z = Sigma_Z3,
                   trait_names = tn3,
                   backend = function(beta, se, z, Sigma_Z, trait_names,
                                      two_sided) {
                     list(p = 0.01, subset = c("T2", "T1"))
                   })
  expect_named(res, c("asset_p", "asset_best_subset",
                      "asset_positive_subset", "asset_negative_subset",
                      "status", "settings"))
  expect_identical(res$status, "ok")
  expect_identical(res$asset_p, 0.01)
  # best subset is normalised to canonical trait_names order
  expect_identical(res$asset_best_subset, c("T1", "T2"))
  expect_identical(res$asset_positive_subset, c("T1", "T2"))
  expect_identical(res$asset_negative_subset, character())

  # z may be supplied directly (no se)
  res_z <- asset_cmp(beta3, z = beta3 / se3, Sigma_Z = Sigma_Z3,
                     trait_names = tn3,
                     backend = function(...) list(p = 0.5, subset = "T3"))
  expect_identical(res_z$status, "ok")
  expect_identical(res_z$asset_best_subset, "T3")

  # data.frame raw shape (h.forest-style): one row per subset, p column
  # plus per-trait membership columns; smallest p wins
  raw_df <- data.frame(
    T1 = c(TRUE, TRUE, FALSE), T2 = c(TRUE, FALSE, TRUE),
    T3 = c(NA, NA, NA), p.value = c(0.30, 0.02, 0.40),
    stringsAsFactors = FALSE
  )
  res_df <- asset_cmp(beta3, se = se3, Sigma_Z = Sigma_Z3,
                      trait_names = tn3,
                      backend = function(...) raw_df)
  expect_identical(res_df$status, "ok")
  expect_identical(res_df$asset_p, 0.02)
  expect_identical(res_df$asset_best_subset, "T1")
})

test_that("adapter failure paths: asset_failed and invalid_input", {
  # backend throws -> asset_failed, never an error
  expect_no_error(
    res <- asset_cmp(beta3, se = se3, Sigma_Z = Sigma_Z3,
                     trait_names = tn3,
                     backend = function(...) stop("boom"))
  )
  expect_identical(res$status, "asset_failed")
  expect_true(is.na(res$asset_p))
  expect_identical(res$asset_best_subset, character())
  # backend returns garbage -> asset_failed
  res2 <- asset_cmp(beta3, se = se3, Sigma_Z = Sigma_Z3,
                    trait_names = tn3, backend = function(...) 42)
  expect_identical(res2$status, "asset_failed")
  # malformed inputs -> condped_invalid_input
  expect_error(asset_cmp(unname(beta3), se = se3, Sigma_Z = Sigma_Z3,
                         trait_names = tn3, backend = function(...) NULL),
               class = "condped_invalid_input")
  expect_error(asset_cmp(beta3, Sigma_Z = Sigma_Z3, trait_names = tn3,
                         backend = function(...) NULL),
               class = "condped_invalid_input")
  expect_error(asset_cmp(beta3, se = se3, Sigma_Z = diag(3),
                         trait_names = tn3, backend = function(...) NULL),
               NA) # diagonal Sigma_Z is legal (just not required)
  bad_sigma <- Sigma_Z3
  bad_sigma[1, 2] <- 0.9
  expect_error(asset_cmp(beta3, se = se3, Sigma_Z = bad_sigma,
                         trait_names = tn3, backend = function(...) NULL),
               class = "condped_invalid_input")
  expect_error(asset_cmp(beta3, se = se3, Sigma_Z = Sigma_Z3,
                         trait_names = c("T1", "T2"),
                         backend = function(...) NULL),
               class = "condped_invalid_input")
})

# ---- Anchor 37: correlated-trait input is honored ---------------------------

test_that("Anchor 37: full non-diagonal Sigma_Z reaches the backend unchanged", {
  fb <- mk_recording_backend()
  res <- asset_cmp(beta3, se = se3, Sigma_Z = Sigma_Z3,
                   trait_names = tn3, two_sided = TRUE,
                   backend = fb$backend)
  expect_identical(res$status, "ok")
  call <- fb$store$calls[[1L]]
  expect_identical(call$Sigma_Z, Sigma_Z3)          # arrives unchanged
  expect_true(any(call$Sigma_Z[row(call$Sigma_Z) != col(call$Sigma_Z)] != 0))
  expect_false(isTRUE(all.equal(call$Sigma_Z, diag(3))))
  expect_identical(call$trait_names, tn3)
  expect_identical(call$two_sided, TRUE)
  expect_equal(call$beta, beta3[tn3])
  expect_equal(call$z, (beta3 / se3)[tn3])
})

# ---- Anchor 38: two-sided direction/subset normalisation --------------------

test_that("Anchor 38: positive/negative subsets split by effect sign", {
  beta_mixed <- c(T1 = 1.5, T2 = -1.2, T3 = -0.7)
  backend <- function(...) list(p = 1e-4, subset = c("T3", "T1", "T2"))
  res <- asset_cmp(beta_mixed, se = se3, Sigma_Z = Sigma_Z3,
                   trait_names = tn3, two_sided = TRUE, backend = backend)
  expect_identical(res$asset_best_subset, tn3)      # canonical order
  expect_identical(res$asset_positive_subset, "T1")
  expect_identical(res$asset_negative_subset, c("T2", "T3"))

  fb <- mk_recording_backend()
  res1 <- asset_cmp(beta_mixed, se = se3, Sigma_Z = Sigma_Z3,
                    trait_names = tn3, two_sided = FALSE,
                    backend = fb$backend)
  expect_false(fb$store$calls[[1L]]$two_sided)      # flag forwarded
  # the sign split is by CondPED effect sign, regardless of sidedness
  expect_identical(res1$asset_positive_subset,
                   tn3[sign(beta_mixed[res1$asset_best_subset]) > 0])
  expect_identical(res1$asset_negative_subset,
                   tn3[sign(beta_mixed[res1$asset_best_subset]) < 0])
})

# ---- Anchor 39: ASSET missing -> asset_not_available, no error -------------

test_that("Anchor 39: backend NULL without ASSET gives asset_not_available", {
  skip_if(requireNamespace("ASSET", quietly = TRUE),
          "ASSET is installed; the not-available path does not apply")
  expect_no_error(
    res <- asset_cmp(beta3, se = se3, Sigma_Z = Sigma_Z3,
                     trait_names = tn3, backend = NULL)
  )
  expect_identical(res$status, "asset_not_available")
  expect_true(is.na(res$asset_p))
  expect_identical(res$asset_best_subset, character())
  expect_identical(res$asset_positive_subset, character())
  expect_identical(res$asset_negative_subset, character())

  # the whole orchestrator also degrades gracefully
  expect_no_error(
    out <- pipelines(fit_b, G_b, locus_b, scan_b,
                     trait_names = c("T1", "T2"), asset_backend = NULL)
  )
  expect_identical(out$lead_asset$status, "asset_not_available")
  expect_identical(out$resolved_asset$status, "asset_not_available")
  expect_identical(out$lead_asset$loci[[1L]]$status, "asset_not_available")
  expect_identical(out$condped_full$status, "ok")
})

# ---- Anchor 40: lead_asset uses only the marginal lead SNP -------------------

test_that("Anchor 40: lead_asset feeds the marginal lead-SNP effects", {
  lead_calls <- 1L   # lead pipeline runs first, one call per locus
  call <- fb_a$store$calls[[lead_calls]]
  mg <- estimate_mt_effects(fit_a, G_a, targets = "v1")
  el <- mg$effects_long
  expect_equal(call$beta, stats::setNames(el$beta, el$trait)[c("T1", "T2")])
  expect_equal(call$se, stats::setNames(el$se, el$trait)[c("T1", "T2")])
  # exactly one lead call per locus: no resolution involved in lead_asset
  expect_identical(length(out_a$lead_asset$loci), 1L)
  expect_identical(out_a$lead_asset$loci[[1L]]$lead_snp, "v1")
  expect_identical(out_a$lead_asset$loci[[1L]]$beta, call$beta)
  expect_identical(out_a$lead_asset$status, "ok")
  # total backend calls = 1 lead + n_resolved_signals
  expect_identical(length(fb_a$store$calls),
                   1L + out_a$resolved_asset$n_signals)
})

# ---- Anchor 41: resolved_asset uses final joint effects ---------------------

test_that("Anchor 41: resolved_asset feeds the final joint effects", {
  res_dir <- resolve_locus_signals(fit_a, G_a, locus_a)
  expect_identical(nrow(res_dir$signals), 2L)
  for (k in 1:2) {
    sid <- paste0("chr1:L1::S", k)
    b <- res_dir$beta[res_dir$beta$signal_id == sid, ]
    call <- fb_a$store$calls[[1L + k]]   # after the single lead call
    expect_equal(call$beta, stats::setNames(b$beta, b$trait)[c("T1", "T2")])
    expect_equal(call$se, stats::setNames(b$se, b$trait)[c("T1", "T2")])
  }
  # the joint effects differ from the marginals under LD
  b1 <- res_dir$beta[res_dir$beta$signal_id == "chr1:L1::S1", ]
  joint <- stats::setNames(b1$beta, b1$trait)[c("T1", "T2")]
  mg <- estimate_mt_effects(fit_a, G_a, targets = "v1")
  marg <- stats::setNames(mg$effects_long$beta,
                          mg$effects_long$trait)[c("T1", "T2")]
  expect_gt(max(abs(joint - marg)), 0.05)
  expect_identical(out_a$resolved_asset$status, "ok")
})

# ---- Anchor 42: condped_full returns eta/rho/Rep structures -----------------

test_that("Anchor 42: condped_full returns eta/rho/Rep structures", {
  cf <- out_b$condped_full
  expect_identical(cf$status, "ok")
  # the pleiotropic signal has a 2-trait candidate set
  expect_identical(cf$candidate_sets[["chr1:L1::S1"]], c("T1", "T2"))
  expect_identical(
    cf$signal_table$signal_status[cf$signal_table$signal_id == "chr1:L1::S1"],
    "multi_trait"
  )
  sa <- cf$subset_analysis
  # rho: representation loss per enumerated subset (2^2 = 4 rows)
  expect_gte(nrow(sa$subset_table), 4L)
  expect_true(all(c("representing_set", "representation_loss",
                    "conditional_effect") %in% names(sa$subset_table)))
  expect_true(all(is.finite(sa$subset_table$representation_loss)))
  # eta: conditional effects for the non-degenerate subsets
  eta_rows <- sa$subset_table$set_size %in% c(1L)
  expect_true(all(vapply(sa$subset_table$conditional_effect[eta_rows],
                         function(x) !is.null(x) && length(x) >= 1L,
                         logical(1))))
  # Rep: minimum representative sets + tolerance path are populated
  expect_gt(nrow(sa$minimum_representative_sets), 0L)
  expect_true(all(c("trait_set", "set_size", "representation_loss") %in%
                    names(sa$minimum_representative_sets)))
  expect_gt(nrow(sa$tolerance_path), 0L)
  expect_true("irreducible_modules" %in% names(sa))
})

# ---- Anchor 43: lead vs resolved differ under linked pseudo-pleiotropy ------

test_that("Anchor 43: linked trait-specific variants split lead vs resolved", {
  expect_true(r2_linked > 0.3 && r2_linked < 0.9)   # r2 ~ 0.6 by design
  # lead: the single marginal call sees BOTH traits (pseudo multi-trait)
  lead_res <- out_a$lead_asset$loci[["chr1:L1"]]
  expect_identical(lead_res$asset_best_subset, c("T1", "T2"))
  expect_gt(abs(fb_a$store$calls[[1L]]$z[["T2"]]), 2)
  # resolved: each signal's subset contains only its own trait
  ra <- out_a$resolved_asset$signals
  expect_identical(length(ra), 2L)
  expect_identical(ra[["chr1:L1::S1"]]$representative_snp, "v1")
  expect_identical(ra[["chr1:L1::S2"]]$representative_snp, "v2")
  expect_identical(ra[["chr1:L1::S1"]]$asset_best_subset, "T1")
  expect_identical(ra[["chr1:L1::S2"]]$asset_best_subset, "T2")
  # the two pipelines genuinely disagree on this locus
  expect_false(identical(lead_res$asset_best_subset,
                         ra[["chr1:L1::S1"]]$asset_best_subset))
})

# ---- pipeline-level wiring / validation --------------------------------------

test_that("orchestrator wires pipelines and validates inputs", {
  expect_identical(out_a$status, "ok")
  expect_named(out_a, c("lead_asset", "resolved_asset", "condped_full",
                        "status"))
  # resolved_asset and condped_full share ONE resolution: same signals
  expect_identical(sort(names(out_a$resolved_asset$signals)),
                   sort(out_a$condped_full$signal_table$signal_id))
  expect_error(pipelines(fit_a, G_a, locus_a, scan_a,
                         trait_names = c("T1", "WRONG"),
                         asset_backend = fb_a$backend),
               class = "condped_invalid_input")
  expect_error(pipelines(fit_a, unname(G_a), locus_a, scan_a,
                         trait_names = c("T1", "T2"),
                         asset_backend = fb_a$backend),
               class = "condped_invalid_input")
})

test_that("real ASSET::fast_asset integration smoke test", {
  testthat::skip_if_not_installed("ASSET")
  set.seed(99)
  tn <- paste0("T", 1:3)
  Sigma_Z <- matrix(c(1.0, 0.4, 0.2,
                      0.4, 1.0, 0.3,
                      0.2, 0.3, 1.0), 3, 3,
                    dimnames = list(tn, tn))
  res <- CondPED:::.run_asset_comparison(
    beta = c(T1 = 0.30, T2 = 0.25, T3 = -0.10),
    se = c(T1 = 0.10, T2 = 0.10, T3 = 0.10),
    Sigma_Z = Sigma_Z,
    trait_names = tn,
    sample_size = 1000
  )
  expect_identical(res$status, "ok")
  expect_true(is.finite(res$asset_p))
  expect_true(res$asset_p >= 0 && res$asset_p <= 1)
  # subsets come from the frozen schema and are valid traits
  expect_true(all(res$asset_best_subset %in% tn))
  expect_identical(sort(res$asset_best_subset),
                   sort(union(res$asset_positive_subset,
                              res$asset_negative_subset)))
  # settings record the screening threshold and Neff
  expect_equal(res$settings$scr_pthr, 0.05)
  expect_equal(res$settings$Neff, 1000)
  # real path requires se and sample_size
  expect_error(
    CondPED:::.run_asset_comparison(
      beta = c(T1 = 0.3, T2 = 0.2, T3 = 0.1),
      Sigma_Z = Sigma_Z, trait_names = tn, sample_size = 1000
    ),
    class = "condped_invalid_input"
  )
  expect_error(
    CondPED:::.run_asset_comparison(
      beta = c(T1 = 0.3, T2 = 0.2, T3 = 0.1),
      se = c(T1 = 0.1, T2 = 0.1, T3 = 0.1),
      Sigma_Z = Sigma_Z, trait_names = tn
    ),
    class = "condped_invalid_input"
  )
})

# ---- Sigma_Z source alignment (Stage 6C收口) ----------------------------------

.mk_ld_fixture <- function(n = 500L, seed = 301) {
  set.seed(seed)
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- pmin(x1 + stats::rbinom(n, 1, 0.1), 2)   # linked, r2 ~ 0.6-0.8
  G <- cbind(c1 = x1, c2 = x2,
             matrix(stats::rbinom(n * 3, 2, 0.3), n, 3,
                    dimnames = list(NULL, c("n1", "n2", "n3"))))
  Y <- cbind(x1 * 0.8, x2 * 0.8) +
    matrix(stats::rnorm(n * 2, sd = 0.6), n, 2,
           dimnames = list(NULL, c("T1", "T2")))
  fit <- fit_mt_null(Y, K = diag(n), n_starts = 2L,
                     control = list(maxit = 200L))
  scan <- scan_mt_omnibus(fit, G)
  locus <- list(
    loci = data.frame(locus_id = "chr1:1-5", chromosome = "chr1",
                      lead_snp = "c1", lead_p = 1e-9,
                      stringsAsFactors = FALSE),
    membership = data.frame(locus_id = "chr1:1-5",
                            marker_id = colnames(G),
                            position = seq_len(ncol(G)) * 1000,
                            stringsAsFactors = FALSE)
  )
  list(fit = fit, G = G, scan = scan, locus = locus, Y = Y)
}

test_that("A: lead_asset Sigma_Z equals cov2cor of the marginal covariance", {
  fx <- .mk_ld_fixture()
  out <- pipelines(fx$fit, fx$G, fx$locus, fx$scan,
                   trait_names = c("T1", "T2"),
                   asset_backend = function(...) {
                     list(p = 0.05, subset = c("T1"))
                   })
  lead_entry <- out$lead_asset$loci[[1]]
  marg <- estimate_mt_effects(fx$fit, fx$G, targets = "c1")
  expect_equal(lead_entry$Sigma_Z,
               stats::cov2cor(marg$covariance[, , 1]),
               tolerance = 1e-10)
})

test_that("B: resolved_asset Sigma_Z equals cov2cor of the joint block", {
  fx <- .mk_ld_fixture(seed = 302)
  res <- resolve_locus_signals(fx$fit, fx$G, fx$locus, max_signals = 2L)
  out <- pipelines(fx$fit, fx$G, fx$locus, fx$scan,
                   trait_names = c("T1", "T2"),
                   asset_backend = function(...) {
                     list(p = 0.05, subset = c("T1"))
                   })
  m <- 2L
  for (k in seq_along(res$signals$signal_id)) {
    sid <- res$signals$signal_id[k]
    ord <- res$signals$signal_order[k]
    lid <- res$signals$locus_id[k]
    joint_cov <- res$covariance[[lid]]
    idx <- (ord - 1L) * m + seq_len(m)
    expect_equal(out$resolved_asset$signals[[sid]]$Sigma_Z,
                 stats::cov2cor(joint_cov[idx, idx, drop = FALSE]),
                 tolerance = 1e-10, label = sid)
  }
})

test_that("C: pipelines use effect-covariance correlation, not Sigma_P_ref", {
  fx <- .mk_ld_fixture(seed = 303)
  out <- pipelines(fx$fit, fx$G, fx$locus, fx$scan,
                   trait_names = c("T1", "T2"),
                   asset_backend = function(...) {
                     list(p = 0.05, subset = c("T1"))
                   })
  ref_cor <- stats::cov2cor(fx$fit$Sigma_P_ref)
  marg <- estimate_mt_effects(fx$fit, fx$G, targets = "c1")
  eff_cor <- stats::cov2cor(marg$covariance[, , 1])
  # fixture must actually separate the two candidates
  expect_false(isTRUE(all.equal(ref_cor, eff_cor)))
  used <- out$lead_asset$loci[[1]]$Sigma_Z
  expect_equal(used, eff_cor, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(used, ref_cor)))
})

test_that("D: marker/trait order changes keep Sigma_Z and results aligned", {
  fx <- .mk_ld_fixture(seed = 304)
  backend <- function(...) list(p = 0.05, subset = c("T1"))
  out1 <- pipelines(fx$fit, fx$G, fx$locus, fx$scan,
                    trait_names = c("T1", "T2"), asset_backend = backend)
  perm <- c(3, 1, 4, 2, 5)
  G2 <- fx$G[, perm]
  scan2 <- scan_mt_omnibus(fx$fit, G2)
  locus2 <- fx$locus
  locus2$membership <- fx$locus$membership[
    match(colnames(G2), fx$locus$membership$marker_id), ]
  out2 <- pipelines(fx$fit, G2, locus2, scan2,
                    trait_names = c("T1", "T2"), asset_backend = backend)
  expect_equal(out1$lead_asset$loci[[1]]$Sigma_Z,
               out2$lead_asset$loci[[1]]$Sigma_Z, tolerance = 1e-12)
  expect_equal(out1$lead_asset$loci[[1]]$asset_p,
               out2$lead_asset$loci[[1]]$asset_p)
})

test_that("E: real fast_asset smoke through the pipeline stays ok", {
  testthat::skip_if_not_installed("ASSET")
  fx <- .mk_ld_fixture(seed = 305)
  out <- pipelines(fx$fit, fx$G, fx$locus, fx$scan,
                   trait_names = c("T1", "T2"), asset_backend = NULL)
  expect_identical(out$lead_asset$loci[[1]]$status, "ok")
  expect_true(is.finite(out$lead_asset$loci[[1]]$asset_p))
  expect_true(all(out$lead_asset$loci[[1]]$asset_best_subset %in%
                    c("T1", "T2")))
})
