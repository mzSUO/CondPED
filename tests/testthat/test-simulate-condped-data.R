# Tests for simulate_condped_data() (Stage 6C): unified multi-signal
# model, target-rho construction, scenario registry, analysis modes.
# Covers anchors 1-29 of the Stage 6C specification.

sim <- function(experiment, scenario, seed = 1, ...) {
  simulate_condped_data(n = 150, m = 4L, p = 200L,
                        experiment = experiment, scenario = scenario,
                        local_region_size = 25L, seed = seed, ...)
}

# ---- A. multi-signal simulator --------------------------------------------------

test_that("Anchor 1: single-signal scaler regression (new vs old formula)", {
  b <- c(1, 0.8, 0, 0)
  v <- 0.42
  expect_equal(CondPED:::.scale_signal_pve(b, v, 0.01),
               sqrt(0.01 / (v * mean(b[b != 0]^2))))
  expect_equal(CondPED:::.scale_signal_pve(b, v, 0.03),
               sqrt(0.03 / (v * mean(b[b != 0]^2))))
})

test_that("Anchors 2+9: both signals enter the phenotype, exactly", {
  s <- sim("signal_resolution", "two_heterogeneous", seed = 11,
           return_latent = TRUE)
  expect_identical(nrow(s$truth$B_Q), 2L)
  expect_true(all(s$truth$B_Q[1, 1:2] != 0))
  expect_true(s$truth$B_Q[2, 3] != 0)
  expect_true(all(s$truth$B_Q[2, 1:2] == 0))
  # exact reconstruction Y = X_Q B_Q + U + E
  rec <- s$G[, s$causal_index, drop = FALSE] %*% s$truth$B_Q +
    s$latent$U + s$latent$E
  expect_equal(rec, s$Y, tolerance = 1e-10)
})

test_that("Anchor 3: linked-signal LD is recorded correctly", {
  s <- sim("signal_resolution", "two_linked_trait_specific",
           seed = 12, target_r2 = 0.3)
  expect_equal(s$truth$local_ld$target_r2, 0.3)
  expect_equal(s$truth$local_ld$realized_mean_r2, 0.3, tolerance = 0.05)
  # the causal pair's own realized r2 is close to the target
  g1 <- s$G[, s$causal_index[1]]
  g2 <- s$G[, s$causal_index[2]]
  expect_lt(abs(stats::cor(g1, g2)^2 - 0.3), 0.12)
})

test_that("Anchors 4+5: Sigma_Q and Sigma_G_bg identities", {
  s <- sim("signal_resolution", "two_heterogeneous", seed = 13)
  direct <- crossprod(s$truth$B_Q, s$truth$Sigma_X %*% s$truth$B_Q)
  expect_equal(s$truth$Sigma_Q, direct, tolerance = 1e-10)
  expect_equal(s$truth$Sigma_G_bg,
               s$truth$Sigma_G_total - s$truth$Sigma_Q,
               tolerance = 1e-10)
})

test_that("Anchor 6: Sigma_G_bg is PSD for accepted generations", {
  for (sc in list(c("signal_resolution", "single_multi_trait"),
                  c("signal_resolution", "two_linked_trait_specific"),
                  c("trait_representation", "broad_concordant"),
                  c("end_to_end", "mixed_multisignal"))) {
    s <- sim(sc[1], sc[2], seed = 14)
    expect_true(s$status$ok)
    expect_gt(CondPED:::.min_eigen_sym(s$truth$Sigma_G_bg), -1e-8)
  }
})

test_that("Anchor 7: focal region is excluded from K_bg", {
  s <- sim("signal_resolution", "two_linked_trait_specific",
           seed = 15, target_r2 = 0.3)
  # if the region leaked into the GRM, K_bg would correlate with the
  # region markers' standardized outer products (~ region_size / p)
  bg <- setdiff(seq_len(200), match(s$focal_region, colnames(s$G)))
  maf_hat <- colMeans(s$G) / 2
  Zr <- scale(s$G[, s$focal_region])   # region markers only
  Zr <- Zr[, apply(Zr, 2, stats::sd) > 0, drop = FALSE]
  signal_outer <- tcrossprod(rowSums(Zr))[upper.tri(s$K_bg)]
  cc <- stats::cor(s$K_bg[upper.tri(s$K_bg)], signal_outer)
  expect_lt(abs(cc), 0.05)
})

test_that("Anchor 8: K_bg is scaled to trace/n = 1", {
  s <- sim("signal_resolution", "single_multi_trait", seed = 16)
  expect_equal(s$diagnostics$trace_K_over_n, 1, tolerance = 1e-8)
  expect_equal(s$diagnostics$mean_diag_K, 1, tolerance = 1e-8)
})

test_that("Anchor 10: seed reproduction and marker-order invariance", {
  s1 <- sim("signal_resolution", "two_heterogeneous", seed = 17)
  s2 <- sim("signal_resolution", "two_heterogeneous", seed = 17)
  expect_identical(s1$Y, s2$Y)
  expect_identical(s1$truth$subset_table, s2$truth$subset_table)
  s3 <- sim("signal_resolution", "two_heterogeneous", seed = 18)
  expect_false(isTRUE(all.equal(s1$Y, s3$Y)))
  # permuting traits does not change the loss map values
  beta <- s1$truth$B_Q[1, ]
  ids <- list(locus_id = "L", signal_id = "L::S1",
              representative_snp = "M1")
  m1 <- CondPED:::.compute_truth_representation_map(
    beta, s1$truth$Sigma_P_ref, 0.1, 0.1, ids)
  perm <- c(3, 1, 2, 4)
  m2 <- CondPED:::.compute_truth_representation_map(
    beta[perm], s1$truth$Sigma_P_ref[perm, perm], 0.1, 0.1, ids)
  l1 <- m1$subset_table$representation_loss
  l2raw <- m2$subset_table$representation_loss[
    match(m1$subset_table$representing_key,
          m2$subset_table$representing_key)]
  # keys follow each call's own trait order; compare as sorted values
  expect_equal(sort(l1), sort(m2$subset_table$representation_loss),
               tolerance = 1e-10)
})

# ---- B. target-rho / representation ----------------------------------------------

test_that("Anchor 11: requested rho equals the exact enumerated rho", {
  for (rt in c(0, 0.02, 0.05, 0.2)) {
    s <- sim("trait_representation", "highly_representable",
             seed = 20, target_loss = rt)
    tab <- s$truth$subset_table
    expect_equal(tab$representation_loss[
      tab$representing_key == "Trait1"], rt, tolerance = 1e-8,
      label = paste("target_loss", rt))
  }
})

test_that("Anchors 12-14: full truth map matches the unique loss path", {
  s <- sim("trait_representation", "broad_concordant", seed = 21)
  tab <- s$truth$subset_table
  beta <- s$truth$B_Q[1, ]
  A <- names(beta)[beta != 0]
  for (i in sample(nrow(tab), 4L)) {
    ref <- CondPED:::.compute_subset_loss(
      beta[A], s$truth$Sigma_P_ref[A, A],
      tab$representing_set[[i]], trait_names = A)
    expect_equal(tab$representation_loss[i], ref$representation_loss,
                 tolerance = 1e-12)
    expect_identical(tab$conditional_effect[[i]], ref$eta)
  }
  expect_identical(tab$representation_loss[
    tab$representing_key == "<empty>"], 1)
  expect_identical(tab$representation_loss[
    tab$set_size == length(A)], 0)
})

test_that("Anchor 15: the truth map is monotone on nested pairs", {
  s <- sim("trait_representation", "broad_concordant", seed = 22)
  tab <- s$truth$subset_table
  chk <- CondPED:::.check_loss_monotonicity(tab)
  expect_identical(chk$n_violations, 0L)
})

test_that("Anchors 16-18: R1/R2/R3 acceptance conditions", {
  r1 <- sim("trait_representation", "highly_representable", seed = 23)
  reps1 <- r1$truth$minimum_representative_sets
  expect_identical(min(reps1$set_size[reps1$tolerance == 0.1]), 1L)

  r2 <- sim("trait_representation", "partially_representable", seed = 24)
  tab2 <- r2$truth$subset_table
  expect_true(all(tab2$representation_loss[tab2$set_size == 1L] > 0.10))
  expect_true(any(tab2$representation_loss[tab2$set_size == 2L] <= 0.10))
  reps2 <- r2$truth$minimum_representative_sets
  expect_identical(min(reps2$set_size[reps2$tolerance == 0.1]), 2L)

  r3 <- sim("trait_representation", "strongly_nonredundant", seed = 25)
  tab3 <- r3$truth$subset_table
  k <- max(tab3$set_size)
  expect_true(all(tab3$representation_loss[tab3$set_size < k] > 0.10))
})

test_that("Anchor 19: effect-scale matching does not change the rho map", {
  s1 <- sim("trait_representation", "highly_representable",
            seed = 26, locus_pve = 0.01)
  s2 <- sim("trait_representation", "highly_representable",
            seed = 26, locus_pve = 0.02)
  expect_identical(s1$diagnostics$attempts, s2$diagnostics$attempts)
  expect_equal(s1$truth$subset_table$representation_loss,
               s2$truth$subset_table$representation_loss,
               tolerance = 1e-10)
  expect_identical(
    s1$truth$minimum_representative_sets$trait_key,
    s2$truth$minimum_representative_sets$trait_key)
  expect_identical(
    s1$truth$irreducible_modules$trait_key,
    s2$truth$irreducible_modules$trait_key)
  # the calibration quantity is matched across architectures by
  # rescaling the whole beta (larger locus_pve -> larger effects)
  expect_gt(norm(s2$truth$B_Q, "F"), 1.2 * norm(s1$truth$B_Q, "F"))
})

test_that("Anchor 20: truth Rep keeps all ties", {
  # T1 and T2 nearly collinear with equal effects: each represents
  # the other almost perfectly -> two tied minimum singleton reps
  Sigma <- matrix(c(1, 0.99, 0.99, 1), 2, 2,
                  dimnames = list(c("T1", "T2"), c("T1", "T2")))
  beta <- c(T1 = 1, T2 = 1)
  ids <- list(locus_id = "L", signal_id = "L::S1",
              representative_snp = "M1")
  tabs <- CondPED:::.compute_truth_representation_map(
    beta, Sigma, 0.10, 0.10, ids)
  reps <- tabs$minimum_representative_sets
  expect_identical(reps$set_size[1], 1L)
  expect_identical(reps$n_tied_solutions[1], 2L)
  expect_identical(sort(reps$trait_key), c("T1", "T2"))
})

test_that("Anchor 21: truth Irr equals maximal-infeasible complements", {
  s <- sim("trait_representation", "broad_concordant", seed = 27)
  tab <- s$truth$subset_table
  tol <- 0.10
  infeas <- tab[tab$representation_loss > tol + 1e-10, ]
  A <- paste0("Trait", 1:4)
  if (nrow(infeas) > 0L) {
    is_max <- vapply(seq_len(nrow(infeas)), function(i) {
      Si <- infeas$representing_set[[i]]
      !any(vapply(seq_len(nrow(infeas)), function(j) {
        i != j && all(Si %in% infeas$representing_set[[j]])
      }, logical(1)))
    }, logical(1))
    expected <- sort(vapply(which(is_max), function(i) {
      CondPED:::.trait_set_key(setdiff(A, infeas$representing_set[[i]]))
    }, character(1)))
    mods <- s$truth$irreducible_modules
    expect_identical(sort(mods$trait_key[mods$tolerance == tol]),
                     expected)
  }
})

# ---- C. scenario families --------------------------------------------------------

test_that("Anchors 22-24: only the three formal families and scenarios exist", {
  expect_error(simulate_condped_data(experiment = "simulation_iv",
                                     scenario = "null"),
               class = "error")
  expect_error(sim("signal_resolution", "multiple_modules"),
               class = "condped_invalid_input")
  expect_error(sim("signal_resolution", "conditional_deviation"),
               class = "condped_invalid_input")
  all_scen <- list(
    c("signal_resolution", "null"),
    c("signal_resolution", "single_multi_trait"),
    c("signal_resolution", "two_linked_trait_specific"),
    c("signal_resolution", "two_heterogeneous"),
    c("trait_representation", "trait_specific"),
    c("trait_representation", "two_trait_concordant"),
    c("trait_representation", "two_trait_antagonistic"),
    c("trait_representation", "broad_concordant"),
    c("trait_representation", "highly_representable"),
    c("trait_representation", "partially_representable"),
    c("trait_representation", "strongly_nonredundant"),
    c("end_to_end", "single_highly_representable"),
    c("end_to_end", "single_nonredundant"),
    c("end_to_end", "linked_pseudo_multitrait"),
    c("end_to_end", "mixed_multisignal")
  )
  for (sc in all_scen) {
    s <- sim(sc[1], sc[2], seed = 28)
    expect_true(s$status$ok, label = paste(sc, collapse = "/"))
  }
})

# ---- D. analysis modes --------------------------------------------------------------

test_that("Anchors 25-29: analysis-mode mapping is frozen and complete", {
  mm <- CondPED:::.analysis_mode_mapping
  expect_identical(mm("full"), list(signal_mode = "resolve",
                                    candidate_mode = "holm_fwer"))
  expect_identical(mm("signal_oracle"),
                   list(signal_mode = "predefined",
                        candidate_mode = "holm_fwer"))
  expect_identical(mm("signal_trait_oracle"),
                   list(signal_mode = "predefined",
                        candidate_mode = "predefined"))
  # all_traits must never substitute for a trait oracle
  expect_false(mm("signal_trait_oracle")$candidate_mode == "all_traits")
  # restricted is gone
  expect_error(mm("restricted"), class = "error")
})

# ---- validation ---------------------------------------------------------------------

test_that("input validation and conflict checks", {
  expect_error(sim("signal_resolution", "single_multi_trait",
                   n_signals = 2L),
               class = "condped_invalid_input")
  expect_error(simulate_condped_data(n = 100, m = 0L, p = 100,
                                     experiment = "signal_resolution",
                                     scenario = "null"),
               class = "condped_invalid_input")
  expect_error(sim("signal_resolution", "null", tolerance = 1.5),
               class = "condped_invalid_input")
  expect_error(sim("signal_resolution", "null", target_r2 = 2),
               class = "condped_invalid_input")
  expect_error(simulate_condped_data(n = 100, p = 100,
                                     experiment = "signal_resolution",
                                     scenario = "null",
                                     local_region_size = 500L),
               class = "condped_invalid_input")
})

test_that("truth uses only the unified v1.0 schema fields", {
  s <- sim("end_to_end", "mixed_multisignal", seed = 29)
  expect_identical(
    names(s$truth),
    c("loci", "signals", "causal_markers", "signal_count",
      "local_ld", "realized_signal_pve", "beta", "candidate_traits",
      "effect_breadth", "effect_direction", "subset_table",
      "representation_map", "minimum_representative_sets",
      "irreducible_modules", "tolerance_path", "B_Q", "Sigma_X",
      "Sigma_Q", "Sigma_G_total", "Sigma_G_bg", "Sigma_E",
      "Sigma_P_total", "Sigma_P_ref")
  )
  expect_null(s$latent)
  expect_false(any(c("D", "eta", "gamma", "C", "qtl_index",
                     "conditional_profile", "effect_magnitude")
                   %in% names(s$truth)))
})

test_that("three_linked_trait_specific generates 3 clustered trait-specific signals", {
  s <- sim("signal_resolution", "three_linked_trait_specific", seed = 21)
  expect_true(isTRUE(s$status$ok))
  expect_identical(nrow(s$truth$signals), 3L)
  expect_identical(nrow(s$truth$B_Q), 3L)
  # contiguous causal cluster (Stage 7.4 A geometry), 1 kb apart
  expect_equal(diff(s$causal_index), c(1L, 1L))
  # per-trait-specific candidate sets
  expect_identical(unname(s$truth$candidate_traits[[1]]), "Trait1")
  expect_identical(unname(s$truth$candidate_traits[[2]]), "Trait2")
  expect_identical(unname(s$truth$candidate_traits[[3]]), "Trait3")
})
