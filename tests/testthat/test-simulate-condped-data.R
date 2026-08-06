# Tests for simulate_condped_data() (v1.0 Stage 4).
#
# All checks read population truth only; none of the acceptance logic
# may depend on sample statistics (frozen rule).

sim_archs <- c("null", "candidate_single", "candidate_pair",
               "candidate_dense", "representative_singleton",
               "representative_pair", "representative_full",
               "irreducible_singleton", "irreducible_pair",
               "multiple_modules")

sim <- function(architecture, seed = 1, ...) {
  simulate_condped_data(n = 80, m = 4L, p = 150L,
                        architecture = architecture,
                        correlation = "block", seed = seed, ...)
}

# ---- contract surface ---------------------------------------------------------

test_that("architecture enumeration matches the contract exactly", {
  expect_identical(eval(formals(simulate_condped_data)$architecture),
                   sim_archs)
})

test_that("truth uses only v1.0 field names", {
  s <- sim("candidate_pair")
  expect_identical(
    names(s$truth),
    c("beta", "candidate_traits", "subset_table",
      "minimum_representative_sets", "irreducible_modules",
      "tolerance_path", "locus_pve", "target_loss", "delta",
      "Sigma_causal", "Sigma_G_total", "Sigma_G_bg", "Sigma_E",
      "Sigma_P_total")
  )
  expect_null(s$latent)          # return_latent = FALSE is the default
  expect_false(any(c("D", "eta", "gamma", "C", "conditional_profile",
                     "conditional_partitions", "conditional_components",
                     "B_Q") %in% names(s$truth)))
})

test_that("seed reproduction is value-identical; different seeds differ", {
  s1 <- sim("representative_singleton", seed = 5)
  s2 <- sim("representative_singleton", seed = 5)
  expect_identical(s1$Y, s2$Y)
  expect_identical(s1$truth$beta, s2$truth$beta)
  s3 <- sim("representative_singleton", seed = 6)
  expect_false(isTRUE(all.equal(s1$Y, s3$Y)))
})

test_that("candidate truth equals the nonzero-beta traits", {
  for (a in sim_archs[-1]) {
    s <- sim(a, seed = 3)
    beta <- s$truth$beta[1, ]
    expect_identical(s$truth$candidate_traits,
                     unname(paste0("Trait", seq_len(4))[beta != 0]),
                     label = a)
  }
  s0 <- sim("null")
  expect_identical(s0$truth$candidate_traits, character())
  expect_true(all(s0$truth$beta == 0))
})

test_that("causal variant stays out of K_bg", {
  s <- simulate_condped_data(n = 400, m = 4L, p = 200L,
                             architecture = "candidate_single",
                             correlation = "independent",
                             structured = TRUE, seed = 9)
  zc <- as.vector(scale(s$G[, s$causal_index]))
  outer_z <- outer(zc, zc)
  # if the causal marker leaked into the GRM, K_bg would correlate with
  # zc %o% zc at roughly 1/p; excluded it stays near zero
  cc <- stats::cor(K_vec <- s$K_bg[upper.tri(s$K_bg)],
                   outer_z[upper.tri(outer_z)])
  expect_lt(abs(cc), 0.05)
  expect_equal(mean(diag(s$K_bg)), 1, tolerance = 1e-8)
})

test_that("Sigma_G_bg is PSD and PVE hits the target", {
  for (a in sim_archs[-1]) {
    s <- sim(a, seed = 4)
    expect_gt(CondPED:::.min_eigen_sym(s$truth$Sigma_G_bg), -1e-8)
    lp <- s$truth$locus_pve[s$truth$candidate_traits]
    expect_equal(mean(lp), 0.01, tolerance = 1e-8, label = a)
    # Sigma_P = Sigma_G_bg + Sigma_causal + Sigma_E
    expect_equal(s$truth$Sigma_P_total,
                 s$truth$Sigma_G_bg + s$truth$Sigma_causal + s$truth$Sigma_E,
                 tolerance = 1e-10)
  }
})

test_that("PVE scaling leaves the structural truth unchanged", {
  s1 <- sim("representative_pair", seed = 8, locus_pve = 0.01)
  s2 <- sim("representative_pair", seed = 8, locus_pve = 0.05)
  expect_equal(s1$truth$subset_table$representation_loss,
               s2$truth$subset_table$representation_loss,
               tolerance = 1e-8)
  expect_identical(s1$truth$minimum_representative_sets$trait_key,
                   s2$truth$minimum_representative_sets$trait_key)
  expect_identical(s1$truth$irreducible_modules$trait_key,
                   s2$truth$irreducible_modules$trait_key)
  # but the effect vector really was rescaled
  expect_gt(norm(s2$truth$beta, "2"), 2 * norm(s1$truth$beta, "2"))
})

# ---- closed-form construction (Methods eqs. 39-41) -----------------------------

test_that("closed form: target set has population loss exactly target_loss", {
  for (a in c("representative_singleton", "representative_pair")) {
    s <- sim(a, seed = 11, target_loss = 0.05)
    key <- paste(s$generator$settings$target_representative_set,
                 collapse = "|")
    tab <- s$truth$subset_table
    expect_equal(tab$representation_loss[tab$representing_key == key],
                 0.05, tolerance = 1e-8, label = a)
    expect_equal(s$truth$target_loss, 0.05)
    expect_true(is.finite(s$truth$delta))
  }
})

# ---- per-scene acceptance conditions --------------------------------------------

test_that("representative_singleton: target is the minimum feasible set", {
  s <- sim("representative_singleton", seed = 21)
  reps <- s$truth$minimum_representative_sets
  reps <- reps[reps$tolerance == 0.10, ]
  expect_identical(reps$set_size[1], 1L)
  expect_true("Trait1" %in% reps$trait_key)
  # no smaller (empty) set is ever feasible: loss(empty) = 1 > tau
  tab <- s$truth$subset_table
  expect_gt(tab$representation_loss[tab$representing_key == "<empty>"], 0.10)
})

test_that("representative_pair: pair feasible, all singletons infeasible", {
  s <- sim("representative_pair", seed = 22)
  tab <- s$truth$subset_table
  singles <- tab[tab$set_size == 1L, ]
  expect_true(all(singles$representation_loss > 0.10))
  reps <- s$truth$minimum_representative_sets
  reps <- reps[reps$tolerance == 0.10, ]
  expect_true("Trait1|Trait2" %in% reps$trait_key)
  expect_identical(reps$set_size[1], 2L)
})

test_that("representative_full: every proper subset is infeasible", {
  s <- sim("representative_full", seed = 23)
  tab <- s$truth$subset_table
  proper <- tab[tab$set_size < 4L, ]
  expect_true(all(proper$representation_loss > 0.10))
  reps <- s$truth$minimum_representative_sets
  reps <- reps[reps$tolerance == 0.10, ]
  expect_identical(reps$trait_key, "Trait1|Trait2|Trait3|Trait4")
})

test_that("irreducible_singleton: deleting the target trait is infeasible", {
  s <- sim("irreducible_singleton", seed = 24)
  tab <- s$truth$subset_table
  comp <- paste(setdiff(s$truth$candidate_traits, "Trait1"), collapse = "|")
  expect_gt(tab$representation_loss[tab$representing_key == comp], 0.10)
  mods <- s$truth$irreducible_modules
  expect_true("Trait1" %in% mods$trait_key[mods$tolerance == 0.10])
})

test_that("irreducible_pair: joint deletion infeasible, single deletions feasible", {
  s <- sim("irreducible_pair", seed = 25)
  tab <- s$truth$subset_table
  loss_of <- function(key) tab$representation_loss[tab$representing_key == key]
  A <- s$truth$candidate_traits
  expect_gt(loss_of(paste(setdiff(A, c("Trait1", "Trait2")), collapse = "|")),
            0.10)
  expect_lte(loss_of(paste(setdiff(A, "Trait1"), collapse = "|")), 0.10)
  expect_lte(loss_of(paste(setdiff(A, "Trait2"), collapse = "|")), 0.10)
  mods <- s$truth$irreducible_modules
  expect_true("Trait1|Trait2" %in% mods$trait_key[mods$tolerance == 0.10])
})

test_that("multiple_modules: extraction matches the preset exactly", {
  s <- sim("multiple_modules", seed = 26)
  mods <- s$truth$irreducible_modules
  expect_identical(sort(mods$trait_key[mods$tolerance == 0.10]),
                   sort(c("Trait1", "Trait3|Trait4")))
})

test_that("truth_margin keeps every loss away from the boundary", {
  for (a in c("representative_singleton", "representative_pair",
              "representative_full", "irreducible_singleton",
              "irreducible_pair", "multiple_modules")) {
    s <- sim(a, seed = 30)
    losses <- s$truth$subset_table$representation_loss
    expect_true(all(abs(losses - 0.10) >= 0.02), label = a)
  }
  # margin that target_loss violates is rejected at the interface
  expect_error(sim("representative_singleton", target_loss = 0.09,
                   truth_margin = 0.02),
               class = "condped_invalid_input")
})

test_that("null scene: empty truth, zero effects, legal structure", {
  s <- sim("null", seed = 31)
  expect_true(s$status$ok)
  expect_identical(nrow(s$truth$subset_table), 0L)
  expect_identical(nrow(s$truth$minimum_representative_sets), 0L)
  expect_identical(nrow(s$truth$irreducible_modules), 0L)
  expect_identical(s$truth$locus_pve, setNames(rep(0, 4),
                                               paste0("Trait", 1:4)))
})

test_that("input validation", {
  expect_error(simulate_condped_data(architecture = "conditional_deviation"),
               class = "error")   # removed v0.3 architecture
  expect_error(simulate_condped_data(n = 1, architecture = "null"),
               class = "condped_invalid_input")
  expect_error(simulate_condped_data(n = 80, architecture = "null",
                                     tolerance = 1),
               class = "condped_invalid_input")
  expect_error(sim("null", target_loss = -0.1),
               class = "condped_invalid_input")
  expect_error(sim("null", truth_margin = -1),
               class = "condped_invalid_input")
  expect_error(simulate_condped_data(n = 80, m = 3L,
                                     architecture = "multiple_modules"),
               class = "condped_invalid_input")
})
