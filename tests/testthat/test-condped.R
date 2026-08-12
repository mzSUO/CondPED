# Tests for condped() (Stage 6B-5): signal-level integration.

# Build a two-signal locus: S1 strong marginal; S2 marginally masked by
# LD (suppression) but conditionally strong.
mk_two_signal <- function(n = 600L, m = 3L, seed = 1L) {
  set.seed(seed)
  x1 <- stats::rbinom(n, 2, 0.3)
  x2 <- x1 + stats::rbinom(n, 1, 0.25)
  x2 <- pmin(x2, 2)                       # cor(x1, x2) ~ 0.7
  rho <- stats::cor(x1, x2)
  # joint effects chosen so S2's MARGINAL omnibus is ~0 on every
  # trait (suppression through LD), while its conditional effects are
  # strong on traits 1-2 and weak on trait 3
  b1 <- c(1.0, 0.9, 0.3)
  b2 <- c(-rho, -0.9 * rho, -0.3 * rho + 0.02)
  G <- cbind(S1 = x1, S2 = x2,
             matrix(stats::rbinom(n * 150, 2, 0.3), n, 150))
  colnames(G)[3:ncol(G)] <- paste0("N", 1:150)
  B <- matrix(0, ncol(G), m)
  B[1, ] <- b1
  B[2, ] <- b2
  Y <- G %*% B + matrix(stats::rnorm(n * m, sd = 0.9), n, m,
                        dimnames = list(NULL, paste0("Trait", 1:m)))
  pos <- c(10000, 10500, seq(20000, length.out = 150, by = 1000))
  chr <- rep("chr1", ncol(G))
  list(Y = Y, G = G, position = pos, chromosome = chr)
}

run_cp <- function(dat, control = list(), ...) {
  condped(dat$Y, G = dat$G, K = NULL, chromosome = dat$chromosome,
          position = dat$position,
          control = utils::modifyList(
            list(null_control = list(maxit = 300L), window_bp = 5000),
            control),
          ...)
}

test_that("Anchor 1: one-signal end-to-end", {
  set.seed(21)
  n <- 400; m <- 3
  G <- matrix(stats::rbinom(n * 120, 2, 0.3), n, 120)
  colnames(G) <- paste0("M", 1:120)
  Y <- G[, 1] %*% t(c(1, 0.8, 0)) +
    matrix(stats::rnorm(n * m, sd = 0.9), n, m,
           dimnames = list(NULL, paste0("Trait", 1:m)))
  f <- condped(Y, G = G, K = NULL,
               position = seq_len(120) * 1000,
               control = list(null_control = list(maxit = 300L),
                              window_bp = 3000))
  expect_s3_class(f, "condped_fit")
  expect_identical(f$status$code, "ok")
  expect_identical(
    names(f),
    c("null_fit", "omnibus", "loci", "signals", "effects",
      "candidate_traits", "signal_summary", "subset_analysis",
      "pve", "crossfit", "bootstrap", "settings", "status",
      "diagnostics")
  )
  expect_identical(nrow(f$signals), 1L)
  expect_true(all(c("p_adjusted", "selected") %in% names(f$omnibus$omnibus)))
  expect_null(f$pve)
})

test_that("Anchors 2+8+9+10: two-signal locus, joint betas downstream, ids", {
  dat <- mk_two_signal()
  f <- run_cp(dat)
  sigs <- f$signals$representative_snp
  expect_true(all(c("S1", "S2") %in% sigs))
  # downstream beta equals the final joint model, signal by signal
  joint <- CondPED:::.fit_joint_signal_effects(f$null_fit, dat$G, sigs)
  for (s in sigs) {
    est_b <- f$effects$beta$beta[f$effects$beta$representative_snp == s]
    expect_equal(est_b, unname(joint$beta[s, ]), tolerance = 1e-10,
                 label = s)
  }
  # each signal gets its OWN eta/rho map and Rep/Irr rows
  tab <- f$subset_analysis$subset_table
  sid <- unique(tab$signal_id)
  expect_identical(length(sid), length(sigs))
  st <- f$candidate_traits$signal_table
  for (s in sid) {
    sub <- tab[tab$signal_id == s, ]
    k <- st$n_candidates[st$signal_id == s]
    expect_equal(nrow(sub), 2L^k)       # own rho map, per-signal
    expect_length(unique(sub$locus_id), 1L)
    expect_true(all(c("locus_id", "signal_id") %in% names(sub)))
  }
  reps <- f$subset_analysis$minimum_representative_sets
  mods <- f$subset_analysis$irreducible_modules
  path <- f$subset_analysis$tolerance_path
  for (df in list(reps, mods, path)) {
    expect_true(all(c("locus_id", "signal_id") %in% names(df)))
  }
  # full ID chain: signals -> effects -> candidate_traits -> subset
  expect_true(all(c("locus_id", "signal_id", "representative_snp",
                    "trait", "beta", "se") %in% names(f$effects$beta)))
  expect_true(all(c("locus_id", "signal_id", "representative_snp",
                    "trait") %in% names(f$candidate_traits$trait_table)))
  expect_setequal(unique(f$effects$beta$signal_id),
                  f$signals$signal_id)
  expect_setequal(unique(f$candidate_traits$trait_table$signal_id),
                  f$signals$signal_id)
  # every signal with >=1 candidate appears in the subset analysis
  expect_setequal(unique(tab$signal_id), f$signals$signal_id)
  # no cross-contamination: candidate sets are per-signal
  cs <- f$candidate_traits$candidate_sets
  expect_true(all(grepl("::S", names(cs))))
})

test_that("Anchor 3: marginally masked secondary signal still enters Holm", {
  dat <- mk_two_signal(seed = 3)
  f <- run_cp(dat)
  om <- f$omnibus$omnibus
  # S2 is marginally masked: it never enters the omnibus selection
  expect_false("S2" %in% om$marker_id[om$selected])
  # but conditional resolution finds it and it goes through Holm
  expect_true("S2" %in% f$signals$representative_snp)
  sid_S2 <- f$signals$signal_id[f$signals$representative_snp == "S2"]
  expect_true(sid_S2 %in% f$candidate_traits$signal_table$signal_id)
  expect_gt(f$candidate_traits$signal_table$n_candidates[
    f$candidate_traits$signal_table$signal_id == sid_S2], 0L)
})

test_that("Anchor 4+6+7: within-signal Holm, breadth, direction", {
  dat <- mk_two_signal(seed = 4)
  f <- run_cp(dat)
  tt <- f$candidate_traits$trait_table
  expect_true(all(tt$adjust_scope == "within_signal"))
  for (s in unique(tt$signal_id)) {
    rows <- tt[tt$signal_id == s, ]
    ok <- !is.na(rows$p_raw)
    expect_equal(rows$p_adjusted[ok],
                 stats::p.adjust(rows$p_raw[ok], "holm"), label = s)
  }
  st <- f$candidate_traits$signal_table
  expect_identical(st$effect_breadth, st$n_candidates)
  expect_true(all(st$direction_pattern %in%
                    c("not_applicable", "single_trait", "concordant",
                      "antagonistic", "mixed")))
  # per-signal directions are evaluated independently: S1's candidate
  # effects are all positive, S2's all negative
  s1 <- f$signals$signal_id[f$signals$representative_snp == "S1"]
  s2 <- f$signals$signal_id[f$signals$representative_snp == "S2"]
  cand1 <- tt[tt$signal_id == s1 & tt$candidate, ]
  cand2 <- tt[tt$signal_id == s2 & tt$candidate, ]
  expect_true(all(cand1$effect_sign == 1))
  expect_true(all(cand2$effect_sign == -1))
  # both are internally concordant; direction_pattern values are legal
  expect_true(all(st$direction_pattern %in%
                    c("not_applicable", "single_trait", "concordant",
                      "antagonistic", "mixed")))
})

test_that("Anchor 5: predefined sets differ from all_traits", {
  dat <- mk_two_signal(seed = 5)
  ps <- data.frame(
    locus_id = "chr1:10000-10500",
    representative_snp = c("S1", "S2"),
    stringsAsFactors = FALSE
  )
  ps$conditioning_snps <- list("S2", "S1")
  pre_sets <- list("chr1:10000-10500::S1" = "Trait1",
                   "chr1:10000-10500::S2" = "Trait2")
  mk_pre <- function(mode, sets = NULL) {
    condped(
      dat$Y, G = dat$G, K = NULL,
      signal_mode = "predefined", candidate_mode = mode,
      control = list(null_control = list(maxit = 300L),
                     predefined_signals = ps,
                     predefined_candidate_sets = sets)
    )
  }
  f_pre <- mk_pre("predefined", pre_sets)
  f_all <- mk_pre("all_traits")
  for (s in c("chr1:10000-10500::S1", "chr1:10000-10500::S2")) {
    got_pre <- f_pre$candidate_traits$candidate_sets[[s]]
    got_all <- f_all$candidate_traits$candidate_sets[[s]]
    expect_identical(got_pre, pre_sets[[s]])
    expect_identical(got_all, paste0("Trait", 1:3))
    expect_false(identical(got_pre, got_all))
  }
})

test_that("Anchor 11: the global null is fitted exactly once", {
  dat <- mk_two_signal(seed = 6)
  calls <- 0L
  real_fit <- fit_mt_null
  local_mocked_bindings(
    fit_mt_null = function(...) {
      calls <<- calls + 1L
      real_fit(...)
    },
    .package = "CondPED"
  )
  run_cp(dat)
  expect_identical(calls, 1L)
})

test_that("Anchor 12: no effect_magnitude or Pattern classifier anywhere", {
  dat <- mk_two_signal(seed = 7)
  f <- run_cp(dat)
  flat <- unlist(lapply(f, function(x) {
    if (is.data.frame(x)) names(x) else if (is.list(x)) names(x) else NULL
  }))
  expect_false(any(grepl("effect_magnitude|Pattern|pattern_classifier",
                         flat)))
  expect_true(all(c("locus_id", "signal_id", "representative_snp",
                    "n_signals_in_locus", "effect_breadth",
                    "direction_pattern", "signal_status",
                    "numerical_status") %in% names(f$signal_summary)))
})

test_that("two-signal end-to-end anchor: no row mixing, no shared candidate set", {
  dat <- mk_two_signal(seed = 8)
  f <- run_cp(dat)
  sig_tab <- f$signals
  s1 <- sig_tab$signal_id[sig_tab$representative_snp == "S1"]
  s2 <- sig_tab$signal_id[sig_tab$representative_snp == "S2"]
  expect_identical(length(c(s1, s2)), 2L)
  tt <- f$candidate_traits$trait_table
  expect_true(all(tt$signal_id[tt$representative_snp == "S1"] == s1))
  expect_true(all(tt$signal_id[tt$representative_snp == "S2"] == s2))
  a1 <- f$candidate_traits$candidate_sets[[s1]]
  a2 <- f$candidate_traits$candidate_sets[[s2]]
  expect_false(is.null(a1))
  expect_false(is.null(a2))
  # S2's joint trait-3 effect is weak; S1's is strong -> sets differ
  expect_true("Trait3" %in% a1)
  expect_false("Trait3" %in% a2)
  expect_true(all(c("Trait1", "Trait2") %in% a2))
  # rho maps are computed per signal: different candidate-set sizes
  # give different maps, and shared keys carry different losses
  tab <- f$subset_analysis$subset_table
  m1 <- tab[tab$signal_id == s1, ]
  m2 <- tab[tab$signal_id == s2, ]
  expect_false(nrow(m1) == nrow(m2))   # k=3 vs k=2 maps
  shared <- intersect(m1$representing_key, m2$representing_key)
  l1 <- m1$representation_loss[match(shared, m1$representing_key)]
  l2 <- m2$representation_loss[match(shared, m2$representing_key)]
  expect_false(isTRUE(all.equal(unname(l1), unname(l2))))
})

test_that("predefined mode: oracle signals and candidate sets", {
  dat <- mk_two_signal(seed = 9)
  ps <- data.frame(
    locus_id = "chr1:10000-10500",
    representative_snp = c("S1", "S2"),
    stringsAsFactors = FALSE
  )
  ps$conditioning_snps <- list("S2", "S1")
  f <- condped(
    dat$Y, G = dat$G, K = NULL,
    signal_mode = "predefined",
    candidate_mode = "predefined",
    control = list(
      null_control = list(maxit = 300L),
      predefined_signals = ps,
      predefined_candidate_sets = list(
        "chr1:10000-10500::S1" = c("Trait1", "Trait2"),
        "chr1:10000-10500::S2" = c("Trait1", "Trait3")
      )
    )
  )
  expect_identical(f$status$code, "ok")
  expect_null(f$omnibus)
  expect_identical(nrow(f$signals), 2L)
  # predefined signal effects equal the joint-model blocks
  joint <- CondPED:::.fit_joint_signal_effects(f$null_fit, dat$G,
                                               c("S1", "S2"))
  expect_equal(f$effects$beta$beta[f$effects$beta$representative_snp == "S1"],
               unname(joint$beta["S1", ]), tolerance = 1e-10)
  expect_identical(f$candidate_traits$candidate_sets[[
    "chr1:10000-10500::S2"]], c("Trait1", "Trait3"))
})

test_that("predefined mode requires its inputs; resolve requires position", {
  dat <- mk_two_signal(seed = 10)
  expect_error(run_cp(dat, signal_mode = "predefined"),
               class = "condped_invalid_input")
  expect_error(
    condped(dat$Y, G = dat$G, K = NULL, chromosome = dat$chromosome),
    class = "condped_invalid_input"
  )
})

test_that("empty selection returns a complete empty structure", {
  set.seed(22)
  n <- 200; m <- 2
  G <- matrix(stats::rbinom(n * 100, 2, 0.3), n, 100)
  colnames(G) <- paste0("M", 1:100)
  Y <- matrix(stats::rnorm(n * m), n, m,
              dimnames = list(NULL, paste0("Trait", 1:m)))
  f <- condped(Y, G = G, K = NULL,
               position = seq_len(100) * 1000,
               alpha_omnibus = 1e-12, omnibus_adjust = "bonferroni",
               control = list(null_control = list(maxit = 200L)))
  expect_identical(f$status$code, "empty_selection")
  expect_true(f$status$ok)
  expect_null(f$signals)
  expect_identical(nrow(f$subset_analysis$subset_table), 0L)
})
