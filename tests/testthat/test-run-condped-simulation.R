# Tests for run_condped_simulation() (Stage 7): reproducible runner.

small_grid <- function() {
  data.frame(
    experiment = "signal_resolution",
    scenario = "single_multi_trait",
    analysis_mode = "full",
    comparison_pipeline = "condped_full",
    n = 150, m = 4, p = 120,
    locus_pve = 0.02, secondary_signal_pve = 0.01,
    local_ld = NA_character_, target_r2 = NA_real_,
    correlation = "block", target_loss = NA_real_,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}

test_that("1-3: canonical scenario_id is deterministic and order-invariant", {
  s1 <- list(experiment = "end_to_end", scenario = "x", n = 1000,
             m = 4L, tolerance = 0.1)
  s2 <- list(tolerance = 0.1, m = 4L, n = 1000,
             scenario = "x", experiment = "end_to_end")
  expect_identical(CondPED:::.canonical_scenario_id(s1),
                   CondPED:::.canonical_scenario_id(s2))
  expect_identical(CondPED:::.canonical_scenario_id(s1),
                   CondPED:::.canonical_scenario_id(s1))
  # environment fields never enter the id
  s3 <- c(s1, list(workers = 8L, out_dir = "/tmp/zzz"))
  expect_identical(CondPED:::.canonical_scenario_id(s1),
                   CondPED:::.canonical_scenario_id(s3))
  # a different setting changes the id
  s4 <- s1; s4$n <- 500
  expect_false(identical(CondPED:::.canonical_scenario_id(s1),
                         CondPED:::.canonical_scenario_id(s4)))
})

test_that("4: replicate seed depends only on master_seed/scenario_id/rep_id", {
  a <- CondPED:::.seed_for_rep(99L, "scen|x", 3L)
  b <- CondPED:::.seed_for_rep(99L, "scen|x", 3L)
  c <- CondPED:::.seed_for_rep(99L, "scen|x", 4L)
  d <- CondPED:::.seed_for_rep(100L, "scen|x", 3L)
  expect_identical(a, b)
  expect_false(identical(a, c))
  expect_false(identical(a, d))
})

test_that("5-6: workers 1 vs 2 and sequential vs parallel are identical", {
  dir1 <- tempfile(); dir.create(dir1)
  dir2 <- tempfile(); dir.create(dir2)
  g <- small_grid()
  r1 <- run_condped_simulation(g, reps = 2L, out_dir = dir1,
                               master_seed = 77L, workers = 1L,
                               progress = FALSE)
  r2 <- run_condped_simulation(g, reps = 2L, out_dir = dir2,
                               master_seed = 77L, workers = 2L,
                               backend = "parallel", progress = FALSE)
  for (i in seq_len(nrow(r1$manifest))) {
    o1 <- readRDS(r1$manifest$file[i])
    o2 <- readRDS(r2$manifest$file[i])
    expect_identical(o1$seed, o2$seed)
    expect_identical(o1$truth$beta, o2$truth$beta)
    expect_identical(o1$effect_estimates, o2$effect_estimates)
    expect_identical(o1$evaluation$overall_recovery,
                     o2$evaluation$overall_recovery)
  }
})

test_that("7-8: resume is invariant and saves are atomic", {
  dir1 <- tempfile(); dir.create(dir1)
  g <- small_grid()
  r1 <- run_condped_simulation(g, reps = 2L, out_dir = dir1,
                               master_seed = 88L, progress = FALSE)
  o1 <- readRDS(r1$manifest$file[1])
  r2 <- run_condped_simulation(g, reps = 2L, out_dir = dir1,
                               master_seed = 88L, progress = FALSE)
  expect_true(all(grepl("skipped", r2$manifest$status)))
  o2 <- readRDS(r2$manifest$file[1])
  expect_identical(o1$truth$beta, o2$truth$beta)
  # no temp files left behind
  leftovers <- list.files(dirname(r1$manifest$file[1]),
                          pattern = "[.]tmp$")
  expect_length(leftovers, 0L)
})

test_that("9-10: a failed replicate is captured and does not poison others", {
  dir1 <- tempfile(); dir.create(dir1)
  g <- small_grid()
  # corrupt one scenario so its simulation errors out
  g2 <- g
  g2$scenario <- "does_not_exist"
  g2$experiment <- "signal_resolution"
  both <- rbind(g, g2)
  res <- run_condped_simulation(both, reps = 1L, out_dir = dir1,
                                master_seed = 55L, progress = FALSE)
  expect_identical(res$completed, 1L)
  expect_identical(nrow(res$failed), 1L)
  bad <- res$manifest$file[res$manifest$status != "ok"]
  obj <- readRDS(bad)
  expect_false(obj$status$ok)
  expect_identical(obj$failure_stage, "runtime")
  expect_true(nchar(obj$error_message) > 0L)
  expect_true(all(c("scenario_id", "rep_id", "seed", "settings") %in%
                    names(obj)))
  # fail_policy = "stop" aborts
  expect_error(run_condped_simulation(g2, reps = 1L,
                                      out_dir = tempfile(),
                                      master_seed = 55L,
                                      fail_policy = "stop",
                                      progress = FALSE))
})

test_that("11: the evaluator reads saved replicates", {
  dir1 <- tempfile(); dir.create(dir1)
  res <- run_condped_simulation(small_grid(), reps = 2L,
                                out_dir = dir1, master_seed = 66L,
                                progress = FALSE)
  ev <- evaluate_condped_simulation(res$manifest$file,
                                    group_by = "architecture")
  expect_identical(ev$status$code, "ok")
  expect_identical(nrow(ev$replicate_metrics), 2L)
  expect_true(all(!is.na(ev$replicate_metrics$missed_signal_rate)))
})

paired_grid <- function() {
  g <- small_grid()
  g$experiment <- "end_to_end"
  g$scenario <- "mixed_multisignal"
  g$comparison_pipeline <- "paired"
  g$n <- 300; g$p <- 150
  g$locus_pve <- 0.08; g$secondary_signal_pve <- 0.04
  g
}

test_that("12-13: paired pipelines share one simulation identity", {
  dir1 <- tempfile(); dir.create(dir1)
  g <- paired_grid()
  res <- run_condped_simulation(g, reps = 1L, out_dir = dir1,
                                master_seed = 44L, progress = FALSE)
  obj <- readRDS(res$manifest$file[1])
  expect_identical(obj$comparison_pipeline, "paired")
  expect_true(all(c("lead_asset", "resolved_asset") %in%
                    names(obj$asset_results)))
  # one simulation: both pipelines saw the same truth/Y/G
  expect_identical(obj$truth$loci$locus_id,
                   obj$truth$signals$locus_id[1])
  lead_entry <- obj$asset_results$lead_asset$loci[[1]]
  res_entry <- obj$asset_results$resolved_asset$signals[[1]]
  # one simulation feeds both pipelines; ASSET status is well-formed
  expect_true(lead_entry$status %in%
                c("ok", "asset_not_available", "asset_failed",
                  "rank_deficient"))
  expect_true(!is.null(res_entry$status))
  # regression: full-mode candidate/subset fields must be populated
  # (flattened .pipeline_condped_full() output, signal-id named)
  expect_false(is.null(obj$estimates$candidate_sets))
  expect_true(length(obj$estimates$candidate_sets) > 0L)
  expect_true(all(nzchar(names(obj$estimates$candidate_sets))))
  expect_identical(names(obj$estimates$asset$condped$sets),
                   names(obj$estimates$candidate_sets))
  expect_false(is.null(obj$subset_table))
})

test_that("14: a real ASSET result is readable from the saved replicate", {
  testthat::skip_if_not_installed("ASSET")
  dir1 <- tempfile(); dir.create(dir1)
  g <- paired_grid()
  res <- run_condped_simulation(g, reps = 1L, out_dir = dir1,
                                master_seed = 45L, progress = FALSE)
  obj <- readRDS(res$manifest$file[1])
  st <- obj$asset_results$lead_asset$loci[[1]]$status
  expect_true(st %in% c("ok", "asset_failed", "rank_deficient"))
  if (identical(st, "ok")) {
    expect_true(is.finite(obj$asset_results$lead_asset$loci[[1]]$asset_p))
  }
})

test_that("15: the runner recomputes no statistics (source audit)", {
  # audit the deparsed bodies of the runner's functions; this works both
  # under devtools::test() (source tree) and R CMD check (installed pkg)
  ns <- asNamespace("CondPED")
  runner_fns <- c(
    "run_condped_simulation", ".pilot_grid", ".canonical_scenario_id",
    ".fnv1a_int", ".seed_for_rep", ".scenario_dir_key",
    ".run_one_replicate", ".execute_replicate",
    ".oracle_signals_from_truth", ".effects_from_predefined",
    ".asset_estimates_for_eval", ".save_replicate_atomic"
  )
  src <- paste(unlist(lapply(runner_fns, function(f) {
    deparse(get(f, envir = ns))
  })), collapse = "\n")
  # no formula re-implementation markers; genome-wide BH selection is
  # orchestration (Stage 6B-1), Holm must not be re-implemented here
  expect_false(grepl("crossprod\\(.*Vinv", src))
  expect_false(grepl("pchisq", src))
  expect_false(grepl("method\\s*=\\s*[\"']holm[\"']", src, ignore.case = TRUE))
  expect_false(grepl("effect_magnitude", src))
})
