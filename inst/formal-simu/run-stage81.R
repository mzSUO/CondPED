## Stage 8.1 pilot: all main scenarios, R = 20 replicates each.
## Parameters follow the freeze document
## `inst/formal-simu/05 模拟20260825.md` (secondary PVE 0.020, r2 = 0.3,
## rho = 0.05 for R1/E1, n = 1000, m = 4, p = 1000, h2 = 0.5, block
## correlation, tolerance = 0.10, alpha_omnibus = 0.05).
##
## Pairing (freeze 1.64): within every replicate, Lead+ASSET /
## Resolved+ASSET / CondPED share the same Y / G / truth / seed — the
## formal runner (`run_condped_simulation`, comparison_pipeline =
## "paired") enforces this by construction; per-replicate seeds depend
## only on (master_seed, canonical scenario_id, rep_id).
##
## Checkpoints: one RDS per replicate under
## inst/formal-simu/output/stage81/<scenario_key>/rep_NNNN.rds (atomic
## save by the runner). Failed replicates are retried exactly once and
## recorded in retry_log.csv. Resume: existing ok replicates are skipped.
##
## Launch ONLY via:  bash run_sim.sh stage81 inst/formal-simu/run-stage81.R

devtools::load_all(quiet = TRUE)

master_seed <- 20260825L
reps <- 20L
workers <- as.integer(Sys.getenv("STAGE81_WORKERS", "4"))
out_dir <- "inst/formal-simu/output/stage81"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- grid: freeze main scenarios ------------------------------------------
mk <- function(experiment, scenario, analysis_mode, comparison_pipeline,
               locus_pve = 0.02, secondary_signal_pve = 0.02,
               target_r2 = NA_real_, target_loss = NA_real_) {
  data.frame(
    experiment = experiment, scenario = scenario,
    analysis_mode = analysis_mode, comparison_pipeline = comparison_pipeline,
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = locus_pve, secondary_signal_pve = secondary_signal_pve,
    local_ld = NA_character_, target_r2 = target_r2,
    correlation = "block", target_loss = target_loss,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}
grid <- rbind(
  ## Simulation I
  mk("signal_resolution", "null", "full", "condped_full",
     locus_pve = 0.02),                                   # I-0
  mk("signal_resolution", "single_multi_trait", "full", "condped_full"),  # I-1
  mk("signal_resolution", "two_linked_trait_specific", "full",
     "condped_full", target_r2 = 0.3),                    # I-2
  ## Simulation II-A (signal_oracle)
  mk("trait_representation", "trait_specific", "signal_oracle",
     "condped_full"),                                     # A1
  mk("trait_representation", "two_trait_concordant", "signal_oracle",
     "condped_full"),                                     # A2
  mk("trait_representation", "two_trait_antagonistic", "signal_oracle",
     "condped_full"),                                     # A3
  mk("trait_representation", "broad_concordant", "signal_oracle",
     "condped_full"),                                     # A4
  ## Simulation II-B (signal_trait_oracle), rho = 0.05 for R1
  mk("trait_representation", "highly_representable", "signal_trait_oracle",
     "condped_full", target_loss = 0.05),                 # R1
  mk("trait_representation", "strongly_nonredundant", "signal_trait_oracle",
     "condped_full"),                                     # R3
  ## Simulation III (paired three pipelines)
  mk("end_to_end", "single_highly_representable", "full", "paired",
     target_loss = 0.05),                                 # E1
  mk("end_to_end", "single_nonredundant", "full", "paired"),  # E2
  mk("end_to_end", "linked_pseudo_multitrait", "full", "paired",
     target_r2 = 0.3)                                     # E3
)
write.csv(grid, file.path(out_dir, "grid.csv"), row.names = FALSE)

## ---- peak-RSS sampler (master + direct children workers) ------------------
rss_file <- file.path(out_dir, "peak_rss_kb.txt")
sampler <- parallel::mcparallel({
  master <- Sys.getpid()
  peak <- 0
  repeat {
    rss <- tryCatch({
      pids <- c(master, as.integer(system2("pgrep", c("-P", master),
                                           stdout = TRUE)))
      vals <- as.integer(system2("ps", c("-o", "rss=", "-p",
                                         paste(pids, collapse = ",")),
                                 stdout = TRUE))
      sum(vals)
    }, error = function(e) 0)
    if (is.finite(rss) && rss > peak) {
      peak <- rss
      writeLines(format(peak, scientific = FALSE), rss_file)
    }
    Sys.sleep(5)
  }
}, silent = TRUE)
on.exit(parallel::mccollect(sampler), add = TRUE)

## ---- main run --------------------------------------------------------------
t0 <- Sys.time()
res <- run_condped_simulation(
  grid, reps = reps, out_dir = out_dir,
  master_seed = master_seed, workers = workers, backend = "parallel",
  resume = TRUE, fail_policy = "save_unstable", progress = FALSE
)
wall_main <- difftime(Sys.time(), t0, units = "secs")

## ---- one automatic retry of failed replicates ------------------------------
manifest <- res$manifest
failed <- manifest[manifest$status != "ok", , drop = FALSE]
retry_log <- data.frame()
if (nrow(failed) > 0L) {
  retry_log <- data.frame(
    scenario_id = failed$scenario_id, rep_id = failed$rep_id,
    first_status = failed$status, stringsAsFactors = FALSE
  )
  unlink(failed$file)          # force one clean retry with identical seed
  res2 <- run_condped_simulation(
    grid, reps = reps, out_dir = out_dir,
    master_seed = master_seed, workers = workers, backend = "parallel",
    resume = TRUE, fail_policy = "save_unstable", progress = FALSE
  )
  manifest <- res2$manifest
  retry_log$retry_status <- manifest$status[
    match(paste(retry_log$scenario_id, retry_log$rep_id),
          paste(manifest$scenario_id, manifest$rep_id))]
}
write.csv(retry_log, file.path(out_dir, "retry_log.csv"), row.names = FALSE)

## ---- missing-output handling check (freeze 1.62 item 9) --------------------
## (a) resume must skip all existing ok replicates;
## (b) a deleted ok replicate must be regenerated with identical content.
t_resume <- Sys.time()
res3 <- run_condped_simulation(
  grid, reps = reps, out_dir = out_dir,
  master_seed = master_seed, workers = workers, backend = "parallel",
  resume = TRUE, fail_policy = "save_unstable", progress = FALSE
)
n_skipped_ok <- sum(res3$manifest$status == "skipped_ok")
ok_files <- manifest$file[manifest$status == "ok"]
probe <- ok_files[1L]
probe_before <- readRDS(probe)
unlink(probe)
res4 <- run_condped_simulation(
  grid, reps = reps, out_dir = out_dir,
  master_seed = master_seed, workers = 1L, backend = "sequential",
  resume = TRUE, fail_policy = "save_unstable", progress = FALSE
)
probe_after <- readRDS(probe)
probe_identical <- identical(probe_before$evaluation, probe_after$evaluation) &&
  identical(probe_before$seed, probe_after$seed)
missing_check <- data.frame(
  n_ok_before = length(ok_files),
  n_skipped_ok_on_resume = n_skipped_ok,
  probe_file = probe,
  probe_regenerated = file.exists(probe),
  probe_content_identical = probe_identical
)
write.csv(missing_check, file.path(out_dir, "missing_output_check.csv"),
          row.names = FALSE)

## ---- final manifest from disk (ground truth) -------------------------------
## resume runs report "skipped_ok" for existing replicates; rebuild the
## authoritative status table from the saved per-replicate RDS files.
rep_files <- list.files(out_dir, pattern = "^rep_[0-9]+\\.rds$",
                        recursive = TRUE, full.names = TRUE)
disk_status <- vapply(rep_files, function(f) {
  o <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(o)) return("corrupt")
  if (isTRUE(o$status$ok)) "ok" else paste0("failed:", o$status$code)
}, character(1))
final <- data.frame(file = rep_files, status = disk_status,
                    stringsAsFactors = FALSE)
saveRDS(final, file.path(out_dir, "manifest.rds"))

wall_total <- difftime(Sys.time(), t0, units = "secs")
write.csv(data.frame(metric = c("wall_main_seconds", "wall_total_seconds",
                                "n_replicates", "n_ok", "n_failed",
                                "n_retried", "n_retry_recovered"),
                     value = c(wall_main, wall_total, nrow(final),
                               sum(final$status == "ok"),
                               sum(final$status != "ok"),
                               nrow(retry_log),
                               sum(retry_log$retry_status %in% "ok"))),
          file.path(out_dir, "runtime.csv"), row.names = FALSE)
cat("STAGE81 PILOT RUN DONE  wall:", round(wall_total, 1), "s  ok:",
    sum(final$status == "ok"), "/", nrow(final), "\n")
