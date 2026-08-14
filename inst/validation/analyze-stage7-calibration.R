## Stage 7.1 Phase B analysis: build the full calibration table (section D)
## from the completed calibration grid run. No parameter set is dropped;
## selection uses only the pre-defined informative-window criteria.
devtools::load_all(quiet = TRUE)

out_dir <- "inst/validation/output/calibration"
grid <- readRDS(file.path(out_dir, "calibration_grid.rds"))
manifest <- readRDS(file.path(out_dir, "manifest.rds"))

sid_for <- function(row) {
  settings <- as.list(row[setdiff(names(row), "parameter_set_id")])
  CondPED:::.canonical_scenario_id(settings)
}
grid$scenario_id <- vapply(seq_len(nrow(grid)), function(i) {
  sid_for(grid[i, , drop = FALSE])
}, character(1))

metrics_needed <- c(
  "locus_detection_rate", "conditional_signal_count_exact_recovery",
  "conditional_secondary_signal_power", "extra_signal_rate",
  "candidate_exact_recovery", "rep_exact_family_recovery",
  "lead_trait_breadth_inflation", "resolved_trait_breadth_inflation",
  "signal_resolution_attempt_rate"
)

rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  files <- manifest$file[manifest$scenario_id == g$scenario_id]
  man <- manifest[manifest$scenario_id == g$scenario_id, ]
  fail_rate <- mean(man$status != "ok")
  ok_files <- files[man$status == "ok"]
  ev <- evaluate_condped_simulation(ok_files, metrics = metrics_needed,
                                    group_by = character())
  s <- ev$summary
  r <- ev$replicate_metrics
  getm <- function(x) if (x %in% names(s)) s[[x]][1] else NA_real_
  # per-replicate conditional chain + paired breadth improvement
  ce <- r$candidate_exact_recovery == 1
  rep_given_cand <- if (any(ce, na.rm = TRUE)) {
    mean(r$rep_exact_family_recovery[ce], na.rm = TRUE)
  } else {
    NA_real_
  }
  improvement <- mean(abs(r$lead_trait_breadth_inflation) -
                        abs(r$resolved_trait_breadth_inflation),
                      na.rm = TRUE)
  rows[[i]] <- data.frame(
    scenario = g$scenario,
    parameter_set_id = g$parameter_set_id,
    locus_pve = g$locus_pve,
    secondary_signal_pve = g$secondary_signal_pve,
    target_r2 = g$target_r2,
    n_reps = nrow(man),
    locus_detection = getm("locus_detection_rate"),
    conditional_signal_exact =
      getm("conditional_signal_count_exact_recovery"),
    conditional_secondary_power =
      getm("conditional_secondary_signal_power"),
    extra_signal_rate = getm("extra_signal_rate"),
    candidate_exact = getm("candidate_exact_recovery"),
    rep_exact = getm("rep_exact_family_recovery"),
    rep_exact_given_candidate_exact = rep_given_cand,
    lead_breadth_error = getm("lead_trait_breadth_inflation"),
    resolved_breadth_error = getm("resolved_trait_breadth_inflation"),
    resolution_improvement = improvement,
    runtime_median_s = stats::median(man$runtime),
    failure_rate = fail_rate,
    stringsAsFactors = FALSE
  )
  cat(sprintf("[%d/%d] %s done\n", i, nrow(grid), g$parameter_set_id))
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(out_dir, "calibration_table.csv"),
          row.names = FALSE)
print(tab, row.names = FALSE)
cat("\nCALIBRATION TABLE DONE\n")
