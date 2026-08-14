## Stage 7 pilot summary: read the 20-replicate pilot output and summarise
## per scenario family with the frozen evaluator. No statistics are
## recomputed here; everything comes from saved replicates +
## evaluate_condped_simulation().
devtools::load_all(quiet = TRUE)

pilot_dir <- "inst/validation/output/pilot"
manifest <- readRDS(file.path(pilot_dir, "manifest.rds"))

# scenario is embedded in the canonical scenario_id ("k=v|k=v" string)
manifest$scenario <- sub(".*(?:^|\\|)scenario=([^|]*).*", "\\1",
                         manifest$scenario_id)
# resumed runs mark pre-existing complete replicates as "skipped_ok";
# both statuses denote a usable, complete replicate file
manifest$run_status <- ifelse(manifest$status == "skipped_ok", "ok",
                              manifest$status)

cat("=== replicate counts ===\n")
print(table(manifest$scenario, manifest$run_status))
cat("total files:", nrow(manifest), "\n")

cat("\n=== failure rates per scenario ===\n")
fr <- tapply(manifest$run_status != "ok", manifest$scenario, mean)
print(round(fr, 3))

cat("\n=== runtime / file size ===\n")
runtimes <- manifest$runtime
sizes <- file.info(manifest$file)$size
cat("runtime per replicate (s): median", round(median(runtimes, na.rm = TRUE), 1),
    " max", round(max(runtimes, na.rm = TRUE), 1), "\n")
cat("file size (MB): median", round(median(sizes) / 1e6, 2),
    " max", round(max(sizes) / 1e6, 2),
    " total", round(sum(sizes) / 1e6, 1), "\n")

summarise_family <- function(label, scenarios, metrics) {
  cat("\n=== ", label, " ===\n", sep = "")
  files <- manifest$file[manifest$scenario %in% scenarios &
                           manifest$run_status == "ok"]
  cat("replicates included:", length(files), "\n")
  ev <- evaluate_condped_simulation(files, group_by = "architecture")
  sm <- ev$summary
  keep <- intersect(metrics, names(sm))
  print(sm[, c("architecture", "n", keep), drop = FALSE],
        row.names = FALSE)
  invisible(ev)
}

sim1 <- summarise_family(
  "Simulation I / signal_resolution (mode=full)",
  c("null", "single_multi_trait", "two_linked_trait_specific",
    "two_heterogeneous"),
  c("type1_omnibus", "power_omnibus", "locus_detection_rate",
    "signal_count_exact_recovery", "secondary_signal_power",
    "missed_signal_rate", "extra_signal_rate"))

sim2a <- summarise_family(
  "Simulation II-A / trait_representation (mode=signal_oracle)",
  c("trait_specific", "two_trait_concordant", "two_trait_antagonistic",
    "broad_concordant"),
  c("beta_rmse", "beta_coverage", "candidate_tpr", "candidate_fdp",
    "candidate_exact_recovery", "effect_breadth_bias", "direction_recovery"))

sim2b <- summarise_family(
  "Simulation II-B / representation architectures (mode=signal_trait_oracle)",
  c("highly_representable", "partially_representable",
    "strongly_nonredundant"),
  c("representation_loss_bias", "representation_loss_rmse",
    "representation_threshold_accuracy", "rep_min_cardinality_recovery",
    "rep_exact_family_recovery", "rep_tie_coverage",
    "irr_family_recovery", "irr_missed_rate", "irr_extra_split_rate"))

sim3 <- summarise_family(
  "Simulation III / end_to_end paired comparators (mode=full)",
  c("single_highly_representable", "single_nonredundant",
    "linked_pseudo_multitrait", "mixed_multisignal"),
  c("lead_trait_breadth_inflation", "resolved_trait_breadth_inflation",
    "pseudo_multitrait_interpretation_rate", "asset_trait_set_jaccard",
    "asset_trait_set_precision", "asset_trait_set_recall",
    "asset_direction_recovery", "rep_exact_family_recovery"))

## ASSET backend status audit across saved replicates
cat("\n=== ASSET backend status (paired pipelines) ===\n")
is_paired <- grepl("comparison_pipeline=paired", manifest$scenario_id,
                   fixed = TRUE)
asset_status <- unlist(lapply(manifest$file[
  is_paired & manifest$run_status == "ok"],
  function(f) {
    o <- readRDS(f)
    c(vapply(o$asset_results$lead_asset$loci, `[[`, "", "status"),
      vapply(o$asset_results$resolved_asset$signals, `[[`, "", "status"))
  }))
print(table(asset_status))

cat("\nPILOT SUMMARY DONE\n")
