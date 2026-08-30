## Stage 9 batch 9: Sigma_P = I4 phenotype-correlation sensitivity
## (freeze §1.58). Main-scenario rerun under correlation = "independent":
## two_linked_trait_specific, R=200, paired. Distinguishes LD-induced
## mixing from trait-correlation-induced pseudo-mixing.
## Launch: bash run_sim.sh stage9-i4 inst/formal-simu/run-stage9-i4.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-i4"
grid <- data.frame(
  experiment = "signal_resolution",
  scenario = "two_linked_trait_specific",
  analysis_mode = "full", comparison_pipeline = "paired",
  n = 1000L, m = 4L, p = 1000L,
  locus_pve = 0.02, secondary_signal_pve = 0.02,
  local_ld = NA_character_, target_r2 = 0.3,
  correlation = "independent", target_loss = NA_real_,
  tolerance = 0.10, alpha_omnibus = 0.05,
  stringsAsFactors = FALSE
)

final <- stage9_run_grid(grid, out_dir, reps_total = 200L, workers = 4L)

rows <- do.call(rbind, lapply(final$file[final$status == "ok"], function(f) {
  o <- readRDS(f)
  if (!is.data.frame(o$evaluation)) return(NULL)
  ev <- o$evaluation
  data.frame(rep_id = o$rep_id,
             locus_detection = ev$locus_detection_rate,
             secondary_power = ev$secondary_signal_power,
             overall_recovery = ev$overall_recovery,
             cand_exact = ev$candidate_exact_recovery,
             direction = ev$direction_recovery,
             count_exact = ev$signal_count_exact_recovery,
             extra_rate = ev$extra_signal_rate,
             pseudo_rate = ev$pseudo_multitrait_interpretation_rate,
             stringsAsFactors = FALSE)
}))
summ <- data.frame(
  n = nrow(rows),
  locus_detection = stage9_prop(rows$locus_detection)["estimate"],
  secondary_power = stage9_prop(rows$secondary_power)["estimate"],
  overall_recovery = stage9_prop(rows$overall_recovery)["estimate"],
  cand_exact = stage9_prop(rows$cand_exact)["estimate"],
  direction = stage9_prop(rows$direction)["estimate"],
  count_exact = stage9_prop(rows$count_exact)["estimate"],
  extra_rate = stage9_prop(rows$extra_rate)["estimate"],
  pseudo_rate = stage9_cont(rows$pseudo_rate)["estimate"]
)
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(rows, file.path(out_dir, "per_rep.csv"), row.names = FALSE)

## strict matching-based recovery (Stage 8.2 convention; headline metric)
strict <- stage9_strict_recovery(final$file[final$status == "ok"],
                                 workers = 4L)
write.csv(strict, file.path(out_dir, "strict_recovery_per_rep.csv"),
          row.names = FALSE)
strict_rate <- mean(strict$full_recovery)
strict_sec <- mean(strict$secondary_recovery)

lines <- c(
  "# Stage 9 batch 9: Sigma_P = I4 sensitivity (freeze §1.58)",
  "",
  sprintf("- two_linked under correlation=independent, R=200, %d/%d ok",
          sum(final$status == "ok"), nrow(final)),
  "",
  sprintf("- locus detection %.3f; secondary power %.3f; overall recovery %.3f; cand exact %.3f; direction %.3f",
          summ$locus_detection, summ$secondary_power, summ$overall_recovery,
          summ$cand_exact, summ$direction),
  sprintf("- count_exact %.3f; extra rate %.3f; pseudo-multitrait rate %.3f",
          summ$count_exact, summ$extra_rate, summ$pseudo_rate),
  "",
  "- 对照主场景（block cor, Stage 8.2）：pseudo rate 应在 I4 下明显下降；",
  "  若不降，说明 pseudo-mixing 有 trait-correlation 之外的来源"
)
lines <- c(lines, "",
  sprintf("- 严格匹配口径（Stage 8.2 约定）：full recovery %.3f, secondary %.3f",
          strict_rate, strict_sec))
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
cat("BATCH9 DONE\n")
