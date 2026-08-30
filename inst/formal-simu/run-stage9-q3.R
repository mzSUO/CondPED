## Stage 9 batch 7: K_Q = 3 scalability (freeze §1.61).
## three_linked_trait_specific (registry addition, Stage 7.7), R=200,
## paired pipelines. Only computational stability and statistical validity
## are claimed (per §1.61).
## Launch: bash run_sim.sh stage9-q3 inst/formal-simu/run-stage9-q3.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-q3"
grid <- data.frame(
  experiment = "signal_resolution",
  scenario = "three_linked_trait_specific",
  analysis_mode = "full", comparison_pipeline = "paired",
  n = 1000L, m = 4L, p = 1000L,
  locus_pve = 0.02, secondary_signal_pve = 0.02,
  local_ld = NA_character_, target_r2 = 0.3,
  correlation = "block", target_loss = NA_real_,
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
             overall_recovery = ev$overall_recovery,
             count_exact = ev$signal_count_exact_recovery,
             extra_rate = ev$extra_signal_rate,
             runtime = o$runtime,
             n_signals = o$diagnostics$n_signals,
             stringsAsFactors = FALSE)
}))
summ <- data.frame(
  n = nrow(rows),
  locus_detection = stage9_prop(rows$locus_detection)["estimate"],
  overall_recovery = stage9_prop(rows$overall_recovery)["estimate"],
  count_exact = stage9_prop(rows$count_exact)["estimate"],
  extra_rate = stage9_prop(rows$extra_rate)["estimate"],
  mean_runtime = mean(rows$runtime, na.rm = TRUE),
  mean_n_signals = mean(rows$n_signals, na.rm = TRUE)
)
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(rows, file.path(out_dir, "per_rep.csv"), row.names = FALSE)

lines <- c(
  "# Stage 9 batch 7: K_Q = 3 scalability (freeze §1.61)",
  "",
  sprintf("- three_linked_trait_specific, R=200, %d/%d ok", sum(final$status == "ok"), nrow(final)),
  "",
  sprintf("- locus detection %.3f; overall recovery %.3f; count_exact %.3f; extra rate %.3f",
          summ$locus_detection, summ$overall_recovery, summ$count_exact, summ$extra_rate),
  sprintf("- mean signals estimated %.2f; mean rep runtime %.1fs", summ$mean_n_signals, summ$mean_runtime),
  "",
  "- 定位（§1.61）：K_Q=3 只声明确算稳定与统计有效，不进主文",
  "- 0 numerical failure = 计算稳定；locus detection 与 recovery 为统计有效性的 sanity 证据"
)
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
cat("BATCH7 DONE\n")
