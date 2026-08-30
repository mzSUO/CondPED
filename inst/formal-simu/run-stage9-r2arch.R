## Stage 9 batch 10: R2 partially_representable architecture
## (freeze §1.14), R=200, signal_trait_oracle.
## Launch: bash run_sim.sh stage9-r2arch inst/formal-simu/run-stage9-r2arch.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-r2arch"
grid <- data.frame(
  experiment = "trait_representation",
  scenario = "partially_representable",
  analysis_mode = "signal_trait_oracle",
  comparison_pipeline = "condped_full",
  n = 1000L, m = 4L, p = 1000L,
  locus_pve = 0.02, secondary_signal_pve = 0.02,
  local_ld = NA_character_, target_r2 = NA_real_,
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
             rep_exact = ev$rep_exact_family_recovery,
             rep_card = ev$rep_min_cardinality_recovery,
             rho_mae = ev$representation_loss_mae,
             threshold_acc = ev$representation_threshold_accuracy,
             irr_recovery = ev$irr_family_recovery,
             cand_exact = ev$candidate_exact_recovery,
             stringsAsFactors = FALSE)
}))
summ <- data.frame(
  n = nrow(rows),
  rep_exact = stage9_prop(rows$rep_exact)["estimate"],
  rep_card = stage9_prop(rows$rep_card)["estimate"],
  rho_mae = stage9_cont(rows$rho_mae)["estimate"],
  threshold_acc = stage9_prop(rows$threshold_acc)["estimate"],
  irr_recovery = stage9_prop(rows$irr_recovery)["estimate"],
  cand_exact = stage9_prop(rows$cand_exact)["estimate"]
)
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(rows, file.path(out_dir, "per_rep.csv"), row.names = FALSE)

lines <- c(
  "# Stage 9 batch 10: R2 partially_representable (freeze §1.14)",
  "",
  sprintf("- R=200, %d/%d ok", sum(final$status == "ok"), nrow(final)),
  "",
  sprintf("- rep_exact %.3f; rep_card %.3f; rho MAE %.4f; thr-side acc %.3f; Irr %.3f; cand exact %.3f",
          summ$rep_exact, summ$rep_card, summ$rho_mae, summ$threshold_acc,
          summ$irr_recovery, summ$cand_exact),
  "",
  "- 定位（§1.14/§1.32）：R2 放 supplementary，作为 R1/R3 主对比的中间架构"
)
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
cat("BATCH10 DONE\n")
