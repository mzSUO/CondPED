## Stage 9 batch 4: representation strength rho sensitivity (freeze §1.59.1).
## R1 (highly_representable) at target_loss in {0.02, 0.08}, R=200/cell,
## signal_trait_oracle. §1.59.1: verify the R1/R3 conclusion does not depend
## on the single rho_main = 0.05 point.
## Launch: bash run_sim.sh stage9-rho inst/formal-simu/run-stage9-rho.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-rho"
grid <- do.call(rbind, lapply(c(0.02, 0.08), function(rho) {
  data.frame(
    experiment = "trait_representation",
    scenario = "highly_representable",
    analysis_mode = "signal_trait_oracle",
    comparison_pipeline = "condped_full",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = NA_real_,
    correlation = "block", target_loss = rho,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}))

final <- stage9_run_grid(grid, out_dir, reps_total = 200L, workers = 4L)

grid$key <- vapply(seq_len(nrow(grid)), function(i) stage9_key(grid[i, ]),
                   character(1))
final$dir <- basename(dirname(final$file))
final$rho <- grid$target_loss[match(final$dir, grid$key)]

rows <- do.call(rbind, lapply(final$file[final$status == "ok"], function(f) {
  o <- readRDS(f)
  if (!is.data.frame(o$evaluation)) return(NULL)
  ev <- o$evaluation
  data.frame(rho = o$settings$target_loss,
             rep_exact = ev$rep_exact_family_recovery,
             rep_card = ev$rep_min_cardinality_recovery,
             rho_mae = ev$representation_loss_mae,
             threshold_acc = ev$representation_threshold_accuracy,
             cand_exact = ev$candidate_exact_recovery,
             stringsAsFactors = FALSE)
}))
summ <- do.call(rbind, lapply(split(rows, rows$rho), function(d) {
  data.frame(
    rho = d$rho[1], n = nrow(d),
    rep_exact = stage9_prop(d$rep_exact)["estimate"],
    rep_card = stage9_prop(d$rep_card)["estimate"],
    rho_mae = stage9_cont(d$rho_mae)["estimate"],
    threshold_acc = stage9_prop(d$threshold_acc)["estimate"],
    cand_exact = stage9_prop(d$cand_exact)["estimate"],
    stringsAsFactors = FALSE)
}))
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(rows, file.path(out_dir, "per_rep.csv"), row.names = FALSE)

lines <- c(
  "# Stage 9 batch 4: rho sensitivity (freeze §1.59.1)",
  "",
  sprintf("- R1 target_loss in {0.02, 0.08}, R=200/cell, %d/%d ok",
          sum(final$status == "ok"), nrow(final)),
  "",
  "| rho | rep_exact | rep_card | rho MAE | thr-side acc | cand exact |",
  "|---|---|---|---|---|---|",
  apply(summ, 1, function(r) sprintf("| %.2f | %.3f | %.3f | %.4f | %.3f | %.3f |",
    as.numeric(r[["rho"]]), as.numeric(r[["rep_exact"]]),
    as.numeric(r[["rep_card"]]), as.numeric(r[["rho_mae"]]),
    as.numeric(r[["threshold_acc"]]), as.numeric(r[["cand_exact"]]))),
  "",
  "- reference: main rho=0.05 (Stage 8.3, R=500): rep_exact 0.722, rho MAE 0.035",
  "- 核对点：R1/R3 结论方向不随 rho 点翻转"
)
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
cat("BATCH4 DONE\n")
