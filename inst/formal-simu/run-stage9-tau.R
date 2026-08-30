## Stage 9 batch 5: representation threshold tau sensitivity (freeze §1.59).
## R1 (highly_representable, rho_main=0.05) at tolerance in {0.05, 0.20},
## R=200/cell, signal_trait_oracle.
## Launch: bash run_sim.sh stage9-tau inst/formal-simu/run-stage9-tau.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-tau"
## NOTE (evidence recorded): the tau = 0.05 cell is structurally
## infeasible — R1's acceptance rule requires the truth singleton rho to
## lie strictly below tau, so rho_main = 0.05 with tau = 0.05 collapses to
## the boundary and no generation is accepted (verified: first run's
## tau=0.05 cell failed 200/200 with "invalid_covariance"/unstable; the
## failed replicate files are kept in output/stage9-tau/ as evidence).
## This is itself the tau-sensitivity finding: the R1 architecture only
## exists for tau > rho. The batch therefore runs tau = 0.20 only.
grid <- do.call(rbind, lapply(c(0.20), function(tau) {
  data.frame(
    experiment = "trait_representation",
    scenario = "highly_representable",
    analysis_mode = "signal_trait_oracle",
    comparison_pipeline = "condped_full",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = NA_real_,
    correlation = "block", target_loss = 0.05,
    tolerance = tau, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}))

final <- stage9_run_grid(grid, out_dir, reps_total = 200L, workers = 4L)

grid$key <- vapply(seq_len(nrow(grid)), function(i) stage9_key(grid[i, ]),
                   character(1))
final$dir <- basename(dirname(final$file))
final$tau <- grid$tolerance[match(final$dir, grid$key)]

rows <- do.call(rbind, lapply(final$file[final$status == "ok"], function(f) {
  o <- readRDS(f)
  if (!is.data.frame(o$evaluation)) return(NULL)
  ev <- o$evaluation
  data.frame(tau = o$settings$tolerance,
             rep_exact = ev$rep_exact_family_recovery,
             rep_card = ev$rep_min_cardinality_recovery,
             rho_mae = ev$representation_loss_mae,
             threshold_acc = ev$representation_threshold_accuracy,
             stringsAsFactors = FALSE)
}))
summ <- do.call(rbind, lapply(split(rows, rows$tau), function(d) {
  data.frame(
    tau = d$tau[1], n = nrow(d),
    rep_exact = stage9_prop(d$rep_exact)["estimate"],
    rep_card = stage9_prop(d$rep_card)["estimate"],
    rho_mae = stage9_cont(d$rho_mae)["estimate"],
    threshold_acc = stage9_prop(d$threshold_acc)["estimate"],
    stringsAsFactors = FALSE)
}))
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(rows, file.path(out_dir, "per_rep.csv"), row.names = FALSE)

lines <- c(
  "# Stage 9 batch 5: tau sensitivity (freeze §1.59)",
  "",
  sprintf("- R1 (rho=0.05) at tau in {0.05, 0.20}, R=200/cell, %d/%d ok",
          sum(final$status == "ok"), nrow(final)),
  "",
  "| tau | rep_exact | rep_card | rho MAE | thr-side acc |",
  "|---|---|---|---|---|",
  apply(summ, 1, function(r) sprintf("| %.2f | %.3f | %.3f | %.4f | %.3f |",
    as.numeric(r[["tau"]]), as.numeric(r[["rep_exact"]]),
    as.numeric(r[["rep_card"]]), as.numeric(r[["rho_mae"]]),
    as.numeric(r[["threshold_acc"]]))),
  "",
  "- reference: main tau=0.10 (Stage 8.3, R=500): rep_exact 0.722",
  "- 核对点：tau=0.20 更宽 → rep_exact 应升（边界效应方向性）",
  "- **tau=0.05 格结构性不可行**：R1 要求 truth rho < tau，rho_main=0.05 与",
  "  tau=0.05 重合于边界，acceptance rule 全部拒绝（200/200 failed，",
  "  本目录保留失败 rep 为证据）——tau 敏感性在 tau <= rho_main 处",
  "  不是性能下降而是架构不存在"
)
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
cat("BATCH5 DONE\n")
