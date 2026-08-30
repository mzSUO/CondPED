## Stage 9 batch 8: I-3 two_heterogeneous + E4 mixed_multisignal
## (freeze §1.24 / §1.46), R=200 each, paired pipelines.
## Launch: bash run_sim.sh stage9-i3e4 inst/formal-simu/run-stage9-i3e4.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-i3e4"
grid <- rbind(
  data.frame(
    experiment = "signal_resolution", scenario = "two_heterogeneous",
    analysis_mode = "full", comparison_pipeline = "paired",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = 0.3,
    correlation = "block", target_loss = NA_real_,
    tolerance = 0.10, alpha_omnibus = 0.05, stringsAsFactors = FALSE),
  data.frame(
    experiment = "end_to_end", scenario = "mixed_multisignal",
    analysis_mode = "full", comparison_pipeline = "paired",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = 0.3,
    correlation = "block", target_loss = NA_real_,
    tolerance = 0.10, alpha_omnibus = 0.05, stringsAsFactors = FALSE)
)

final <- stage9_run_grid(grid, out_dir, reps_total = 200L, workers = 4L)

grid$key <- vapply(seq_len(nrow(grid)), function(i) stage9_key(grid[i, ]),
                   character(1))
final$dir <- basename(dirname(final$file))
final$scenario <- grid$scenario[match(final$dir, grid$key)]

rows <- do.call(rbind, lapply(final$file[final$status == "ok"], function(f) {
  o <- readRDS(f)
  if (!is.data.frame(o$evaluation)) return(NULL)
  ev <- o$evaluation
  data.frame(scenario = o$settings$architecture,
             locus_detection = ev$locus_detection_rate,
             overall_recovery = ev$overall_recovery,
             secondary_power = ev$secondary_signal_power,
             cand_exact = ev$candidate_exact_recovery,
             direction = ev$direction_recovery,
             count_exact = ev$signal_count_exact_recovery,
             extra_rate = ev$extra_signal_rate,
             stringsAsFactors = FALSE)
}))
summ <- do.call(rbind, lapply(split(rows, rows$scenario), function(d) {
  data.frame(
    scenario = d$scenario[1], n = nrow(d),
    locus_detection = stage9_prop(d$locus_detection)["estimate"],
    secondary_power = stage9_prop(d$secondary_power)["estimate"],
    overall_recovery = stage9_prop(d$overall_recovery)["estimate"],
    cand_exact = stage9_prop(d$cand_exact)["estimate"],
    direction = stage9_prop(d$direction)["estimate"],
    count_exact = stage9_prop(d$count_exact)["estimate"],
    extra_rate = stage9_prop(d$extra_rate)["estimate"],
    stringsAsFactors = FALSE)
}))
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(rows, file.path(out_dir, "per_rep.csv"), row.names = FALSE)

lines <- c(
  "# Stage 9 batch 8: I-3 two_heterogeneous + E4 mixed_multisignal (freeze §1.24 / §1.46)",
  "",
  sprintf("- R=200 each, paired, %d/%d ok", sum(final$status == "ok"), nrow(final)),
  "",
  "| scenario | locus det | secondary power | overall recovery | cand exact | direction | count_exact | extra rate |",
  "|---|---|---|---|---|---|---|---|",
  apply(summ, 1, function(r) sprintf("| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |",
    r[["scenario"]], as.numeric(r[["locus_detection"]]),
    as.numeric(r[["secondary_power"]]), as.numeric(r[["overall_recovery"]]),
    as.numeric(r[["cand_exact"]]), as.numeric(r[["direction"]]),
    as.numeric(r[["count_exact"]]), as.numeric(r[["extra_rate"]]))),
  "",
  "- 定位：I-3 证明 locus 内 multi-trait 与 trait-specific signal 共存；E4 只作展示",
  "- E4 的 extra rate 预期偏高（novel-locus FDR 传播，Stage 7.5/8.5 已定性）"
)
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
cat("BATCH8 DONE\n")
