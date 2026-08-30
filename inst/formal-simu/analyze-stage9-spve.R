## Stage 9 batch 1 report (standalone; the run completed 600/600 ok).
## Regenerates summary.csv / per_rep.csv / report.md from the manifest.
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-spve"
grid <- read.csv(file.path(out_dir, "grid.csv"))
final <- readRDS(file.path(out_dir, "manifest.rds"))

grid$key <- vapply(seq_len(nrow(grid)), function(i) stage9_key(grid[i, ]),
                   character(1))
final$scenario_id_dir <- basename(dirname(final$file))
final$spve <- grid$secondary_signal_pve[match(final$scenario_id_dir, grid$key)]

rows <- do.call(rbind, lapply(final$file[final$status == "ok"], function(f) {
  o <- readRDS(f)
  if (!is.data.frame(o$evaluation)) return(NULL)
  ev <- o$evaluation
  data.frame(spve = o$settings$secondary_signal_pve,
             locus_detection = ev$locus_detection_rate,
             overall_recovery = ev$overall_recovery,
             secondary_power = ev$secondary_signal_power,
             cand_exact = ev$candidate_exact_recovery,
             direction = ev$direction_recovery,
             count_exact = ev$signal_count_exact_recovery,
             extra_rate = ev$extra_signal_rate,
             stringsAsFactors = FALSE)
}))
summ <- do.call(rbind, lapply(split(rows, rows$spve), function(d) {
  data.frame(
    spve = d$spve[1], n = nrow(d),
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
  "# Stage 9 batch 1: secondary-signal PVE sensitivity (freeze §1.56)",
  "",
  sprintf("- spve grid {0.0075, 0.015, 0.0175}, R=200/cell, paired, %d/%d ok",
          sum(final$status == "ok"), nrow(final)),
  "",
  "| spve | locus det | secondary power | overall recovery | cand exact | direction | count_exact | extra rate |",
  "|---|---|---|---|---|---|---|---|",
  apply(summ, 1, function(r) sprintf("| %.4f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |",
    as.numeric(r[["spve"]]), as.numeric(r[["locus_detection"]]),
    as.numeric(r[["secondary_power"]]), as.numeric(r[["overall_recovery"]]),
    as.numeric(r[["cand_exact"]]), as.numeric(r[["direction"]]),
    as.numeric(r[["count_exact"]]), as.numeric(r[["extra_rate"]]))),
  "",
  sprintf("- reference: main spve=0.020 at R=500 (Stage 8.2): secondary 0.778, full 0.778, extra 0.486"),
  "- 结论核对点：secondary power 应随 spve 单调上升"
)
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
cat("BATCH1 REPORT DONE\n")
