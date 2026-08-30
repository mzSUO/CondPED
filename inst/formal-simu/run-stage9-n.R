## Stage 9 batch 2: sample-size sensitivity (freeze §1.55).
## two_linked_trait_specific at n in {500, 750, 1250}, R=200/cell, paired,
## deterministic seeds, workers=4.
## Launch: bash run_sim.sh stage9-n inst/formal-simu/run-stage9-n.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-n"
grid <- do.call(rbind, lapply(c(500L, 750L, 1250L), function(nn) {
  data.frame(
    experiment = "signal_resolution",
    scenario = "two_linked_trait_specific",
    analysis_mode = "full", comparison_pipeline = "paired",
    n = nn, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = 0.3,
    correlation = "block", target_loss = NA_real_,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}))

final <- stage9_run_grid(grid, out_dir, reps_total = 200L, workers = 4L)

grid$key <- vapply(seq_len(nrow(grid)), function(i) stage9_key(grid[i, ]),
                   character(1))
final$dir <- basename(dirname(final$file))
final$n <- grid$n[match(final$dir, grid$key)]

rows <- do.call(rbind, lapply(final$file[final$status == "ok"], function(f) {
  o <- readRDS(f)
  if (!is.data.frame(o$evaluation)) return(NULL)
  ev <- o$evaluation
  data.frame(n = o$settings$n,
             locus_detection = ev$locus_detection_rate,
             secondary_power = ev$secondary_signal_power,
             overall_recovery = ev$overall_recovery,
             cand_exact = ev$candidate_exact_recovery,
             direction = ev$direction_recovery,
             count_exact = ev$signal_count_exact_recovery,
             extra_rate = ev$extra_signal_rate,
             runtime = o$runtime,
             stringsAsFactors = FALSE)
}))
summ <- do.call(rbind, lapply(split(rows, rows$n), function(d) {
  data.frame(
    n = d$n[1], n_ok = nrow(d),
    locus_detection = stage9_prop(d$locus_detection)["estimate"],
    secondary_power = stage9_prop(d$secondary_power)["estimate"],
    overall_recovery = stage9_prop(d$overall_recovery)["estimate"],
    cand_exact = stage9_prop(d$cand_exact)["estimate"],
    direction = stage9_prop(d$direction)["estimate"],
    count_exact = stage9_prop(d$count_exact)["estimate"],
    extra_rate = stage9_prop(d$extra_rate)["estimate"],
    mean_runtime = mean(d$runtime, na.rm = TRUE),
    stringsAsFactors = FALSE)
}))
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
  "# Stage 9 batch 2: sample-size sensitivity (freeze §1.55)",
  "",
  sprintf("- n in {500, 750, 1250}, R=200/cell, paired, %d/%d ok",
          sum(final$status == "ok"), nrow(final)),
  "",
  "| n | locus det | secondary power | overall recovery | cand exact | direction | count_exact | extra rate | mean rep runtime (s) |",
  "|---|---|---|---|---|---|---|---|---|",
  apply(summ, 1, function(r) sprintf("| %d | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.1f |",
    as.integer(r[["n"]]), as.numeric(r[["locus_detection"]]),
    as.numeric(r[["secondary_power"]]), as.numeric(r[["overall_recovery"]]),
    as.numeric(r[["cand_exact"]]), as.numeric(r[["direction"]]),
    as.numeric(r[["count_exact"]]), as.numeric(r[["extra_rate"]]),
    as.numeric(r[["mean_runtime"]]))),
  "",
  "- reference: main n=1000 at R=500 (Stage 8.2): full/secondary 0.778",
  "- 核对点：recovery 随 n 单调不降；runtime 随 n 增长可接受"
)
lines <- c(lines, "",
  sprintf("- 严格匹配口径（Stage 8.2 约定）：full recovery %.3f, secondary %.3f",
          strict_rate, strict_sec))
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
cat("BATCH2 DONE\n")
