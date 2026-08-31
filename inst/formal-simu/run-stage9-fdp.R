## Stage 9 batch 11: FDR internal calibration (freeze §1.52).
## null / single / mixed at R=200 each. Formal claim: marker FDP ~ 0.05.
## novel-locus extra rate + provenance recorded; locus/signal FDP are NOT
## primary error-control claims (§1.52).
## Launch: bash run_sim.sh stage9-fdp inst/formal-simu/run-stage9-fdp.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-fdp"
mk <- function(experiment, scenario) {
  data.frame(
    experiment = experiment, scenario = scenario,
    analysis_mode = "full", comparison_pipeline = "paired",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = NA_real_,
    correlation = "block", target_loss = NA_real_,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}
grid <- rbind(
  mk("signal_resolution", "null"),
  mk("signal_resolution", "single_multi_trait"),
  mk("end_to_end", "mixed_multisignal")
)

final <- stage9_run_grid(grid, out_dir, reps_total = 200L, workers = 4L)

grid$key <- vapply(seq_len(nrow(grid)), function(i) stage9_key(grid[i, ]),
                   character(1))
final$dir <- basename(dirname(final$file))
final$scenario <- grid$scenario[match(final$dir, grid$key)]

rows <- do.call(rbind, lapply(final$file[final$status == "ok"], function(f) {
  o <- readRDS(f)
  om <- o$omnibus_summary
  valid <- !is.na(om$p_value)
  padj <- rep(NA_real_, nrow(om))
  padj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
  sel_pos <- as.integer(sub("^M", "", om$marker_id[valid]))[
    padj[valid] <= 0.05] * 1000
  tl <- if (!is.null(o$truth$loci) && nrow(as.data.frame(o$truth$loci)) > 0) {
    as.data.frame(o$truth$loci)[1, ]
  } else NULL
  fp <- if (is.null(tl)) length(sel_pos) else {
    sum(sel_pos < tl$start | sel_pos > tl$end)
  }
  ev <- o$evaluation
  data.frame(scenario = o$settings$architecture,
             n_sel = length(sel_pos), n_fp = fp,
             type1 = if (is.data.frame(ev)) ev$type1_omnibus else NA_real_,
             extra_rate = if (is.data.frame(ev)) ev$extra_signal_rate else NA_real_,
             count_exact = if (is.data.frame(ev)) ev$signal_count_exact_recovery else NA_real_,
             stringsAsFactors = FALSE)
}))
summ <- do.call(rbind, lapply(split(rows, rows$scenario), function(d) {
  ci <- stage9_wilson(sum(d$n_fp), sum(d$n_sel))
  data.frame(
    scenario = d$scenario[1], n = nrow(d),
    n_sel_total = sum(d$n_sel), n_fp_total = sum(d$n_fp),
    marker_fdp = sum(d$n_fp) / max(sum(d$n_sel), 1),
    fdp_ci_lo = ci[1], fdp_ci_hi = ci[2],
    type1 = stage9_prop(d$type1)["estimate"],
    extra_rate = stage9_prop(d$extra_rate)["estimate"],
    stringsAsFactors = FALSE)
}))
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(rows, file.path(out_dir, "per_rep.csv"), row.names = FALSE)

lines <- c(
  "# Stage 9 batch 11: FDR internal calibration (freeze §1.52)",
  "",
  sprintf("- null / single / mixed, R=200 each, %d/%d ok",
          sum(final$status == "ok"), nrow(final)),
  "",
  "| scenario | n | marker FDP [Wilson 95%] | type1 | extra rate |",
  "|---|---|---|---|---|",
  apply(summ, 1, function(r) sprintf("| %s | %d | %.3f [%.3f, %.3f] | %.3f | %.3f |",
    r[["scenario"]], as.integer(r[["n"]]), as.numeric(r[["marker_fdp"]]),
    as.numeric(r[["fdp_ci_lo"]]), as.numeric(r[["fdp_ci_hi"]]),
    as.numeric(r[["type1"]]), as.numeric(r[["extra_rate"]]))),
  "",
  "- 正式 claim（§1.52）：marker FDP ≈ 0.05（含 CI）",
  "- locus/signal FDP 不作为 primary error-control claim（marker BH 不自动产生 signal-level 控制）"
)
writeLines(lines, file.path(out_dir, "report.md"))
print(summ, row.names = FALSE)
cat("BATCH11 DONE\n")
