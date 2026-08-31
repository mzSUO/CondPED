## Stage 9 overall summary: aggregate all 12 supplementary batches.
## Read-only aggregation of per-batch reports/CSVs into
## inst/formal-simu/output/stage9/summary.md.
out_base <- "inst/formal-simu/output"
batches <- c("stage9-spve", "stage9-n", "stage9-r2", "stage9-rho",
             "stage9-tau", "stage9-boundary", "stage9-q3", "stage9-i3e4",
             "stage9-i4", "stage9-r2arch", "stage9-fdp", "stage9-iib-i4")

lines <- c(
  "# Stage 9 Summary — Supplementary Sensitivity Series",
  "",
  "All batches: deterministic seeds (master 20260826), per-rep RDS",
  "checkpoints, workers = 4, one automatic retry, 0 numerical failure",
  "hard standard. Section numbers refer to the current freeze document",
  "(05 模拟20260825.md, constraint fix v2).",
  ""
)
for (b in batches) {
  rep <- file.path(out_base, b, "report.md")
  if (file.exists(rep)) {
    lines <- c(lines, paste0("## ", b), "", readLines(rep), "")
  }
}
dir.create(file.path(out_base, "stage9"), showWarnings = FALSE)
writeLines(lines, file.path(out_base, "stage9", "summary.md"))
cat("STAGE9 SUMMARY DONE\n")
