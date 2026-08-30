## Stage 9 batch 6: boundary sensitivity (freeze §1.60).
## rho* in {0.08, 0.10, 0.12} around tau = 0.10 on R1, R=200/cell.
## rho* >= tau is known to be structurally unreachable for R1
## (invalid_covariance, Stage 7.6 P4); such cells are recorded as
## generation-infeasible — themselves the boundary finding. P(rho_hat
## crosses tau) is computed for feasible cells.
## Launch: bash run_sim.sh stage9-boundary inst/formal-simu/run-stage9-boundary.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-boundary"
tau <- 0.10
grid <- do.call(rbind, lapply(c(0.08, 0.10, 0.12), function(rho) {
  data.frame(
    experiment = "trait_representation",
    scenario = "highly_representable",
    analysis_mode = "signal_trait_oracle",
    comparison_pipeline = "condped_full",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = NA_real_,
    correlation = "block", target_loss = rho,
    tolerance = tau, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}))

final <- stage9_run_grid(grid, out_dir, reps_total = 200L, workers = 4L)

grid$key <- vapply(seq_len(nrow(grid)), function(i) stage9_key(grid[i, ]),
                   character(1))
final$dir <- basename(dirname(final$file))
final$rho <- grid$target_loss[match(final$dir, grid$key)]

rows <- do.call(rbind, lapply(final$file, function(f) {
  o <- readRDS(f)
  data.frame(rho = o$settings$target_loss, rep_id = o$rep_id,
             ok = isTRUE(o$status$ok),
             code = o$status$code, stringsAsFactors = FALSE)
}))
gen <- do.call(rbind, lapply(split(rows, rows$rho), function(d) {
  data.frame(rho = d$rho[1], n = nrow(d), n_ok = sum(d$ok),
             gen_feasible_rate = mean(d$ok),
             fail_codes = paste(unique(d$code[!d$ok]), collapse = ";"),
             stringsAsFactors = FALSE)
}))

## singleton rho_hat and P(cross tau) on feasible cells
cross <- do.call(rbind, lapply(final$file[final$status == "ok"], function(f) {
  o <- readRDS(f)
  tmap <- as.data.frame(o$truth$representation_map)
  s1 <- tmap[tmap$set_size == 1, ]
  truth_key <- s1$representing_key[which.min(s1$representation_loss)]
  st <- o$subset_table
  est_rho <- st$representation_loss[st$representing_key == truth_key][1]
  data.frame(rho = o$settings$target_loss,
             est_rho = if (length(est_rho)) est_rho else NA_real_,
             cross = isTRUE(est_rho >= tau), stringsAsFactors = FALSE)
}))
csum <- do.call(rbind, lapply(split(cross, cross$rho), function(d) {
  data.frame(rho = d$rho[1], n = nrow(d),
             abs_rho_minus_tau = abs(d$rho[1] - tau),
             P_cross = mean(d$cross, na.rm = TRUE), stringsAsFactors = FALSE)
}))
write.csv(gen, file.path(out_dir, "generation.csv"), row.names = FALSE)
write.csv(csum, file.path(out_dir, "crossing.csv"), row.names = FALSE)

lines <- c(
  "# Stage 9 batch 6: boundary sensitivity (freeze §1.60)",
  "",
  sprintf("- rho* in {0.08, 0.10, 0.12}, tau = 0.10, R=200/cell, %d/%d ok",
          sum(final$status == "ok"), nrow(final)),
  "",
  "## generation feasibility",
  "",
  "| rho* | n | n_ok | feasible rate | fail codes |",
  "|---|---|---|---|---|",
  apply(gen, 1, function(r) sprintf("| %.2f | %d | %d | %.3f | %s |",
    as.numeric(r[["rho"]]), as.integer(r[["n"]]), as.integer(r[["n_ok"]]),
    as.numeric(r[["gen_feasible_rate"]]), r[["fail_codes"]])),
  "",
  "## P(rho_hat crosses tau) on feasible cells",
  "",
  "| rho* | n | |rho*-tau| | P(cross) |",
  "|---|---|---|---|",
  apply(csum, 1, function(r) sprintf("| %.2f | %d | %.2f | %.3f |",
    as.numeric(r[["rho"]]), as.integer(r[["n"]]),
    as.numeric(r[["abs_rho_minus_tau"]]), as.numeric(r[["P_cross"]]))),
  "",
  "- 参考点（Stage 8.3 主场景，rho*=0.05, tau=0.10, n=500）：P(cross) = 0.278",
  "- 注意：rho* >= tau 的 cell 若 generation-infeasible，属于结构性事实：",
  "  R1 架构按定义要求 truth rho < tau，边界外侧的真值不可生成——",
  "  这正是 boundary sensitivity 要报告的机制边界"
)
writeLines(lines, file.path(out_dir, "report.md"))
print(gen, row.names = FALSE)
print(csum, row.names = FALSE)
cat("BATCH6 DONE\n")
