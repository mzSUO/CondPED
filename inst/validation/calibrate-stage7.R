## Stage 7.1 Phase B: targeted calibration grid (B1-B3).
## 49 parameter sets x 30 replicates, fixed calibration master seed.
## Statistical thresholds are never touched; only data-generating
## parameters (PVE multipliers, target_r2) vary. The full grid is
## preserved; no parameter set is dropped post hoc.
devtools::load_all(quiet = TRUE)

base <- CondPED:::.pilot_grid()
row_for <- function(scen) base[base$scenario == scen, , drop = FALSE]

mult <- c(1.0, 1.25, 1.5, 2.0)
grid_rows <- list()
add <- function(scen, lp_mult, sec_mult = NA_real_, tr2 = NA_real_) {
  r <- row_for(scen)
  r$locus_pve <- r$locus_pve * lp_mult
  if (!is.na(sec_mult)) r$secondary_signal_pve <- r$secondary_signal_pve * sec_mult
  if (!is.na(tr2)) r$target_r2 <- tr2
  r$parameter_set_id <- sprintf("%s|lp%.2f|sp%.2f|tr2%.1f", scen, lp_mult,
                                sec_mult, ifelse(is.na(tr2), -1, tr2))
  grid_rows[[length(grid_rows) + 1L]] <<- r
}

# B1: Simulation I PVE calibration
for (ml in mult) add("single_multi_trait", ml)
for (ml in mult) for (ms in mult) {
  add("two_heterogeneous", ml, ms)
  add("two_linked_trait_specific", ml, ms)
}
# B2: linked_pseudo_multitrait: target_r2 x secondary PVE multiplier
for (tr2 in c(0.3, 0.5, 0.7)) for (ms in c(1.0, 1.5, 2.0)) {
  add("linked_pseudo_multitrait", 1.0, ms, tr2)
}
# B3: single_nonredundant global scaling (truth rho map invariant)
for (ml in mult) add("single_nonredundant", ml)

cal_grid <- do.call(rbind, grid_rows)
rownames(cal_grid) <- NULL
psid <- cal_grid$parameter_set_id
cal_grid$parameter_set_id <- NULL   # label kept out of scenario_id

out_dir <- "inst/validation/output/calibration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(data.frame(parameter_set_id = psid, cal_grid),
        file.path(out_dir, "calibration_grid.rds"))

cat("calibration grid rows:", nrow(cal_grid), "\n")
t0 <- Sys.time()
res <- run_condped_simulation(
  cal_grid, reps = 30L, out_dir = out_dir,
  master_seed = 20260813L,
  workers = 4L, backend = "parallel", progress = FALSE
)
cat("WALL_SECONDS:", round(difftime(Sys.time(), t0, units = "secs"), 1), "\n")
saveRDS(res$manifest, file.path(out_dir, "manifest.rds"))
cat("DONE rows:", nrow(res$manifest), "\n")
print(table(res$manifest$status))
