## Stage 7.3: run mixed_multisignal simulation for extra-signal diagnosis.
## Uses the formal runner so output is compatible with
## diagnose-stage73-mixed-extra.R.
devtools::load_all(quiet = TRUE)

base <- CondPED:::.pilot_grid()
row <- base[base$scenario == "mixed_multisignal", , drop = FALSE]
row$locus_pve <- 0.02
row$secondary_signal_pve <- 0.0075

out_dir <- "inst/validation/output/stage73/mixed_extra"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
res <- run_condped_simulation(
  row, reps = 100L, out_dir = out_dir,
  master_seed = 20260815L,
  workers = 4L, backend = "parallel", progress = FALSE
)
wall <- difftime(Sys.time(), t0, units = "secs")
cat("WALL_SECONDS:", round(wall, 1), "\n")
saveRDS(res$manifest, file.path(out_dir, "manifest.rds"))
cat("DONE rows:", nrow(res$manifest), "\n")
print(table(res$manifest$status))
