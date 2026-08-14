## Stage 7.2 tasks 4-6: sanity checks for single_highly_representable and
## mixed_multisignal, plus the linked_pseudo_multitrait redesign grid.
## 22 parameter sets x 30 reps. No statistical threshold is touched.
devtools::load_all(quiet = TRUE)

base <- CondPED:::.pilot_grid()
row_for <- function(scen) base[base$scenario == scen, , drop = FALSE]

grid_rows <- list()
add <- function(scen, lp, sp, tr2 = NA_real_, tag) {
  r <- row_for(scen)
  r$locus_pve <- lp
  r$secondary_signal_pve <- sp
  if (!is.na(tr2)) r$target_r2 <- tr2
  r$parameter_set_id <- tag
  grid_rows[[length(grid_rows) + 1L]] <<- r
}

# task 4: single_highly_representable, baseline vs locus_pve = 0.02
add("single_highly_representable", 0.01, 0.005, tag = "shr|baseline")
add("single_highly_representable", 0.02, 0.005, tag = "shr|lp0.02")

# task 5: mixed_multisignal, baseline vs lp 0.02 / sp 0.0075
add("mixed_multisignal", 0.01, 0.005, tag = "mixed|baseline")
add("mixed_multisignal", 0.02, 0.0075, tag = "mixed|lp0.02_sp0.0075")

# task 6: linked_pseudo_multitrait redesign grid
# truth always: signal1 -> Trait1 only, signal2 -> Trait2 only.
# common signal scale s: locus_pve = 0.01 * s
# beta2/beta1 ratio r: secondary_signal_pve = locus_pve * r^2
for (tr2 in c(0.4, 0.6, 0.7)) {
  for (s in c(1.0, 1.5, 2.0)) {
    for (ratio in c(0.75, 1.0)) {
      lp <- 0.01 * s
      sp <- lp * ratio^2
      add("linked_pseudo_multitrait", lp, sp, tr2,
          tag = sprintf("lpm|tr2%.1f|s%.1f|r%.2f", tr2, s, ratio))
    }
  }
}

g <- do.call(rbind, grid_rows)
rownames(g) <- NULL
psid <- g$parameter_set_id
g$parameter_set_id <- NULL

out_dir <- "inst/validation/output/stage72"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(data.frame(parameter_set_id = psid, g),
        file.path(out_dir, "grid.rds"))

cat("grid rows:", nrow(g), "\n")
t0 <- Sys.time()
res <- run_condped_simulation(
  g, reps = 30L, out_dir = out_dir,
  master_seed = 20260814L,
  workers = 4L, backend = "parallel", progress = FALSE
)
cat("WALL_SECONDS:", round(difftime(Sys.time(), t0, units = "secs"), 1), "\n")
saveRDS(res$manifest, file.path(out_dir, "manifest.rds"))
cat("DONE rows:", nrow(res$manifest), "\n")
print(table(res$manifest$status))
