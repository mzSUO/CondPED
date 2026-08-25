## Stage 8.1 pilot analysis: verify the 10 freeze pilot checks
## (freeze doc 05 模拟20260825.md section 1.62) against the per-replicate
## RDS checkpoints and write pilot_report.md. Read-only analysis.
devtools::load_all(quiet = TRUE)

out_dir <- "inst/formal-simu/output/stage81"
grid <- read.csv(file.path(out_dir, "grid.csv"))
final <- readRDS(file.path(out_dir, "manifest.rds"))
runtime <- read.csv(file.path(out_dir, "runtime.csv"))
missing_check <- read.csv(file.path(out_dir, "missing_output_check.csv"))
retry_log <- tryCatch({
  rl <- read.csv(file.path(out_dir, "retry_log.csv"))
  rl
}, error = function(e) data.frame())
peak_kb <- suppressWarnings(
  as.numeric(readLines(file.path(out_dir, "peak_rss_kb.txt"))))

## map scenario dir keys back to grid scenarios (row-wise, keeping types:
## apply() would coerce numerics to character and break the canonical id)
sid_of <- function(row) {
  settings <- as.list(row)
  settings$master_seed <- 20260825L
  CondPED:::.canonical_scenario_id(settings)
}
grid$scenario_id <- vapply(seq_len(nrow(grid)), function(i) {
  sid_of(grid[i, , drop = FALSE])
}, character(1))
grid$key <- vapply(grid$scenario_id, CondPED:::.scenario_dir_key,
                   character(1))
final$scenario_dir <- basename(dirname(final$file))
final <- merge(final, grid[, c("key", "scenario", "experiment",
                               "analysis_mode", "comparison_pipeline")],
               by.x = "scenario_dir", by.y = "key")

read_rep <- function(f) tryCatch(readRDS(f), error = function(e) NULL)

checks <- list()
add_check <- function(item, pass, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    item = item, result = if (isTRUE(pass)) "PASS" else "FAIL",
    evidence = evidence, stringsAsFactors = FALSE
  )
}

## ---- 1. realized LD (I-2, E3 target r2 = 0.3) ------------------------------
ld_rows <- final[final$scenario %in% c("two_linked_trait_specific",
                                       "linked_pseudo_multitrait") &
                   final$status == "ok", ]
r2s <- unlist(lapply(ld_rows$file, function(f) {
  o <- read_rep(f)
  if (is.null(o) || is.null(o$truth$local_ld)) return(NULL)
  o$truth$local_ld$realized_mean_r2
}))
med_r2 <- stats::median(r2s)
add_check("1. realized LD (target 0.3)",
          abs(med_r2 - 0.3) < 0.05,
          sprintf("median realized mean r2 = %.3f over %d linked reps (range %.3f-%.3f)",
                  med_r2, length(r2s), min(r2s), max(r2s)))

## ---- 2. signal PVE ----------------------------------------------------------
## realized_marginal_signal_pve is the per-signal mean over ALL traits
## (including zero-effect traits); the frozen calibration target is the
## per-ACTIVE-trait PVE, so normalize by m / (#active traits). scale_match
## scenarios (R1/R3/E1/E2) are calibrated on the Q-form by design and are
## reported separately.
pve_rows <- lapply(final$file[final$status == "ok"], function(f) {
  o <- read_rep(f)
  if (is.null(o)) return(NULL)
  rp <- o$truth$realized_signal_pve
  if (is.null(rp) || nrow(rp) == 0L) return(NULL)
  n_act <- vapply(o$truth$candidate_traits[rp$signal_id], length, integer(1))
  data.frame(arch = o$settings$architecture, target = rp$target_pve,
             realized = rp$realized_marginal_signal_pve,
             n_active = n_act, stringsAsFactors = FALSE)
})
pve_df <- do.call(rbind, pve_rows)
m_traits <- 4L
scale_match_arch <- c("highly_representable", "strongly_nonredundant",
                      "single_highly_representable", "single_nonredundant")
fx <- pve_df[!pve_df$arch %in% scale_match_arch, ]
sm <- pve_df[pve_df$arch %in% scale_match_arch, ]
fx_adj <- fx$realized * m_traits / fx$n_active
add_check("2. signal PVE (target 0.02, per-active-trait)",
          abs(stats::median(fx_adj) - 0.02) < 0.002,
          sprintf(paste0("fixed-dirs scenarios: median realized per-active-trait ",
                         "PVE = %.4f (n=%d signals); scale-match scenarios by design: ",
                         "median realized total PVE = %.4f (Q-form matched)"),
                  stats::median(fx_adj), nrow(fx), stats::median(sm$realized)))

## ---- 3. Sigma_G_bg PSD -------------------------------------------------------
min_eigs <- vapply(final$file[final$status == "ok"], function(f) {
  o <- read_rep(f)
  if (is.null(o) || is.null(o$truth$Sigma_G_bg)) return(NA_real_)
  CondPED:::.min_eigen_sym(o$truth$Sigma_G_bg)
}, numeric(1))
add_check("3. Sigma_G^bg PSD", all(min_eigs >= -1e-8, na.rm = TRUE),
          sprintf("min eigenvalue over %d reps = %.3e", sum(!is.na(min_eigs)),
                  min(min_eigs, na.rm = TRUE)))

## ---- 4. truth rho map --------------------------------------------------------
rho_bad <- 0L
rho_seen <- 0L
for (f in final$file[final$status == "ok"]) {
  o <- read_rep(f)
  if (is.null(o) || is.null(o$truth$representation_map)) next
  rm_ <- o$truth$representation_map
  rho_seen <- rho_seen + 1L
  rl <- rm_$representation_loss
  if (any(!is.finite(rl)) || any(rl < -1e-8)) rho_bad <- rho_bad + 1L
}
add_check("4. truth rho map",
          rho_bad == 0L,
          sprintf("%d reps with representation_map, %d with non-finite/negative rho",
                  rho_seen, rho_bad))

## ---- 5. Rep/Irr truth --------------------------------------------------------
r13 <- final[final$scenario %in% c("highly_representable",
                                   "strongly_nonredundant") &
               final$status == "ok", ]
rep_ok <- 0L; rep_bad <- 0L
for (i in seq_len(nrow(r13))) {
  o <- read_rep(r13$file[i])
  if (is.null(o)) { rep_bad <- rep_bad + 1L; next }
  mrs <- o$truth$minimum_representative_sets
  if (is.null(mrs) || nrow(mrs) == 0L) { rep_bad <- rep_bad + 1L; next }
  minsz <- min(mrs$set_size)
  want <- if (r13$scenario[i] == "highly_representable") 1L else 4L
  if (minsz == want) rep_ok <- rep_ok + 1L else rep_bad <- rep_bad + 1L
}
add_check("5. Rep/Irr truth (R1 min Rep=1, R3 min Rep=4)",
          rep_bad == 0L && rep_ok == nrow(r13),
          sprintf("R1/R3 reps checked: %d, truth min-Rep size correct: %d",
                  nrow(r13), rep_ok))

## ---- 6. ASSET correlation matrix --------------------------------------------
paired <- final[final$comparison_pipeline == "paired" & final$status == "ok", ]
asset_ok <- 0L; asset_bad <- 0L; asset_not_run <- 0L
for (f in paired$file) {
  o <- read_rep(f)
  if (is.null(o)) { asset_bad <- asset_bad + 1L; next }
  ar <- o$asset_results
  if (is.null(ar) || is.null(ar$lead_asset) || is.null(ar$resolved_asset)) {
    ## discovery found no locus -> ASSET never ran (a discovery outcome,
    ## not an ASSET/correlation-matrix failure)
    asset_not_run <- asset_not_run + 1L
    next
  }
  ## adapter validates Sigma_Z (symmetric, unit diagonal, PSD) internally;
  ## status ok <=> a valid trait-correlation matrix was accepted by ASSET
  if (identical(ar$lead_asset$status, "ok") &&
      identical(ar$resolved_asset$status, "ok")) {
    asset_ok <- asset_ok + 1L
  } else {
    asset_bad <- asset_bad + 1L
  }
}
add_check("6. ASSET correlation matrix",
          asset_bad == 0L,
          sprintf("paired reps: %d; ASSET ran: %d (both adapters ok); no-locus (ASSET not run): %d; adapter failures: %d",
                  nrow(paired), asset_ok, asset_not_run, asset_bad))

## ---- 7. runtime ---------------------------------------------------------------
wall_main <- runtime$value[runtime$metric == "wall_main_seconds"]
wall_tot <- runtime$value[runtime$metric == "wall_total_seconds"]
add_check("7. runtime", is.finite(wall_main),
          sprintf("main grid wall = %.0fs, total incl. retry/resume/probe = %.0fs (240 reps, 4 workers)",
                  wall_main, wall_tot))

## ---- 8. signal matching --------------------------------------------------------
mm_bad <- 0L; mm_seen <- 0L
for (f in final$file[final$status == "ok"]) {
  o <- read_rep(f)
  if (is.null(o)) next
  ev <- o$evaluation
  if (is.null(ev) || !is.data.frame(ev)) next
  mm_seen <- mm_seen + 1L
  r2m <- ev$representative_to_causal_r2
  if (length(r2m) && any(is.finite(r2m)) &&
      (any(r2m[is.finite(r2m)] < -1e-8) ||
       any(r2m[is.finite(r2m)] > 1 + 1e-8))) {
    mm_bad <- mm_bad + 1L
  }
}
add_check("8. signal matching",
          mm_bad == 0L,
          sprintf("%d evaluated reps; representative-to-causal r2 out-of-range violations: %d",
                  mm_seen, mm_bad))

## ---- 9. missing-output handling -------------------------------------------------
add_check("9. missing-output handling",
          isTRUE(missing_check$n_skipped_ok_on_resume >= 200) &&
            isTRUE(missing_check$probe_regenerated) &&
            isTRUE(missing_check$probe_content_identical),
          sprintf("resume skipped %d ok reps; deleted probe regenerated: %s; content identical: %s",
                  missing_check$n_skipped_ok_on_resume,
                  missing_check$probe_regenerated,
                  missing_check$probe_content_identical))

## ---- 10. numerical failure rate --------------------------------------------------
n_fail <- sum(final$status != "ok")
add_check("10. numerical failure rate == 0",
          n_fail == 0L,
          sprintf("%d/%d replicates ok; retries needed: %d (recovered: %d)",
                  sum(final$status == "ok"), nrow(final),
                  nrow(retry_log),
                  if (nrow(retry_log)) sum(retry_log$retry_status == "ok") else 0L))

## ---- report -----------------------------------------------------------------------
chk <- do.call(rbind, checks)
peak_mb <- if (length(peak_kb) && is.finite(peak_kb)) peak_kb / 1024 else NA_real_
per_scen <- aggregate(status ~ scenario + experiment, data = final,
                      function(s) sprintf("%d/%d ok", sum(s == "ok"),
                                          length(s)))

lines <- c(
  "# Stage 8.1 Pilot Report",
  "",
  sprintf("- grid: %d scenarios x R = 20 reps (freeze doc 05 模拟20260825.md)",
          nrow(grid)),
  sprintf("- launcher: bash run_sim.sh stage81 inst/formal-simu/run-stage81.R"),
  sprintf("- main grid wall time: %.0f s; total (incl. retry/resume/probe): %.0f s",
          wall_main, wall_tot),
  sprintf("- peak RSS (master + workers): %.0f MB", peak_mb),
  sprintf("- replicates: %d ok / %d total", sum(final$status == "ok"),
          nrow(final)),
  "",
  "## Pilot checks (freeze 1.62)",
  "",
  "| # | check | result | evidence |",
  "|---|-------|--------|----------|",
  sprintf("| %s | %s | %s | %s |",
          seq_len(nrow(chk)), chk$item, chk$result, chk$evidence),
  "",
  "## Per-scenario completion",
  "",
  sprintf("- %s", paste(per_scen$scenario, per_scen$status,
                        sep = ": ", collapse = "; ")),
  ""
)
writeLines(lines, file.path(out_dir, "pilot_report.md"))
print(chk, row.names = FALSE)
cat("\npeak RSS MB:", round(peak_mb), "\n")
cat("PILOT ANALYSIS DONE\n")
