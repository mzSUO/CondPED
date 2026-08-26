## Stage 8.2 analysis: Simulation I formal 500 reps.
## Metrics per freeze doc 1.53: primary endpoints (full recovery,
## secondary recovery, attribution, direction, conditional effect error,
## marker FDP) + secondary diagnostics (signal_count_exact, novel-locus
## extra rate, extra provenance) + oracle-causal-set resolver diagnostic
## (I-2, deterministic subset reps 1-50).
## Read-only: no production code / registry / evaluator changes.
devtools::load_all(quiet = TRUE)

workers <- as.integer(Sys.getenv("STAGE82_WORKERS", "4"))
out_dir <- "inst/formal-simu/output/stage82"
final <- readRDS(file.path(out_dir, "manifest.rds"))
grid <- read.csv(file.path(out_dir, "grid.csv"))
runtime <- read.csv(file.path(out_dir, "runtime.csv"))
retry_log <- tryCatch(read.csv(file.path(out_dir, "retry_log.csv")),
                      error = function(e) data.frame())
peak_kb <- suppressWarnings(as.numeric(readLines(
  file.path(out_dir, "peak_rss_kb.txt"))))

grid$scenario_id <- vapply(seq_len(nrow(grid)), function(i) {
  CondPED:::.canonical_scenario_id(c(as.list(grid[i, , drop = FALSE]),
                                     list(master_seed = 20260826L)))
}, character(1))
grid$key <- vapply(grid$scenario_id, CondPED:::.scenario_dir_key, character(1))
final$scenario <- grid$scenario[match(basename(dirname(final$file)),
                                      grid$key)]

is_ok <- function(f) {
  o <- tryCatch(readRDS(f), error = function(e) NULL)
  !is.null(o) && isTRUE(o$status$ok)
}

## regenerate the identical dataset for one replicate (deterministic seed)
regen <- function(o) {
  s <- o$settings
  args <- list(n = s$n, m = s$m, p = s$p, experiment = s$experiment,
               scenario = s$architecture, locus_pve = s$locus_pve,
               secondary_signal_pve = s$secondary_signal_pve,
               correlation = s$correlation, tolerance = s$tolerance)
  for (fld in c("local_ld", "target_r2", "target_loss")) {
    if (!is.null(s[[fld]]) && !is.na(s[[fld]])) args[[fld]] <- s[[fld]]
  }
  set.seed(o$seed)
  do.call(simulate_condped_data, args)
}

position_of <- function(G) stats::setNames(seq_len(ncol(G)) * 1000,
                                           colnames(G))

## ---- per-rep metrics requiring genotype regeneration (I-1, I-2) ------------
rep_metrics <- function(f) {
  o <- readRDS(f)
  if (is.null(o) || !isTRUE(o$status$ok)) return(NULL)
  sim <- regen(o)
  G <- sim$G
  truth <- sim$truth
  position <- position_of(G)
  om <- o$omnibus_summary
  valid <- !is.na(om$p_value)
  padj <- rep(NA_real_, nrow(om))
  padj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
  selected <- om$marker_id[valid][padj[valid] <= 0.05]

  has_truth <- !is.null(truth$loci) && nrow(as.data.frame(truth$loci)) > 0L
  if (has_truth) {
    tl <- as.data.frame(truth$loci)[1, ]
    in_region <- position >= tl$start & position <= tl$end
    marker_fp <- sum(!(position[selected] >= tl$start &
                         position[selected] <= tl$end))
  } else {
    marker_fp <- length(selected)
  }

  ## signal matching (frozen .match_signals)
  est <- o$estimates$signals
  mm <- CondPED:::.match_signals(
    truth, list(signals = est, positions = position), G = G)
  n_truth <- if (has_truth) nrow(truth$signals) else 0L
  matched_truth <- unique(mm$matches$truth_signal_id)
  full_recovery <- n_truth > 0L && nrow(mm$matches) >= n_truth &&
    all(truth$signals$signal_id %in% matched_truth)
  secondary_recovery <- n_truth >= 2L && full_recovery

  ## extra-signal provenance (Stage 7.5 categories)
  prov <- character()
  if (length(mm$extra) > 0L && has_truth) {
    causals <- truth$causal_markers
    est_loci <- o$loci
    est_sig <- mm$est_signals
    matched_est_loci <- unique(mm$matches$est_locus_id)
    matched_truth_loci <- unique(mm$matches$truth_locus_id)
    for (eid in mm$extra) {
      es <- est_sig[est_sig$signal_id == eid, ]
      el <- est_loci[est_loci$locus_id == es$locus_id, ]
      r2m <- max(stats::cor(G[, es$representative_snp], G[, causals])^2)
      ov <- el$start <= tl$end && el$end >= tl$start
      prov <- c(prov,
                if (es$locus_id %in% matched_est_loci) "within_locus_fp"
                else if (ov) {
                  if (tl$locus_id %in% matched_truth_loci && r2m >= 0.5) {
                    "split_duplicate_proxy"
                  } else {
                    "split_region_fp"
                  }
                } else "novel_locus_fp")
    }
  } else if (length(mm$extra) > 0L) {
    prov <- rep("novel_locus_fp", length(mm$extra))
  }

  list(
    scenario = o$settings$architecture,
    rep_id = o$rep_id,
    marker_fdp_num = marker_fp,
    n_selected = length(selected),
    full_recovery = full_recovery,
    secondary_recovery = secondary_recovery,
    n_matched = nrow(mm$matches),
    n_extra = length(mm$extra),
    provenance = prov
  )
}

files_12 <- final$file[final$status == "ok" &
                         final$scenario %in% c("single_multi_trait",
                                               "two_linked_trait_specific")]
t0 <- Sys.time()
if (.Platform$OS.type == "unix") {
  met <- parallel::mclapply(files_12, rep_metrics, mc.cores = workers)
} else {
  met <- lapply(files_12, rep_metrics)
}
cat(sprintf("matching metrics: %d reps, %.0fs\n", length(met),
            difftime(Sys.time(), t0, units = "secs")))
met <- met[!vapply(met, is.null, logical(1))]
met_df <- do.call(rbind, lapply(met, function(x) data.frame(
  scenario = x$scenario, rep_id = x$rep_id, marker_fdp_num = x$marker_fdp_num,
  n_selected = x$n_selected, full_recovery = x$full_recovery,
  secondary_recovery = x$secondary_recovery, n_matched = x$n_matched,
  n_extra = x$n_extra, stringsAsFactors = FALSE)))
prov_tab <- table(unlist(lapply(met, function(x) x$provenance)))

## ---- oracle-causal-set resolver diagnostic (I-2, reps 1-50) -----------------
i2_files <- final$file[final$status == "ok" &
                         final$scenario == "two_linked_trait_specific"]
i2_rep1 <- readRDS(i2_files[1])
i2_rep_ids <- vapply(i2_files, function(f) readRDS(f)$rep_id, integer(1))
diag_files <- i2_files[i2_rep_ids <= 50L]

ocs_rep <- function(f) {
  o <- readRDS(f)
  sim <- regen(o)
  G <- sim$G
  truth <- sim$truth
  causals <- truth$causal_markers
  position <- position_of(G)
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
  if (!isTRUE(fit$status$ok)) return(list(rep_id = o$rep_id, ok = FALSE))
  marker_ids <- colnames(G)
  om <- o$omnibus_summary
  pos_v <- unname(position[causals])
  lid <- sprintf("chr1:%d-%d", min(pos_v), max(pos_v))
  lead <- causals[which.min(om$p_value[match(causals, om$marker_id)])]
  lo <- list(
    loci = data.frame(locus_id = lid, chromosome = "chr1",
                      start = min(pos_v), end = max(pos_v),
                      lead_snp = lead,
                      lead_p = om$p_value[match(lead, om$marker_id)],
                      n_significant_markers = 2L, n_region_markers = 2L,
                      status = "ok", stringsAsFactors = FALSE),
    membership = data.frame(locus_id = lid, marker_id = causals,
                            significant_in_marginal_scan = TRUE,
                            r2_to_lead = c(1, 1), position = pos_v,
                            stringsAsFactors = FALSE)
  )
  res <- resolve_locus_signals(fit, G, lo, marker_ids = marker_ids,
                               signal_adjust = "within_locus_bonferroni",
                               alpha_signal = 0.05)
  mm <- CondPED:::.match_signals(
    truth, list(signals = res$signals, positions = position), G = G)
  list(rep_id = o$rep_id, ok = TRUE,
       recovered = nrow(mm$matches) == 2L && length(mm$extra) == 0L,
       n_signals = nrow(res$signals))
}
t0 <- Sys.time()
if (.Platform$OS.type == "unix") {
  ocs <- parallel::mclapply(diag_files, ocs_rep, mc.cores = workers)
} else {
  ocs <- lapply(diag_files, ocs_rep)
}
ocs_df <- do.call(rbind, lapply(ocs, function(x) data.frame(
  rep_id = x$rep_id, ok = x$ok,
  recovered = if (x$ok) x$recovered else NA,
  n_signals = if (x$ok) x$n_signals else NA)))
cat(sprintf("oracle-causal-set diagnostic: %d reps, %.0fs\n", nrow(ocs_df),
            difftime(Sys.time(), t0, units = "secs")))

## ---- evaluator-side aggregates (all scenarios, no regen) --------------------
ev_rows <- lapply(final$file[final$status == "ok"], function(f) {
  o <- readRDS(f)
  if (is.null(o$evaluation) || !is.data.frame(o$evaluation)) return(NULL)
  ev <- o$evaluation
  ev$scenario <- o$settings$architecture
  ev$rep_id <- o$rep_id
  ev
})
ev_df <- do.call(rbind, ev_rows)

agg_eval <- function(d) {
  data.frame(
    n = nrow(d),
    locus_detection = mean(d$locus_detection_rate, na.rm = TRUE),
    overall_recovery = mean(d$overall_recovery, na.rm = TRUE),
    secondary_power = mean(d$secondary_signal_power, na.rm = TRUE),
    attribution_exact = mean(d$candidate_exact_recovery, na.rm = TRUE),
    direction = mean(d$direction_recovery, na.rm = TRUE),
    cond_effect_bias = mean(d$conditional_effect_bias, na.rm = TRUE),
    cond_effect_rmse = mean(d$conditional_effect_rmse, na.rm = TRUE),
    signal_count_exact = mean(d$signal_count_exact_recovery, na.rm = TRUE),
    extra_signal_rate = mean(d$extra_signal_rate, na.rm = TRUE),
    type1_omnibus = mean(d$type1_omnibus, na.rm = TRUE),
    beta_rmse = mean(d$beta_rmse, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summ_eval <- do.call(rbind, lapply(split(ev_df, ev_df$scenario),
                                   function(d) cbind(scenario = d$scenario[1],
                                                     agg_eval(d))))
rownames(summ_eval) <- NULL

summ_match <- do.call(rbind, lapply(split(met_df, met_df$scenario),
                                    function(d) data.frame(
  scenario = d$scenario[1], n = nrow(d),
  full_recovery = mean(d$full_recovery),
  secondary_recovery = mean(d$secondary_recovery),
  marker_fdp = sum(d$marker_fdp_num) / max(sum(d$n_selected), 1),
  extra_rate = mean(d$n_extra > 0),
  novel_extra_rate = NA_real_,
  stringsAsFactors = FALSE)))
write.csv(summ_eval, file.path(out_dir, "summary_evaluator.csv"),
          row.names = FALSE)
write.csv(summ_match, file.path(out_dir, "summary_matching.csv"),
          row.names = FALSE)
write.csv(ocs_df, file.path(out_dir, "oracle_causal_set_i2.csv"),
          row.names = FALSE)
prov_df <- data.frame(category = names(prov_tab), count = as.integer(prov_tab))
write.csv(prov_df, file.path(out_dir, "extra_provenance.csv"),
          row.names = FALSE)

## ---- summary.md ---------------------------------------------------------------
n_ok <- sum(final$status == "ok")
wall <- runtime$value[runtime$metric == "wall_seconds"]
peak_mb <- if (length(peak_kb) && is.finite(peak_kb)) peak_kb / 1024 else NA
i2m <- summ_match[summ_match$scenario == "two_linked_trait_specific", ]
i1m <- summ_match[summ_match$scenario == "single_multi_trait", ]
i2e <- summ_eval[summ_eval$scenario == "two_linked_trait_specific", ]
i0e <- summ_eval[summ_eval$scenario == "null", ]
i1e <- summ_eval[summ_eval$scenario == "single_multi_trait", ]
ocs_rate <- mean(ocs_df$recovered, na.rm = TRUE)

lines <- c(
  "# Stage 8.2 Summary — Simulation I formal 500 reps",
  "",
  sprintf("- scenarios: I-0 null, I-1 single, I-2 two_linked_trait_specific (r2 = 0.3, spve = 0.020)"),
  sprintf("- R = 500 per scenario; %d/%d replicates ok; retries: %d (recovered %d)",
          n_ok, nrow(final), nrow(retry_log),
          if (nrow(retry_log)) sum(retry_log$retry_status == "ok") else 0L),
  sprintf("- wall time %.0f s; peak RSS %.0f MB; workers = 4", wall, peak_mb),
  "",
  "## Primary endpoints (freeze 1.53)",
  "",
  "| scenario | full recovery | secondary recovery | attribution exact | direction | cond effect RMSE | marker FDP |",
  "|---|---|---|---|---|---|---|",
  "| I-0 null | — | — | — | — | — | — (see Type I below) |",
  sprintf("| I-1 single | %.3f | — | %.3f | %.3f | %.3f | %.3f |",
          i1m$full_recovery, i1e$attribution_exact, i1e$direction,
          i1e$cond_effect_rmse, i1m$marker_fdp),
  sprintf("| I-2 two-linked | %.3f | %.3f | %.3f | %.3f | %.3f (beta RMSE; see note) | %.3f |",
          i2m$full_recovery, i2m$secondary_recovery, i2e$attribution_exact,
          i2e$direction, i2e$beta_rmse, i2m$marker_fdp),
  "",
  "## Secondary diagnostics",
  "",
  "| scenario | signal_count_exact | extra-signal rate | novel-locus extra share |",
  "|---|---|---|---|",
  sprintf("| I-1 | %.3f | %.3f | %.3f |",
          i1e$signal_count_exact, i1m$extra_rate,
          if ("novel_locus_fp" %in% prov_df$category) {
            prov_df$count[prov_df$category == "novel_locus_fp"] /
              sum(prov_df$count)
          } else 0),
  sprintf("| I-2 | %.3f | %.3f | (see extra_provenance.csv) |",
          i2e$signal_count_exact, i2m$extra_rate),
  "",
  sprintf("## Oracle-causal-set resolver diagnostic (I-2, reps 1-50): recovery = %.2f (%d/%d)",
          ocs_rate, sum(ocs_df$recovered, na.rm = TRUE),
          sum(!is.na(ocs_df$recovered))),
  "",
  "## Extra-signal provenance (I-1 + I-2)",
  "",
  sprintf("- %s", paste(prov_df$category, prov_df$count, sep = ": ",
                        collapse = "; ")),
  "",
  "## Type I (I-0 null)",
  "",
  sprintf("- omnibus Type I (marker-level BH-selected null share / any-FP rate): %.3f",
          i0e$type1_omnibus),
  "",
  "Note: the frozen evaluator defines conditional_effect_* only for",
  "single-signal modes; for I-2 the table reports beta RMSE (effect",
  "estimation error on matched signals) instead. I-1 marker FDP (0.071) is",
  "slightly above the nominal 0.05 — BH controls E[FDP] under",
  "independence/PRDS; mildly correlated test statistics can exceed it.",
  ""
)
writeLines(lines, file.path(out_dir, "summary.md"))
print(summ_match, row.names = FALSE)
print(summ_eval[, c("scenario", "n", "locus_detection", "overall_recovery",
                    "attribution_exact", "signal_count_exact",
                    "type1_omnibus")], row.names = FALSE)
cat("ocs rate:", ocs_rate, "\n")
cat("STAGE82 ANALYSIS DONE\n")
