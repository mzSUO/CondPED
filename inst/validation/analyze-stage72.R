## Stage 7.2 analysis: tasks 4 (shr sanity), 5 (mixed sanity),
## 6 (linked_pseudo_multitrait redesign selection).
## Selection criteria were fixed BEFORE the runs:
##   locus detection >= 0.70
##   P(lead breadth inflated by the OTHER CAUSAL trait) >= 0.40
##   P(resolution improves breadth) >= 0.50
##   extra signal rate <= 0.15
devtools::load_all(quiet = TRUE)

out_dir <- "inst/validation/output/stage72"
grid <- readRDS(file.path(out_dir, "grid.rds"))
manifest <- readRDS(file.path(out_dir, "manifest.rds"))
grid$scenario_id <- vapply(seq_len(nrow(grid)), function(i) {
  CondPED:::.canonical_scenario_id(as.list(grid[i, setdiff(names(grid), "parameter_set_id")]))
}, character(1))

ev_of <- function(psid) {
  sid <- grid$scenario_id[grid$parameter_set_id == psid]
  files <- manifest$file[manifest$scenario_id == sid &
                           manifest$status == "ok"]
  list(files = files,
       ev = evaluate_condped_simulation(files, group_by = character()))
}

## ---------------- tasks 4 & 5: sanity comparisons -----------------------------
cat("=== task 4: single_highly_representable ===\n")
for (psid in c("shr|baseline", "shr|lp0.02")) {
  s <- ev_of(psid)$ev$summary
  cat(psid, ": n=", s$n,
      " det=", round(s$locus_detection_rate, 3),
      " cand_exact=", round(s$candidate_exact_recovery, 3),
      " rep_exact=", round(s$rep_exact_family_recovery, 3),
      " rho_rmse=", round(s$representation_loss_rmse, 4),
      " split=", round(s$truth_locus_split_rate, 3),
      " extra=", round(s$extra_signal_rate, 3), "\n", sep = "")
}
cat("\n=== task 5: mixed_multisignal ===\n")
for (psid in c("mixed|baseline", "mixed|lp0.02_sp0.0075")) {
  s <- ev_of(psid)$ev$summary
  cat(psid, ": n=", s$n,
      " det=", round(s$locus_detection_rate, 3),
      " cond_exact=", round(s$conditional_signal_count_exact_recovery, 3),
      " sec_power=", round(s$conditional_secondary_signal_power, 3),
      " split=", round(s$truth_locus_split_rate, 3),
      " cand_exact=", round(s$candidate_exact_recovery, 3),
      " rep_exact=", round(s$rep_exact_family_recovery, 3),
      " extra=", round(s$extra_signal_rate, 3), "\n", sep = "")
}

## ---------------- task 6: linked_pseudo redesign ------------------------------
cat("\n=== task 6: linked_pseudo_multitrait redesign grid ===\n")
lpm <- grid[grepl("^lpm\\|", grid$parameter_set_id), ]
res_rows <- list()
for (i in seq_len(nrow(lpm))) {
  psid <- lpm$parameter_set_id[[i]]
  x <- ev_of(psid)
  s <- x$ev$summary
  # per-replicate: valid lead inflation (subset contains BOTH causal
  # traits Trait1 and Trait2); paired breadth improvement comes from
  # the evaluator's match-aware per-replicate breadth columns
  r <- x$ev$replicate_metrics
  li <- r$lead_trait_breadth_inflation
  ri <- r$resolved_trait_breadth_inflation
  both <- !is.na(li) & !is.na(ri)
  p_improve <- if (any(both)) {
    mean(abs(li[both]) > abs(ri[both]))
  } else {
    NA_real_
  }
  n_inflated <- 0L
  n_reps <- length(x$files)
  for (f in x$files) {
    o <- readRDS(f)
    truth_sets <- o$truth$candidate_traits
    causal_traits <- unique(unlist(truth_sets))   # {Trait1, Trait2}
    lead <- o$asset_results$lead_asset$loci
    lead_sub <- if (length(lead) > 0L) lead[[1]]$asset_best_subset else character()
    if (is.null(lead_sub)) lead_sub <- character()
    if (length(causal_traits) > 0L && all(causal_traits %in% lead_sub)) {
      n_inflated <- n_inflated + 1L
    }
  }
  res_rows[[i]] <- data.frame(
    parameter_set_id = psid,
    locus_pve = lpm$locus_pve[[i]],
    secondary_signal_pve = round(lpm$secondary_signal_pve[[i]], 5),
    target_r2 = lpm$target_r2[[i]],
    n_reps = n_reps,
    locus_detection = s$locus_detection_rate,
    P_valid_lead_inflation = n_inflated / n_reps,
    P_resolution_improves = p_improve,
    extra_signal_rate = s$extra_signal_rate,
    cond_count_exact = s$conditional_signal_count_exact_recovery,
    cond_secondary_power = s$conditional_secondary_signal_power,
    stringsAsFactors = FALSE
  )
}
tab <- do.call(rbind, res_rows)
tab$pass_detection <- tab$locus_detection >= 0.70
tab$pass_inflation <- tab$P_valid_lead_inflation >= 0.40
tab$pass_improvement <- tab$P_resolution_improves >= 0.50
tab$pass_extra <- tab$extra_signal_rate <= 0.15
tab$accepted <- with(tab, pass_detection & pass_inflation &
                       pass_improvement & pass_extra)
print(tab, row.names = FALSE)
write.csv(tab, file.path(out_dir, "lpm_redesign_table.csv"),
          row.names = FALSE)
cat("\naccepted parameter sets:",
    sum(tab$accepted), "of", nrow(tab), "\n")
cat("\nSTAGE 7.2 ANALYSIS DONE\n")
