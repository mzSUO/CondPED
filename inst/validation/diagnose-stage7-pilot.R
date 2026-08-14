## Stage 7.1 Phase A diagnostics: layer-by-layer failure decomposition,
## linked_pseudo_multitrait LD-mixing audit, Sim II-A per-trait audit.
## Reads the frozen 20-replicate pilot output; recomputes nothing except
## genotype regeneration from the saved seed for LD diagnostics.
devtools::load_all(quiet = TRUE)

pilot_dir <- "inst/validation/output/pilot"
manifest <- readRDS(file.path(pilot_dir, "manifest.rds"))
manifest$scenario <- sub(".*(?:^|\\|)scenario=([^|]*).*", "\\1",
                         manifest$scenario_id)

eval_files <- function(files) {
  ev <- evaluate_condped_simulation(files, group_by = "architecture")
  ev$replicate_metrics
}

## ---------------------------------------------------------------- A2 ----
## Layer chain per replicate: L0 detection, L1 count exact, L2 matching,
## L3 candidate exact, L4 direction, L5 Rep cardinality, L6 Rep family.
layer_chain <- function(label, scenarios) {
  cat("\n=== A2 layer decomposition:", label, "===\n")
  out <- list()
  for (sc in scenarios) {
    files <- manifest$file[manifest$scenario == sc]
    r <- eval_files(files)
    L0 <- !is.na(r$locus_detection_rate) & r$locus_detection_rate > 0
    L1 <- r$conditional_signal_count_exact_recovery == 1
    L2 <- r$signal_coverage == 1 & r$extra_signal_rate == 0
    L3 <- r$candidate_exact_recovery == 1
    L4 <- r$direction_recovery == 1
    L5 <- r$rep_min_cardinality_recovery == 1
    L6 <- r$rep_exact_family_recovery == 1
    pc <- function(num, den) {
      den <- den & !is.na(num)
      if (sum(den) == 0) return(NA_real_)
      mean(num[den], na.rm = TRUE)
    }
    chain <- c(
      n = nrow(r),
      P_L0 = mean(L0),
      P_L1_given_L0 = pc(L1, L0),
      P_L2_given_L1 = pc(L2, L0 & L1),
      P_L3_given_L2 = pc(L3, L0 & L1 & L2),
      P_L4_given_L3 = pc(L4, L3),
      P_L5_given_L3 = pc(L5, L3),
      P_L6_given_L3 = pc(L6, L3),
      P_L3_uncond = mean(L3[!is.na(L3)]),
      P_L6_uncond = mean(L6[!is.na(L6)])
    )
    out[[sc]] <- chain
    cat("--", sc, "--\n")
    print(round(chain, 3))
  }
  invisible(out)
}

a2_sim1 <- layer_chain("Simulation I",
                       c("single_multi_trait", "two_heterogeneous",
                         "two_linked_trait_specific"))
a2_sim3 <- layer_chain("Simulation III",
                       c("single_highly_representable",
                         "single_nonredundant",
                         "linked_pseudo_multitrait", "mixed_multisignal"))

## conditional vs unconditional resolution metrics (Simulation I)
cat("\n=== A1 conditional vs unconditional (Simulation I) ===\n")
files <- manifest$file[manifest$scenario %in%
                         c("single_multi_trait", "two_heterogeneous",
                           "two_linked_trait_specific")]
ev <- evaluate_condped_simulation(files, group_by = "architecture")
keep <- c("signal_resolution_attempt_rate",
          "signal_count_exact_recovery",
          "conditional_signal_count_exact_recovery",
          "secondary_signal_power", "conditional_secondary_signal_power",
          "missed_signal_rate", "conditional_missed_signal_rate",
          "extra_signal_rate", "locus_detection_rate", "power_omnibus")
print(ev$summary[, c("architecture", "n",
                     intersect(keep, names(ev$summary)))],
      row.names = FALSE)

## ---------------------------------------------------------------- A3 ----
cat("\n=== A3 linked_pseudo_multitrait per-replicate diagnostic ===\n")
files <- manifest$file[manifest$scenario == "linked_pseudo_multitrait"]

regen_G <- function(o) {
  s <- o$settings
  args <- list(n = s$n, m = s$m, p = s$p,
               experiment = s$experiment, scenario = s$architecture,
               locus_pve = s$locus_pve,
               secondary_signal_pve = s$secondary_signal_pve,
               correlation = s$correlation, tolerance = s$tolerance)
  for (f in c("local_ld", "target_r2", "target_loss")) {
    if (!is.null(s[[f]]) && !is.na(s[[f]])) args[[f]] <- s[[f]]
  }
  set.seed(o$seed)   # runner sets the replicate seed before simulating
  do.call(simulate_condped_data, args)
}

a3_rows <- list()
for (f in files) {
  o <- readRDS(f)
  sim <- regen_G(o)
  G <- sim$G
  truth <- o$truth
  cm <- truth$causal_markers
  r2_pairs <- if (length(cm) >= 2L) {
    stats::cor(G[, cm])^2
  } else {
    matrix(NA_real_)
  }
  lead <- o$asset_results$lead_asset$loci
  res <- o$asset_results$resolved_asset$signals
  lead_e <- if (length(lead) > 0L) lead[[1]] else NULL
  lead_snp <- if (!is.null(lead_e)) lead_e$lead_snp else NA_character_
  if (is.null(lead_snp)) lead_snp <- NA_character_
  r2_to_causal <- if (!is.na(lead_snp) && length(cm) > 0L) {
    stats::setNames(stats::cor(G[, lead_snp], G[, cm])^2, cm)
  } else {
    stats::setNames(rep(NA_real_, length(cm)), cm)
  }
  t_breadth <- vapply(truth$candidate_traits, length, integer(1))
  lead_w <- if (!is.null(lead_e) && !is.null(lead_e$asset_best_subset)) {
    length(lead_e$asset_best_subset)
  } else {
    0L
  }
  res_w <- vapply(res, function(x) {
    if (is.null(x$asset_best_subset)) 0L else length(x$asset_best_subset)
  }, integer(1))
  # breadth error vs the truth signal each estimate represents:
  # lead is matched to the primary (order-1) truth signal
  lead_err <- lead_w - max(t_breadth)
  res_err <- if (length(res_w) > 0L) {
    res_w - t_breadth[seq_along(res_w)]
  } else {
    integer()
  }
  a3_rows[[length(a3_rows) + 1L]] <- data.frame(
    rep_id = o$rep_id,
    r2_causal_pair = if (length(cm) >= 2L) r2_pairs[1, 2] else NA_real_,
    lead_snp = lead_snp,
    r2_lead_to_causal1 = unname(r2_to_causal[1]),
    r2_lead_to_causal2 = if (length(cm) >= 2L) unname(r2_to_causal[2])
                         else NA_real_,
    lead_breadth = lead_w,
    truth_breadths = paste(t_breadth, collapse = ","),
    resolved_breadths = paste(res_w, collapse = ","),
    lead_breadth_error = lead_err,
    resolved_breadth_error = paste(res_err, collapse = ","),
    paired_improvement = if (length(res_err) > 0L) {
      abs(lead_err) - mean(abs(res_err))
    } else {
      NA_real_
    },
    lead_subset = if (!is.null(lead_e)) {
      paste(lead_e$asset_best_subset, collapse = ",")
    } else {
      ""
    },
    truth_sets = paste(vapply(truth$candidate_traits, paste, character(1),
                              collapse = "+"), collapse = " / "),
    stringsAsFactors = FALSE
  )
}
a3 <- do.call(rbind, a3_rows)
print(a3, row.names = FALSE)
cat("\nsummary:\n")
cat("empirical r2 between causal signals: median",
    round(median(a3$r2_causal_pair, na.rm = TRUE), 3), "\n")
cat("P(lead breadth inflated):",
    mean(a3$lead_breadth_error > 0), "\n")
cat("mean lead breadth error:", round(mean(a3$lead_breadth_error), 3), "\n")
cat("mean paired improvement (|lead err| - |resolved err|):",
    round(mean(a3$paired_improvement, na.rm = TRUE), 3), "\n")

## ---------------------------------------------------------------- B4 ----
cat("\n=== B4 Simulation II-A per-trait audit ===\n")
files2a <- manifest$file[manifest$scenario %in%
                           c("trait_specific", "two_trait_concordant",
                             "two_trait_antagonistic", "broad_concordant")]
b4 <- list()
for (f in files2a) {
  o <- readRDS(f)
  sc <- sub(".*scenario=([^|]*).*", "\\1", o$scenario_id)
  tt <- o$candidate_trait_estimates$trait_table
  cs <- o$estimates$candidate_sets
  for (sid in names(o$truth$candidate_traits)) {
    true_traits <- o$truth$candidate_traits[[sid]]
    est_traits <- if (!is.null(cs[[sid]])) cs[[sid]] else character()
    all_traits <- unique(c(true_traits, est_traits,
                           if (!is.null(tt)) tt$trait else character()))
    for (tr in all_traits) {
      padj <- if (!is.null(tt) && "p_adjusted" %in% names(tt)) {
        tt$p_adjusted[tt$signal_id == sid & tt$trait == tr][1]
      } else {
        NA_real_
      }
      b4[[length(b4) + 1L]] <- data.frame(
        file = basename(f), signal_id = sid,
        scenario = sc, trait = tr,
        is_true = tr %in% true_traits,
        selected = tr %in% est_traits,
        p_adjusted = padj, stringsAsFactors = FALSE
      )
    }
  }
}
b4 <- do.call(rbind, b4)
for (sc in unique(b4$scenario)) {
  sub <- b4[b4$scenario == sc, ]
  cat("\n--", sc, "--\n")
  per_trait <- do.call(rbind, lapply(split(sub, sub$trait), function(x) {
    data.frame(
      trait = x$trait[1],
      TPR = mean(x$selected[x$is_true]),
      FDR_side = mean(x$selected[!x$is_true]),
      med_padj_true = stats::median(x$p_adjusted[x$is_true], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  print(per_trait, row.names = FALSE)
  # missed true traits per (replicate, signal)
  bysig <- split(sub, interaction(sub$file, sub$signal_id, drop = TRUE))
  missed <- vapply(bysig, function(x) sum(x$is_true & !x$selected),
                   integer(1))
  cat("mean missed true traits per signal:", round(mean(missed), 3),
      " overall TPR:", round(mean(sub$selected[sub$is_true]), 3),
      " per-trait false-selection rate:",
      round(mean(sub$selected[!sub$is_true]), 3), "\n")
}

cat("\nPHASE A DIAGNOSTICS DONE\n")
