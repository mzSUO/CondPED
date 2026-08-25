## Stage 7.3: extra-signal source decomposition for mixed_multisignal.
## Adapted from diagnose-stage73-mixed-extra.R to use the stage73 output.
devtools::load_all(quiet = TRUE)

out_dir <- "inst/validation/output/stage73/mixed_extra"
manifest <- readRDS(file.path(out_dir, "manifest.rds"))

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
  set.seed(o$seed)
  do.call(simulate_condped_data, args)
}

overlaps <- function(chr1, s1, e1, chr2, s2, e2) {
  !is.na(chr1) && !is.na(chr2) && chr1 == chr2 && s1 <= e2 && e1 >= s2
}

files <- manifest$file[manifest$status == "ok"]
cats <- character()
rep_has_extra <- logical()
for (f in files) {
  o <- readRDS(f)
  G <- regen_G(o)$G
  mm <- CondPED:::.match_signals(o$truth, o$estimates, G = G)
  extra_ids <- mm$extra
  rep_has_extra <- c(rep_has_extra, length(extra_ids) > 0L)
  if (length(extra_ids) == 0L) next
  est_sig <- o$estimates$signals
  est_loci <- o$loci
  truth_loci <- as.data.frame(o$truth$loci)
  causals <- o$truth$causal_markers
  matched_est_loci <- unique(mm$matches$est_locus_id)
  matched_truth_loci <- unique(mm$matches$truth_locus_id)
  for (eid in extra_ids) {
    es <- est_sig[est_sig$signal_id == eid, ]
    el <- est_loci[est_loci$locus_id == es$locus_id, ]
    r2_max <- max(stats::cor(G[, es$representative_snp],
                             G[, causals])^2)
    ov_truth <- vapply(seq_len(nrow(truth_loci)), function(ti) {
      overlaps(el$chromosome, el$start, el$end,
               truth_loci$chromosome[ti], truth_loci$start[ti],
               truth_loci$end[ti])
    }, logical(1))
    if (es$locus_id %in% matched_est_loci) {
      cats <- c(cats, "within_locus_fp")
    } else if (any(ov_truth)) {
      tid <- truth_loci$locus_id[ov_truth][1]
      if (tid %in% matched_truth_loci && r2_max >= 0.5) {
        cats <- c(cats, "split_duplicate_proxy")
      } else {
        cats <- c(cats, "split_region_fp")
      }
    } else {
      cats <- c(cats, "novel_locus_fp")
    }
  }
}

cat("\n== mixed_multisignal (lp0.02/sp0.0075) ==\n")
cat("total replicates:", length(files), "\n")
cat("replicates with extra signals:", sum(rep_has_extra), "\n")
cat("extra signal rate (replicate level):",
    round(mean(rep_has_extra), 3), "\n")
cat("total extra signals:", length(cats), "\n")
if (length(cats) > 0L) {
  cat("\nproportions:\n")
  print(round(sort(prop.table(table(cats)), decreasing = TRUE), 3))
  cat("\ncounts:\n")
  print(table(cats))
}
cat("\nMIXED EXTRA DECOMPOSITION DONE\n")
