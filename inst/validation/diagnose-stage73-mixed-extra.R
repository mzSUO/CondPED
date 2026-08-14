## Stage 7.3 A: per-replicate source decomposition of extra (unmatched)
## estimated signals in mixed_multisignal, baseline vs lp0.02/sp0.0075.
## Categories:
##   split_duplicate_proxy : extra in an est locus overlapping an
##     already-matched truth locus, representative in high LD (r2>=0.5)
##     with a causal marker (alternative representative produced by
##     locus splitting / 1-1 matching)
##   within_locus_fp       : extra in the SAME est locus as a matched
##     signal, low LD to causals (local conditional false positive)
##   split_region_fp       : extra in a different est locus overlapping
##     a truth locus, low LD to causals
##   novel_locus_fp        : extra in an est locus overlapping no truth
##     locus (genome-wide false discovery)
devtools::load_all(quiet = TRUE)

out_dir <- "inst/validation/output/stage72"
grid <- readRDS(file.path(out_dir, "grid.rds"))
manifest <- readRDS(file.path(out_dir, "manifest.rds"))
grid$scenario_id <- vapply(seq_len(nrow(grid)), function(i) {
  CondPED:::.canonical_scenario_id(as.list(grid[i, setdiff(names(grid), "parameter_set_id")]))
}, character(1))

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

decompose <- function(psid) {
  sid <- grid$scenario_id[grid$parameter_set_id == psid]
  files <- manifest$file[manifest$scenario_id == sid & manifest$status == "ok"]
  cats <- character()
  for (f in files) {
    o <- readRDS(f)
    G <- regen_G(o)$G
    mm <- CondPED:::.match_signals(o$truth, o$estimates, G = G)
    extra_ids <- mm$extra
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
  cats
}

for (psid in c("mixed|baseline", "mixed|lp0.02_sp0.0075")) {
  cats <- decompose(psid)
  cat("\n==", psid, "==  total extra signals:", length(cats), "\n")
  print(round(sort(prop.table(table(cats)), decreasing = TRUE), 3))
  print(table(cats))
}
cat("\nA DONE\n")
