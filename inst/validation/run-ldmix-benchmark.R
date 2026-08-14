## Stage 7.3 B: controlled LD-mixing mechanistic benchmark.
## Truth is always: signal1 -> Trait1, signal2 -> Trait2, two LD-linked
## signals inside ONE predefined locus. No genome-wide discovery decides
## entry: locus membership is fixed from the simulation truth. The same
## Y/G feeds both comparators (marginal lead + ASSET vs conditional
## resolved signals + ASSET). No CondPED statistic is modified.
##
## This benchmark measures the LD-mixing phenomenon and the contribution
## of SNP-space conditional resolution ONLY; it is not a genome-wide
## discovery performance statement.
devtools::load_all(quiet = TRUE)

master_seed <- 20260815L
reps <- 30L
out_dir <- "inst/validation/output/ldmix-benchmark"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

combos <- expand.grid(target_r2 = c(0.3, 0.5, 0.7),
                      scale = c(1.0, 1.5, 2.0),
                      ratio = c(0.75, 1.0))
combos$locus_pve <- 0.01 * combos$scale
combos$secondary_signal_pve <- combos$locus_pve * combos$ratio^2
combos$psid <- sprintf("tr2%.1f|s%.1f|r%.2f", combos$target_r2,
                       combos$scale, combos$ratio)

one_rep <- function(combo, rep_id) {
  psid <- combo$psid
  seed <- CondPED:::.seed_for_rep(master_seed, psid, rep_id)
  set.seed(seed)
  sim <- simulate_condped_data(
    n = 1000, m = 4, p = 1000,
    experiment = "end_to_end", scenario = "linked_pseudo_multitrait",
    locus_pve = combo$locus_pve,
    secondary_signal_pve = combo$secondary_signal_pve,
    target_r2 = combo$target_r2,
    correlation = "block", tolerance = 0.10
  )
  if (!isTRUE(sim$status$ok)) {
    return(data.frame(psid = psid, rep_id = rep_id, status = "sim_failed"))
  }
  G <- sim$G
  truth <- sim$truth
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
  scan <- scan_mt_omnibus(fit, G)
  om <- scan$omnibus

  # predefined locus membership from truth (bypasses discovery)
  position <- stats::setNames(seq_len(ncol(G)) * 1000, colnames(G))
  tl <- as.data.frame(truth$loci)[1, ]
  members <- colnames(G)[position >= tl$start & position <= tl$end]
  lead <- members[which.min(om$p_value[match(members, om$marker_id)])]
  r2_lead <- stats::cor(G[, members], G[, lead])^2
  locus_object <- list(
    loci = data.frame(
      locus_id = tl$locus_id, chromosome = tl$chromosome,
      start = tl$start, end = tl$end,
      lead_snp = lead,
      lead_p = om$p_value[match(lead, om$marker_id)],
      n_significant_markers = length(members),
      n_region_markers = length(members), status = "ok",
      stringsAsFactors = FALSE
    ),
    membership = data.frame(
      locus_id = tl$locus_id, marker_id = members,
      significant_in_marginal_scan = NA,
      r2_to_lead = as.numeric(r2_lead),
      position = unname(position[members]),
      stringsAsFactors = FALSE
    )
  )

  pipes <- CondPED:::.run_comparison_pipelines(
    fit, G, locus_object, scan,
    trait_names = fit$trait_names, tolerance = 0.10
  )

  lead_e <- pipes$lead_asset$loci[[tl$locus_id]]
  lead_sub <- if (is.null(lead_e$asset_best_subset)) character() else
    lead_e$asset_best_subset
  res <- pipes$resolved_asset$signals
  truth_sets <- truth$candidate_traits          # named by truth signal id
  causals <- truth$causal_markers

  # match each resolved signal to the truth signal whose causal marker
  # has max r2 with the resolved representative
  res_info <- lapply(res, function(x) {
    r2c <- stats::cor(G[, x$representative_snp], G[, causals])^2
    best <- which.max(r2c)
    list(subset = if (is.null(x$asset_best_subset)) character() else
      x$asset_best_subset,
      truth_set = truth_sets[[best]],
      r2_causal = as.numeric(r2c[best]))
  })
  res_exact <- vapply(res_info, function(x) {
    setequal(x$subset, x$truth_set)
  }, logical(1))
  res_err <- vapply(res_info, function(x) {
    length(x$subset) - length(x$truth_set)
  }, numeric(1))
  lead_err <- length(lead_sub) -
    length(truth_sets[[1]])   # lead vs primary (order-1) signal
  t_breadth <- vapply(truth_sets, length, integer(1))

  data.frame(
    psid = psid, rep_id = rep_id, status = "ok",
    empirical_r2_causal = if (length(causals) >= 2L) {
      stats::cor(G[, causals[1]], G[, causals[2]])^2
    } else {
      NA_real_
    },
    lead_has_both_causal_traits = all(unique(unlist(truth_sets)) %in%
                                        lead_sub),
    lead_breadth = length(lead_sub),
    lead_breadth_error = lead_err,
    n_resolved = length(res),
    signal_count_recovered = length(res) == length(truth_sets),
    res_exact_trait_rate = if (length(res_exact)) mean(res_exact) else NA_real_,
    res_mean_breadth_error = if (length(res_err)) mean(res_err) else NA_real_,
    paired_improvement = if (length(res_err)) {
      abs(lead_err) - mean(abs(res_err))
    } else {
      NA_real_
    },
    median_rep_to_causal_r2 = if (length(res_info)) {
      stats::median(vapply(res_info, `[[`, numeric(1), "r2_causal"))
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

tasks <- expand.grid(ci = seq_len(nrow(combos)), rep_id = seq_len(reps))
run_task <- function(ti) {
  combo <- combos[tasks$ci[ti], ]
  tryCatch(one_rep(combo, tasks$rep_id[ti]),
           error = function(e) {
             data.frame(psid = combo$psid, rep_id = tasks$rep_id[ti],
                        status = paste("error:", conditionMessage(e)))
           })
}
t0 <- Sys.time()
if (.Platform$OS.type == "unix") {
  out <- parallel::mclapply(seq_len(nrow(tasks)), run_task,
                            mc.cores = 4L)
} else {
  out <- lapply(seq_len(nrow(tasks)), run_task)
}
cat("WALL_SECONDS:", round(difftime(Sys.time(), t0, units = "secs"), 1), "\n")
res <- do.call(rbind, out)
saveRDS(res, file.path(out_dir, "benchmark_rows.rds"))
print(table(res$status))
cat("BENCHMARK RUN DONE\n")
