## Stage 7.3: three-layer paired diagnostic for two_linked_trait_specific.
## One simulated dataset per replicate feeds three analysis layers:
##   3A full       : scan -> BH -> define_associated_loci -> resolve -> joint
##   3B oracle-loc : same dataset, truth-defined locus, same resolver
##   3C oracle-pri : same dataset, condition on true primary causal SNP,
##                   test true secondary causal SNP with the conditional scan
## No CondPED statistic is modified.
devtools::load_all(quiet = TRUE)

master_seed <- 20260815L
reps <- 100L
workers <- 4L
out_dir <- "inst/validation/output/stage73/two_linked"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenario <- "two_linked_trait_specific"
experiment <- "signal_resolution"
locus_pve <- 0.02
secondary_signal_pve <- 0.0075
alpha_omnibus <- 0.05
alpha_signal <- 0.05
window_bp <- 5000L

one_rep <- function(rep_id) {
  tryCatch({
    psid <- sprintf("%s|lp%.3f_sp%.4f", scenario, locus_pve,
                    secondary_signal_pve)
    seed <- CondPED:::.seed_for_rep(master_seed, psid, rep_id)
    set.seed(seed)

    sim <- simulate_condped_data(
      n = 1000L, m = 4L, p = 1000L,
      experiment = experiment, scenario = scenario,
      locus_pve = locus_pve,
      secondary_signal_pve = secondary_signal_pve,
      correlation = "block", tolerance = 0.10
    )
    if (!isTRUE(sim$status$ok)) {
      return(list(rep_id = rep_id, status = "sim_failed"))
  }
  G <- sim$G
  truth <- sim$truth
  marker_ids <- colnames(G)
  position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
  causals <- truth$causal_markers
  truth_signals <- truth$signals
  n_truth <- nrow(truth_signals)

  fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
  if (!isTRUE(fit$status$ok)) {
    return(list(rep_id = rep_id, status = "fit_failed"))
  }
  scan <- scan_mt_omnibus(fit, G)
  om <- scan$omnibus
  valid <- !is.na(om$p_value)
  p_adj <- rep(NA_real_, nrow(om))
  p_adj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
  om$padj <- p_adj
  selected_markers <- om$marker_id[valid][p_adj[valid] <= alpha_omnibus]

  ## ---------- 3A: full pipeline ------------------------------------------
  full <- list(locus_detected = FALSE, n_loci = 0L,
               truth_locus_split = FALSE, n_signals = 0L,
               signal_count_exact = FALSE, secondary_recovered = FALSE,
               n_extra = 0L, n_matched = 0L,
               loci_detail = NULL, sig_markers = character())
  if (length(selected_markers) > 0L) {
    loci_obj <- define_associated_loci(
      scan, G, chromosome = rep("chr1", ncol(G)),
      position = unname(position),
      selected_markers = selected_markers,
      marker_ids = marker_ids,
      method = "physical", window_bp = window_bp
    )
    full$locus_detected <- nrow(loci_obj$loci) > 0L
    full$n_loci <- nrow(loci_obj$loci)
    full$sig_markers <- selected_markers
    full$loci_detail <- loci_obj$loci

    resolution <- resolve_locus_signals(fit, G, loci_obj,
                                        marker_ids = marker_ids,
                                        signal_adjust = "within_locus_bonferroni",
                                        alpha_signal = alpha_signal)
    est_signals <- resolution$signals
    full$n_signals <- nrow(est_signals)
    full$signal_count_exact <- nrow(est_signals) == n_truth

    ## signal matching
    estimates <- list(signals = est_signals, positions = position)
    mm <- CondPED:::.match_signals(truth, estimates, G = G)
    full$n_matched <- nrow(mm$matches)
    full$n_extra <- length(mm$extra)
    matched_truth_ids <- unique(mm$matches$truth_signal_id)
    full$secondary_recovered <- any(matched_truth_ids != truth_signals$signal_id[1])

    ## truth locus split: does the truth locus map to >1 estimated locus?
    truth_lid <- truth$loci$locus_id[1]
    est_lids_for_truth <- unique(mm$matches$est_locus_id[
      mm$matches$truth_locus_id == truth_lid])
    full$truth_locus_split <- length(est_lids_for_truth) > 1L
  }

  ## ---------- 3B: oracle-locus pipeline -----------------------------------
  tl <- as.data.frame(truth$loci)[1, ]
  oracle_members <- marker_ids[position >= tl$start & position <= tl$end]
  oracle_lead <- oracle_members[which.min(om$p_value[match(oracle_members,
                                                           om$marker_id)])]
  r2_lead <- stats::cor(G[, oracle_members], G[, oracle_lead])^2
  oracle_locus <- list(
    loci = data.frame(
      locus_id = tl$locus_id, chromosome = tl$chromosome,
      start = tl$start, end = tl$end,
      lead_snp = oracle_lead,
      lead_p = om$p_value[match(oracle_lead, om$marker_id)],
      n_significant_markers = length(oracle_members),
      n_region_markers = length(oracle_members), status = "ok",
      stringsAsFactors = FALSE
    ),
    membership = data.frame(
      locus_id = tl$locus_id, marker_id = oracle_members,
      significant_in_marginal_scan = oracle_members %in% selected_markers,
      r2_to_lead = as.numeric(r2_lead),
      position = unname(position[oracle_members]),
      stringsAsFactors = FALSE
    )
  )
  oracle_res <- resolve_locus_signals(fit, G, oracle_locus,
                                      marker_ids = marker_ids,
                                      signal_adjust = "within_locus_bonferroni",
                                      alpha_signal = alpha_signal)
  oracle_est <- oracle_res$signals
  oracle <- list(
    n_signals = nrow(oracle_est),
    signal_count_exact = nrow(oracle_est) == n_truth,
    secondary_recovered = FALSE, n_extra = 0L, n_matched = 0L
  )
  estimates_o <- list(signals = oracle_est, positions = position)
  mm_o <- CondPED:::.match_signals(truth, estimates_o, G = G)
  oracle$n_matched <- nrow(mm_o$matches)
  oracle$n_extra <- length(mm_o$extra)
  matched_truth_o <- unique(mm_o$matches$truth_signal_id)
  oracle$secondary_recovered <- any(matched_truth_o != truth_signals$signal_id[1])

  ## ---------- 3C: oracle-primary diagnostic -------------------------------
  pri_causal <- causals[1]
  sec_causal <- causals[2]
  X_C <- G[, match(pri_causal, marker_ids), drop = FALSE]
  proj <- CondPED:::.build_conditional_projection(fit, X_C)
  cond_scan <- CondPED:::.conditional_mt_scan(
    proj, G[, match(sec_causal, marker_ids), drop = FALSE],
    marker_ids = sec_causal, return_effects = FALSE
  )
  ctab <- cond_scan$conditional
  M_remaining <- length(oracle_members) - 1L   # candidates after primary
  p_adj_sec <- min(1, M_remaining * ctab$p_value)
  oracle_primary <- list(
    Q = ctab$Q, p_raw = ctab$p_value, p_adj = p_adj_sec,
    pass = p_adj_sec <= alpha_signal,
    df = ctab$df, status = ctab$status
  )

  ## ---------- locus-splitting detail (3A split replicates) ---------------
  split_detail <- NULL
  if (isTRUE(full$truth_locus_split) && !is.null(full$loci_detail)) {
    causal_pos <- position[causals]
    sig_in_region <- om[om$marker_id %in% selected_markers &
                        om$marker_id %in% oracle_members, ]
    split_detail <- list(
      causal_ids = causals,
      causal_pos = unname(causal_pos),
      causal_dist = abs(diff(unname(causal_pos))),
      causal_r2 = stats::cor(G[, causals[1]], G[, causals[2]])^2,
      n_sig_in_region = nrow(sig_in_region),
      sig_marker_ids = sig_in_region$marker_id,
      sig_marker_pos = unname(position[sig_in_region$marker_id]),
      sig_marker_padj = sig_in_region$padj,
      est_locus_ids = full$loci_detail$locus_id,
      est_locus_start = full$loci_detail$start,
      est_locus_end = full$loci_detail$end,
      causal1_in_est = full$loci_detail$locus_id[
        full$loci_detail$start <= causal_pos[1] &
          full$loci_detail$end >= causal_pos[1]],
      causal2_in_est = full$loci_detail$locus_id[
        full$loci_detail$start <= causal_pos[2] &
          full$loci_detail$end >= causal_pos[2]]
    )
  }

    list(
      rep_id = rep_id, seed = seed, status = "ok",
      empirical_r2_causal = stats::cor(G[, causals[1]], G[, causals[2]])^2,
      causal_dist = abs(diff(unname(position[causals]))),
      n_sig_genomewide = length(selected_markers),
      full = full, oracle = oracle, oracle_primary = oracle_primary,
      split_detail = split_detail
    )
  }, error = function(e) {
    list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
  })
}

t0 <- Sys.time()
if (.Platform$OS.type == "unix") {
  out <- parallel::mclapply(seq_len(reps), one_rep, mc.cores = workers)
} else {
  out <- lapply(seq_len(reps), one_rep)
}
wall <- difftime(Sys.time(), t0, units = "secs")
cat("WALL_SECONDS:", round(wall, 1), "\n")

## robust ok check: x must be a list with character status == "ok"
is_ok <- function(x) {
  is.list(x) && !is.null(x$status) && is.character(x$status) &&
    length(x$status) == 1L && x$status == "ok"
}
ok <- vapply(out, is_ok, logical(1))
cat("OK replicates:", sum(ok), "/", length(out), "\n")

saveRDS(out, file.path(out_dir, "paired_rows.rds"))

## ---------- summary -------------------------------------------------------
rows <- out[ok]
fn <- function(f) vapply(rows, function(x) f(x), logical(1))
nn <- function(f) vapply(rows, function(x) f(x), numeric(1))

summ <- data.frame(
  metric = c(
    "full_locus_detection_rate",
    "truth_locus_split_rate",
    "signal_count_exact_full",
    "signal_count_exact_oracle_locus",
    "secondary_recovery_full",
    "secondary_recovery_oracle_locus",
    "secondary_conditional_power_oracle_primary",
    "extra_signal_rate_full",
    "extra_signal_rate_oracle_locus"
  ),
  value = c(
    mean(fn(function(x) x$full$locus_detected)),
    mean(fn(function(x) x$full$truth_locus_split)),
    mean(fn(function(x) x$full$signal_count_exact)),
    mean(fn(function(x) x$oracle$signal_count_exact)),
    mean(fn(function(x) x$full$secondary_recovered)),
    mean(fn(function(x) x$oracle$secondary_recovered)),
    mean(fn(function(x) x$oracle_primary$pass)),
    mean(nn(function(x) x$full$n_extra) > 0),
    mean(nn(function(x) x$oracle$n_extra) > 0)
  )
)
print(summ, row.names = FALSE)

## paired transitions
full_s <- fn(function(x) x$full$secondary_recovered)
oracle_s <- fn(function(x) x$oracle$secondary_recovered)
primary_s <- fn(function(x) x$oracle_primary$pass)
trans <- data.frame(
  transition = c("full_fail_oracle_success",
                 "oracle_fail_primary_success",
                 "all_success", "all_fail"),
  count = c(
    sum(!full_s & oracle_s),
    sum(!oracle_s & primary_s),
    sum(full_s & oracle_s & primary_s),
    sum(!full_s & !oracle_s & !primary_s)
  )
)
print(trans, row.names = FALSE)

write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(trans, file.path(out_dir, "transitions.csv"), row.names = FALSE)
cat("TWO-LINKED PAIRED DIAGNOSTIC DONE\n")
