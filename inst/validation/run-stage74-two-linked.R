## Stage 7.4 B+D: three-layer paired diagnostic for
## two_linked_trait_specific after the same-locus causal geometry fix
## (Stage 7.4 A), with the corrected within-locus secondary recovery
## metric (Stage 7.4 B):
##   secondary_recovered = the PRIMARY and the SECONDARY truth signal are
##   each matched by a distinct estimated signal (with 1-1 greedy matching
##   this is n_matched == 2). The legacy metric (any matched truth signal
##   other than signal 1) is still reported for comparison.
## Layers: 3A full pipeline, 3B oracle-locus, 3C oracle-primary.
## No CondPED statistic is modified.
devtools::load_all(quiet = TRUE)

master_seed <- 20260815L   # same master seed as Stage 7.3 two_linked
reps <- 100L
workers <- as.integer(Sys.getenv("STAGE74_WORKERS", "16"))
secondary_signal_pve <- as.numeric(Sys.getenv("STAGE74_SPVE", "0.0075"))
out_dir <- "inst/validation/output/stage74/two_linked"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenario <- "two_linked_trait_specific"
experiment <- "signal_resolution"
locus_pve <- 0.02
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

  ## corrected metric helper: both truth signals matched (1-1 greedy =>
  ## distinct estimated signals); legacy metric kept alongside.
  recovery <- function(mm) {
    matched_truth <- unique(mm$matches$truth_signal_id)
    list(
      n_matched = nrow(mm$matches),
      n_extra = length(mm$extra),
      secondary_recovered = nrow(mm$matches) >= n_truth &&
        all(truth_signals$signal_id %in% matched_truth),
      secondary_recovered_legacy = any(matched_truth !=
                                         truth_signals$signal_id[1])
    )
  }

  ## ---------- 3A: full pipeline ------------------------------------------
  full <- list(locus_detected = FALSE, n_loci = 0L,
               truth_locus_split = FALSE, n_signals = 0L,
               signal_count_exact = FALSE, secondary_recovered = FALSE,
               secondary_recovered_legacy = FALSE,
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

    estimates <- list(signals = est_signals, positions = position)
    mm <- CondPED:::.match_signals(truth, estimates, G = G)
    rec <- recovery(mm)
    full$n_matched <- rec$n_matched
    full$n_extra <- rec$n_extra
    full$secondary_recovered <- rec$secondary_recovered
    full$secondary_recovered_legacy <- rec$secondary_recovered_legacy

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
  estimates_o <- list(signals = oracle_est, positions = position)
  mm_o <- CondPED:::.match_signals(truth, estimates_o, G = G)
  rec_o <- recovery(mm_o)
  oracle <- list(
    n_signals = nrow(oracle_est),
    signal_count_exact = nrow(oracle_est) == n_truth,
    secondary_recovered = rec_o$secondary_recovered,
    secondary_recovered_legacy = rec_o$secondary_recovered_legacy,
    n_extra = rec_o$n_extra, n_matched = rec_o$n_matched
  )

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
  M_remaining <- length(oracle_members) - 1L
  p_adj_sec <- min(1, M_remaining * ctab$p_value)
  oracle_primary <- list(
    Q = ctab$Q, p_raw = ctab$p_value, p_adj = p_adj_sec,
    pass = p_adj_sec <= alpha_signal,
    df = ctab$df, status = ctab$status
  )

    list(
      rep_id = rep_id, seed = seed, status = "ok",
      empirical_r2_causal = stats::cor(G[, causals[1]], G[, causals[2]])^2,
      causal_dist = abs(diff(unname(position[causals]))),
      n_sig_genomewide = length(selected_markers),
      full = full, oracle = oracle, oracle_primary = oracle_primary
    )
  }, error = function(e) {
    list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
  })
}

rds <- file.path(out_dir, sprintf("paired_rows_sp%.4f.rds",
                                  secondary_signal_pve))
t0 <- Sys.time()
if (file.exists(rds)) {
  cat("SKIP: existing", rds, "\n")
  out <- readRDS(rds)
} else {
  if (.Platform$OS.type == "unix") {
    out <- parallel::mclapply(seq_len(reps), one_rep, mc.cores = workers)
  } else {
    out <- lapply(seq_len(reps), one_rep)
  }
  saveRDS(out, rds)
}
wall <- difftime(Sys.time(), t0, units = "secs")
cat("WALL_SECONDS:", round(wall, 1), "\n")

is_ok <- function(x) {
  is.list(x) && !is.null(x$status) && is.character(x$status) &&
    length(x$status) == 1L && x$status == "ok"
}
ok <- vapply(out, is_ok, logical(1))
cat("OK replicates:", sum(ok), "/", length(out), "\n")
if (any(!ok)) {
  print(table(vapply(out[!ok], function(x) as.character(x$status),
                     character(1))))
}

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
    "secondary_recovery_full_legacy",
    "secondary_recovery_oracle_locus",
    "secondary_recovery_oracle_locus_legacy",
    "secondary_conditional_power_oracle_primary",
    "extra_signal_rate_full",
    "extra_signal_rate_oracle_locus",
    "median_causal_dist_bp",
    "mean_empirical_r2_causal"
  ),
  value = c(
    mean(fn(function(x) x$full$locus_detected)),
    mean(fn(function(x) x$full$truth_locus_split)),
    mean(fn(function(x) x$full$signal_count_exact)),
    mean(fn(function(x) x$oracle$signal_count_exact)),
    mean(fn(function(x) x$full$secondary_recovered)),
    mean(fn(function(x) x$full$secondary_recovered_legacy)),
    mean(fn(function(x) x$oracle$secondary_recovered)),
    mean(fn(function(x) x$oracle$secondary_recovered_legacy)),
    mean(fn(function(x) x$oracle_primary$pass)),
    mean(nn(function(x) x$full$n_extra) > 0),
    mean(nn(function(x) x$oracle$n_extra) > 0),
    stats::median(nn(function(x) x$causal_dist)),
    mean(nn(function(x) x$empirical_r2_causal))
  )
)
print(summ, row.names = FALSE)

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

tag <- sprintf("%.4f", secondary_signal_pve)
write.csv(summ, file.path(out_dir, sprintf("summary_sp%s.csv", tag)),
          row.names = FALSE)
write.csv(trans, file.path(out_dir, sprintf("transitions_sp%s.csv", tag)),
          row.names = FALSE)
cat("STAGE74 TWO-LINKED PAIRED DIAGNOSTIC DONE\n")
