## Stage 7.4 G: marker/locus/signal-level FDP decomposition for
## mixed_multisignal after the same-locus causal geometry fix.
## Per replicate: full pipeline (scan -> BH -> loci -> resolve), then
##   marker level : each BH-selected marker classified as
##     truth_region (inside truth locus interval) /
##     causal_proxy (outside, but r2 >= 0.5 to a causal) /
##     null_marker (outside, low LD)
##   locus level  : estimated locus overlaps the truth locus or not
##   signal level : extra (unmatched) signals decomposed as
##     within_locus_fp / split_duplicate_proxy / split_region_fp /
##     novel_locus_fp (Stage 7.3 A categories)
## No CondPED statistic is modified.
devtools::load_all(quiet = TRUE)

master_seed <- 20260816L
reps <- 100L
workers <- as.integer(Sys.getenv("STAGE74_WORKERS", "16"))
secondary_signal_pve <- as.numeric(Sys.getenv("STAGE74_SPVE", "0.0075"))
out_dir <- "inst/validation/output/stage74/fdp_mixed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenario <- "mixed_multisignal"
experiment <- "end_to_end"
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
      locus_pve = locus_pve, secondary_signal_pve = secondary_signal_pve,
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

    fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
    if (!isTRUE(fit$status$ok)) {
      return(list(rep_id = rep_id, status = "fit_failed"))
    }
    scan <- scan_mt_omnibus(fit, G)
    om <- scan$omnibus
    valid <- !is.na(om$p_value)
    p_adj <- rep(NA_real_, nrow(om))
    p_adj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
    selected <- om$marker_id[valid][p_adj[valid] <= alpha_omnibus]

    tl <- as.data.frame(truth$loci)[1, ]
    in_region <- function(mid) {
      position[mid] >= tl$start & position[mid] <= tl$end
    }
    r2_to_causals <- function(mid) {
      max(stats::cor(G[, mid], G[, causals])^2)
    }

    ## ---- marker level ----
    marker_cat <- character()
    if (length(selected) > 0L) {
      marker_cat <- vapply(selected, function(mid) {
        if (in_region(mid)) "truth_region"
        else if (r2_to_causals(mid) >= 0.5) "causal_proxy"
        else "null_marker"
      }, character(1))
    }

    ## ---- locus / signal level ----
    locus_cat <- character()
    sig_cat <- character()
    n_signals <- 0L
    if (length(selected) > 0L) {
      loci_obj <- define_associated_loci(
        scan, G, chromosome = rep("chr1", ncol(G)),
        position = unname(position),
        selected_markers = selected,
        marker_ids = marker_ids,
        method = "physical", window_bp = window_bp
      )
      if (nrow(loci_obj$loci) > 0L) {
        locus_cat <- vapply(seq_len(nrow(loci_obj$loci)), function(i) {
          ov <- loci_obj$loci$start[i] <= tl$end &&
            loci_obj$loci$end[i] >= tl$start
          if (ov) "truth_overlap" else "novel_locus"
        }, character(1))

        resolution <- resolve_locus_signals(
          fit, G, loci_obj, marker_ids = marker_ids,
          signal_adjust = "within_locus_bonferroni",
          alpha_signal = alpha_signal
        )
        est_signals <- resolution$signals
        n_signals <- nrow(est_signals)
        estimates <- list(signals = est_signals, positions = position)
        mm <- CondPED:::.match_signals(truth, estimates, G = G)
        extra_ids <- mm$extra
        if (length(extra_ids) > 0L) {
          matched_est_loci <- unique(mm$matches$est_locus_id)
          matched_truth_loci <- unique(mm$matches$truth_locus_id)
          sig_cat <- vapply(extra_ids, function(eid) {
            es <- est_signals[est_signals$signal_id == eid, ]
            el <- loci_obj$loci[loci_obj$loci$locus_id == es$locus_id, ]
            r2_max <- r2_to_causals(es$representative_snp)
            ov_truth <- el$start <= tl$end && el$end >= tl$start
            if (es$locus_id %in% matched_est_loci) {
              "within_locus_fp"
            } else if (ov_truth) {
              if (tl$locus_id %in% matched_truth_loci && r2_max >= 0.5) {
                "split_duplicate_proxy"
              } else {
                "split_region_fp"
              }
            } else {
              "novel_locus_fp"
            }
          }, character(1))
        }
      }
    }

    list(
      rep_id = rep_id, status = "ok",
      n_selected = length(selected),
      marker_cat = marker_cat,
      locus_cat = locus_cat,
      sig_cat = sig_cat,
      n_signals = n_signals
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
cat("WALL_SECONDS:", round(difftime(Sys.time(), t0, units = "secs"), 1), "\n")

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
saveRDS(out, file.path(out_dir, "fdp_rows.rds"))

rows <- out[ok]
mc <- table(unlist(lapply(rows, function(x) x$marker_cat)))
lc <- table(unlist(lapply(rows, function(x) x$locus_cat)))
sc <- table(unlist(lapply(rows, function(x) x$sig_cat)))
n_sel <- vapply(rows, function(x) x$n_selected, numeric(1))

cat("\n== marker level ==\n")
cat("selected markers total:", sum(n_sel),
    " mean per rep:", round(mean(n_sel), 2), "\n")
print(mc)
if (length(mc)) print(round(prop.table(mc), 3))
cat("marker FDP (null_marker / selected):",
    round(if ("null_marker" %in% names(mc)) mc[["null_marker"]] / sum(mc) else 0, 3), "\n")

cat("\n== locus level ==\n")
print(lc)
if (length(lc)) print(round(prop.table(lc), 3))
cat("locus FDP (novel_locus / loci):",
    round(if ("novel_locus" %in% names(lc)) lc[["novel_locus"]] / sum(lc) else 0, 3), "\n")

cat("\n== signal level (extra signals only) ==\n")
print(sc)
if (length(sc)) print(round(prop.table(sc), 3))
n_sig <- vapply(rows, function(x) x$n_signals, numeric(1))
n_extra <- vapply(rows, function(x) length(x$sig_cat), numeric(1))
cat("replicates with extra signals:", round(mean(n_extra > 0), 3), "\n")
cat("signal FDP (extra / estimated):",
    round(sum(n_extra) / max(sum(n_sig), 1), 3), "\n")
cat("\nSTAGE74 FDP DECOMPOSITION DONE\n")
