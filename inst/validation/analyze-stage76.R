## Stage 7.6: method-vs-simulation-vs-evaluator attribution audit.
## Pure diagnostic. No CondPED statistic, threshold, locus definition,
## matching, or frozen simulation parameter is modified. Deterministic
## seeds throughout; two_linked Part 1 reuses the Stage 7.4/7.5 seed so
## all layers are strictly paired with earlier stages.
##
## Parts (cached per part; rerun skips completed parts):
##   P0 problem attribution table (static evidence synthesis)
##   P1 oracle decomposition + geometry sensitivity (two_linked 0.0175)
##   P2 FDP null/single separation (+ reuse stage75 mixed for multi)
##   P3 LD-mixing gradient r2 in {0, 0.3, 0.6, 0.9}
##   P4 attribution rho/tau boundary pilot
##   P5 PVE ratio descriptive + final recommendation
devtools::load_all(quiet = TRUE)

workers <- as.integer(Sys.getenv("STAGE76_WORKERS", "4"))
out_dir <- "inst/validation/output/stage76"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

alpha_omnibus <- 0.05
alpha_signal <- 0.05
alpha_trait <- 0.05
window_bp <- 5000L
tolerance <- 0.10

is_ok <- function(x) {
  is.list(x) && !is.null(x$status) && is.character(x$status) &&
    length(x$status) == 1L && x$status == "ok"
}
plapply <- function(X, FUN, ...) {
  if (.Platform$OS.type == "unix") {
    parallel::mclapply(X, FUN, ..., mc.cores = workers)
  } else {
    lapply(X, FUN, ...)
  }
}

## ==========================================================================
## P0: problem attribution table (evidence synthesis from Stages 7.2-7.5)
## ==========================================================================
p0_csv <- file.path(out_dir, "stage76_problem_attribution.csv")
if (!file.exists(p0_csv)) {
  p0 <- data.frame(
    problem = c(
      "single_multi_trait Rep ~0.54 after candidate exact",
      "truth_locus_split in two-signal scenarios",
      "causal distance 49kb vs window 5kb",
      "cluster fix split 0.21->0.02",
      "secondary recovery metric n_matched=1 flaw",
      "secondary PVE 0.0075->0.0175, power 0.26->0.85",
      "full count_exact 0.37 vs oracle-locus 0.74",
      "full->oracle loss from novel-locus extras",
      "oracle-locus 23/100 failures",
      "marker FDP 0.047 vs locus 0.579 / signal 0.503",
      "92.7% extras from novel-locus FP",
      "ASSET P_valid_mixing 0.22",
      "conditional contamination reduction ~54%",
      "Stage 7.3 ASSET missing",
      "Stage 7.3.5 rank-deficient on zero-information marker",
      "REML BFGS overflow failures",
      "oracle-locus not a strict upper bound (M=49)"
    ),
    attribution = c(
      "B/C: truth rho 0.0588 vs tau 0.10 boundary; metric interpretation",
      "B: simulation geometry (resolved by Stage 7.4 A)",
      "B: simulation geometry bug (49kb ends of region)",
      "B: confirms geometry was the cause",
      "C: evaluator metric flaw (fixed in Stage 7.4 B)",
      "B: under-powered simulation parameter (calibrated)",
      "C+E-metric: count_exact conflates recovery and FDR control",
      "C: discovery-layer BH FPs propagate 1:1 into signals",
      "mixed: 12 screening/multiplicity (B-eval design), 11 power (B)",
      "C: marker BH does not imply locus/signal FDR control",
      "C: same propagation chain",
      "mechanism confirmed (real LD mixing); environment fixed",
      "method works as designed (conditional decomposition)",
      "E: dependency missing (resolved)",
      "D: numerical (fixed, scale-aware zero floor)",
      "D: numerical (fixed, finiteness guard)",
      "C: oracle benchmark design (candidate-set size mismatch)"
    ),
    change_production = c(rep("no", 13), "n/a", "fixed", "fixed", "no"),
    change_simulation = c(
      "no (documented)", "done (7.4 A)", "done (7.4 A)", "done (7.4 A)",
      "no", "done (7.4 C)", "no", "no", "redefine oracle (7.6)",
      "no", "no", "no", "no", "n/a", "no", "no", "redefine oracle (7.6)"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(p0, p0_csv, row.names = FALSE)
  cat("P0 done\n")
}

## ==========================================================================
## P1: oracle decomposition + geometry sensitivity (two_linked @ 0.0175)
## ==========================================================================
p1_rds <- file.path(out_dir, "p1_rows.rds")
if (!file.exists(p1_rds)) {
  p1_rep <- function(rep_id) {
    tryCatch({
      scenario <- "two_linked_trait_specific"
      psid <- sprintf("%s|lp%.3f_sp%.4f", scenario, 0.02, 0.0175)
      seed <- CondPED:::.seed_for_rep(20260815L, psid, rep_id)
      set.seed(seed)
      sim <- simulate_condped_data(
        n = 1000L, m = 4L, p = 1000L,
        experiment = "signal_resolution", scenario = scenario,
        locus_pve = 0.02, secondary_signal_pve = 0.0175,
        correlation = "block", tolerance = tolerance
      )
      if (!isTRUE(sim$status$ok)) {
        return(list(rep_id = rep_id, status = "sim_failed"))
      }
      G <- sim$G
      truth <- sim$truth
      marker_ids <- colnames(G)
      pos0 <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
      causals <- truth$causal_markers
      truth_signals <- truth$signals
      pri_id <- truth_signals$signal_id[1]
      n_truth <- nrow(truth_signals)

      fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
      if (!isTRUE(fit$status$ok)) {
        return(list(rep_id = rep_id, status = "fit_failed"))
      }
      scan <- scan_mt_omnibus(fit, G)
      om <- scan$omnibus
      valid <- !is.na(om$p_value)
      padj <- rep(NA_real_, nrow(om))
      padj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
      selected <- om$marker_id[valid][padj[valid] <= alpha_omnibus]

      tl <- as.data.frame(truth$loci)[1, ]
      in_region <- pos0 >= tl$start & pos0 <= tl$end
      region_members <- marker_ids[in_region]

      resolve_match <- function(lo) {
        res <- resolve_locus_signals(fit, G, lo, marker_ids = marker_ids,
                                     signal_adjust = "within_locus_bonferroni",
                                     alpha_signal = alpha_signal)
        est <- res$signals
        mm <- CondPED:::.match_signals(
          truth, list(signals = est, positions = pos0), G = G)
        matched_truth <- unique(mm$matches$truth_signal_id)
        list(
          n_signals = nrow(est), n_matched = nrow(mm$matches),
          n_extra = length(mm$extra),
          P1 = pri_id %in% matched_truth,
          S1 = nrow(mm$matches) >= n_truth &&
            all(truth_signals$signal_id %in% matched_truth),
          count_exact = nrow(est) == n_truth
        )
      }
      mk_locus <- function(members, positions) {
        pos_v <- unname(positions[members])
        # locus id must stay in the parseable chr:start-end form so that
        # .match_signals can derive the interval (frozen matching rule)
        lid <- sprintf("chr1:%d-%d", min(pos_v), max(pos_v))
        lead <- members[which.min(om$p_value[match(members, om$marker_id)])]
        list(
          loci = data.frame(
            locus_id = lid, chromosome = "chr1", start = min(pos_v),
            end = max(pos_v), lead_snp = lead,
            lead_p = om$p_value[match(lead, om$marker_id)],
            n_significant_markers = length(members),
            n_region_markers = length(members), status = "ok",
            stringsAsFactors = FALSE
          ),
          membership = data.frame(
            locus_id = lid, marker_id = members,
            significant_in_marginal_scan = members %in% selected,
            r2_to_lead = as.numeric(stats::cor(G[, members], G[, lead])^2),
            position = pos_v, stringsAsFactors = FALSE
          )
        )
      }

      ## Oracle 1: oracle-primary (true primary conditioned, test secondary)
      X_C <- G[, causals[1], drop = FALSE]
      proj <- CondPED:::.build_conditional_projection(fit, X_C)
      cs <- CondPED:::.conditional_mt_scan(
        proj, G[, causals[2], drop = FALSE],
        marker_ids = causals[2], return_effects = FALSE
      )
      M_region <- length(region_members) - 1L
      op_pass <- isTRUE(cs$conditional$p_value <= alpha_signal / M_region)
      op_pass_M1 <- isTRUE(cs$conditional$p_value <= alpha_signal)

      ## Oracle 2: oracle-causal-set (only the 2 true causals as candidates)
      oc <- resolve_match(mk_locus(causals, pos0))

      ## Oracle 3: oracle-locus (all 50 region markers)
      ol <- resolve_match(mk_locus(region_members, pos0))

      ## full pipeline
      full <- list(n_signals = 0L, n_matched = 0L, n_extra = 0L,
                   P1 = FALSE, S1 = FALSE, count_exact = FALSE)
      if (length(selected) > 0L) {
        loci_obj <- define_associated_loci(
          scan, G, chromosome = rep("chr1", ncol(G)),
          position = unname(pos0), selected_markers = selected,
          marker_ids = marker_ids, method = "physical",
          window_bp = window_bp
        )
        if (nrow(loci_obj$loci) > 0L) {
          full <- resolve_match(loci_obj)
          ## full true-locus candidate size (fairness comparison vs 50)
          memb <- loci_obj$membership
          tl_lids <- loci_obj$loci$locus_id[
            loci_obj$loci$start <= tl$end & loci_obj$loci$end >= tl$start]
          full$true_locus_n_candidates <- sum(memb$locus_id %in% tl_lids)
        }
      }
      if (is.null(full$true_locus_n_candidates)) {
        full$true_locus_n_candidates <- 0L
      }

      ## geometry sensitivity: emulate causal spacing D kb by relabeling
      ## positions AFTER causal 1 (genotypes and LD unchanged; only window
      ## grouping is affected). Documents locus-construction sensitivity.
      geom <- lapply(c(1, 5, 10, 25), function(D) {
        pos_D <- pos0
        after_c1 <- which(marker_ids %in% region_members) >
          match(causals[1], marker_ids)
        shift <- (D - 1) * 1000
        pos_D[after_c1] <- pos_D[after_c1] + shift
        tl_start <- min(pos_D[region_members])
        tl_end <- max(pos_D[region_members])
        split <- NA
        n_loci_truth <- 0L
        if (length(selected) > 0L) {
          lo <- define_associated_loci(
            scan, G, chromosome = rep("chr1", ncol(G)),
            position = unname(pos_D), selected_markers = selected,
            marker_ids = marker_ids, method = "physical",
            window_bp = window_bp
          )
          lids <- lo$loci
          lids_t <- lids[lids$start <= tl_end & lids$end >= tl_start, ,
                         drop = FALSE]
          n_loci_truth <- nrow(lids_t)
          c1_in <- lids$locus_id[lids$start <= pos_D[causals[1]] &
                                   lids$end >= pos_D[causals[1]]]
          c2_in <- lids$locus_id[lids$start <= pos_D[causals[2]] &
                                   lids$end >= pos_D[causals[2]]]
          split <- length(intersect(c1_in, c2_in)) == 0L &&
            (length(c1_in) + length(c2_in)) > 0L
          rm_ <- tryCatch(resolve_match(lo), error = function(e) NULL)
          S1_D <- if (is.null(rm_)) NA else rm_$S1
        } else {
          S1_D <- NA
        }
        data.frame(D = D, n_loci_overlap_truth = n_loci_truth,
                   causals_split = split, full_S1 = S1_D)
      })

      list(
        rep_id = rep_id, status = "ok",
        r2_causal = as.numeric(stats::cor(G[, causals[1]], G[, causals[2]])^2),
        op_pass = op_pass, op_pass_M1 = op_pass_M1,
        oc = oc, ol = ol, full = full,
        geom = do.call(rbind, geom)
      )
    }, error = function(e) {
      list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
    })
  }
  t0 <- Sys.time()
  p1_rows <- plapply(seq_len(100L), p1_rep)
  ok <- vapply(p1_rows, is_ok, logical(1))
  cat(sprintf("P1: ok %d/100, wall %.1fs\n", sum(ok),
              difftime(Sys.time(), t0, units = "secs")))
  if (any(!ok)) {
    print(table(vapply(p1_rows[!ok], function(x) as.character(x$status),
                       character(1))))
  }
  saveRDS(p1_rows, p1_rds)
} else {
  p1_rows <- readRDS(p1_rds)
  cat("P1: cached\n")
}
p1_rows <- p1_rows[vapply(p1_rows, is_ok, logical(1))]

p1_csv <- file.path(out_dir, "stage76_oracle_decomposition.csv")
if (!file.exists(p1_csv)) {
  n1 <- length(p1_rows)
  od <- data.frame(
    layer = c("oracle_primary", "oracle_primary_M1", "oracle_causal_set",
              "oracle_locus", "full"),
    recovery_count = c(
      sum(vapply(p1_rows, function(x) x$op_pass, logical(1))),
      sum(vapply(p1_rows, function(x) x$op_pass_M1, logical(1))),
      sum(vapply(p1_rows, function(x) x$oc$S1, logical(1))),
      sum(vapply(p1_rows, function(x) x$ol$S1, logical(1))),
      sum(vapply(p1_rows, function(x) x$full$S1, logical(1)))
    ),
    denom = n1
  )
  od$rate <- od$recovery_count / od$denom
  cand <- vapply(p1_rows, function(x) x$full$true_locus_n_candidates,
                 numeric(1))
  attr(od, "note") <- sprintf(
    "median full true-locus candidates: %s; oracle-locus candidates: 50",
    stats::median(cand))
  write.csv(od, p1_csv, row.names = FALSE)
  cat("median full true-locus candidates:", stats::median(cand),
      " (oracle-locus = 50)\n")
  print(od, row.names = FALSE)
}

p1g_csv <- file.path(out_dir, "stage76_geometry_sensitivity.csv")
if (!file.exists(p1g_csv)) {
  gg <- do.call(rbind, lapply(p1_rows, function(x) x$geom))
  gs <- do.call(rbind, lapply(split(gg, gg$D), function(d) {
    data.frame(
      causal_distance_kb = d$D[1],
      n = nrow(d),
      split_rate = mean(d$causals_split, na.rm = TRUE),
      full_S1_rate = mean(d$full_S1, na.rm = TRUE),
      mean_loci_overlap_truth = mean(d$n_loci_overlap_truth)
    )
  }))
  r2v <- vapply(p1_rows, function(x) x$r2_causal, numeric(1))
  gs$mean_causal_r2 <- mean(r2v)
  write.csv(gs, p1g_csv, row.names = FALSE)
  print(gs, row.names = FALSE)
}

## ==========================================================================
## P2: FDP null / single separation (multi reused from stage75)
## ==========================================================================
p2_rds <- file.path(out_dir, "p2_rows.rds")
if (!file.exists(p2_rds)) {
  p2_rep <- function(cfg, rep_id) {
    tryCatch({
      psid <- sprintf("stage76_fdp|%s", cfg$scenario)
      seed <- CondPED:::.seed_for_rep(20260817L, psid, rep_id)
      set.seed(seed)
      sim <- simulate_condped_data(
        n = 1000L, m = 4L, p = 1000L,
        experiment = "signal_resolution", scenario = cfg$scenario,
        locus_pve = 0.02, secondary_signal_pve = 0.0175,
        correlation = "block", tolerance = tolerance
      )
      if (!isTRUE(sim$status$ok)) {
        return(list(rep_id = rep_id, status = "sim_failed"))
      }
      G <- sim$G
      truth <- sim$truth
      marker_ids <- colnames(G)
      position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
      has_truth <- !is.null(truth$loci) && nrow(as.data.frame(truth$loci)) > 0
      if (has_truth) {
        tl <- as.data.frame(truth$loci)[1, ]
        in_region <- position >= tl$start & position <= tl$end
      } else {
        in_region <- rep(FALSE, length(marker_ids))
      }
      fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
      if (!isTRUE(fit$status$ok)) {
        return(list(rep_id = rep_id, status = "fit_failed"))
      }
      scan <- scan_mt_omnibus(fit, G)
      om <- scan$omnibus
      valid <- !is.na(om$p_value)
      padj <- rep(NA_real_, nrow(om))
      padj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
      selected <- om$marker_id[valid][padj[valid] <= alpha_omnibus]
      fp_markers <- sum(!in_region[match(selected, marker_ids)])
      n_novel <- 0L
      n_loci <- 0L
      n_signals <- 0L
      n_extra <- 0L
      if (length(selected) > 0L) {
        loci_obj <- define_associated_loci(
          scan, G, chromosome = rep("chr1", ncol(G)),
          position = unname(position), selected_markers = selected,
          marker_ids = marker_ids, method = "physical",
          window_bp = window_bp
        )
        n_loci <- nrow(loci_obj$loci)
        if (n_loci > 0L) {
          if (has_truth) {
            ov <- loci_obj$loci$start <= tl$end &
              loci_obj$loci$end >= tl$start
            n_novel <- sum(!ov)
          } else {
            n_novel <- n_loci
          }
          resolution <- resolve_locus_signals(
            fit, G, loci_obj, marker_ids = marker_ids,
            signal_adjust = "within_locus_bonferroni",
            alpha_signal = alpha_signal
          )
          n_signals <- nrow(resolution$signals)
          mm <- CondPED:::.match_signals(
            truth, list(signals = resolution$signals,
                        positions = position), G = G)
          n_extra <- length(mm$extra)
        }
      }
      list(rep_id = rep_id, status = "ok", scenario = cfg$scenario,
           n_selected = length(selected), fp_markers = fp_markers,
           n_loci = n_loci, n_novel_loci = n_novel,
           n_signals = n_signals, n_extra = n_extra)
    }, error = function(e) {
      list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
    })
  }
  t0 <- Sys.time()
  p2_rows <- list()
  for (sc in c("null", "single_multi_trait")) {
    rows <- plapply(seq_len(50L), function(i) p2_rep(list(scenario = sc), i))
    ok <- vapply(rows, is_ok, logical(1))
    cat(sprintf("P2 %s: ok %d/50\n", sc, sum(ok)))
    p2_rows <- c(p2_rows, rows[ok])
  }
  cat(sprintf("P2 wall %.1fs\n", difftime(Sys.time(), t0, units = "secs")))
  saveRDS(p2_rows, p2_rds)
} else {
  p2_rows <- readRDS(p2_rds)
  cat("P2: cached\n")
}

p2_csv <- file.path(out_dir, "stage76_fdr_propagation.csv")
if (!file.exists(p2_csv)) {
  mk <- function(rows, label) {
    n <- length(rows)
    s <- function(f) vapply(rows, f, numeric(1))
    data.frame(
      scenario = label, n = n,
      marker_fp_total = sum(s(function(x) x$fp_markers)),
      selected_total = sum(s(function(x) x$n_selected)),
      marker_fdp = sum(s(function(x) x$fp_markers)) /
        max(sum(s(function(x) x$n_selected)), 1),
      reps_with_fp = mean(s(function(x) x$n_selected) > 0),
      novel_loci_total = sum(s(function(x) x$n_novel_loci)),
      loci_total = sum(s(function(x) x$n_loci)),
      locus_fdp = sum(s(function(x) x$n_novel_loci)) /
        max(sum(s(function(x) x$n_loci)), 1),
      extra_signals_total = sum(s(function(x) x$n_extra)),
      signals_total = sum(s(function(x) x$n_signals)),
      signal_fdp = sum(s(function(x) x$n_extra)) /
        max(sum(s(function(x) x$n_signals)), 1)
    )
  }
  fdp_null <- mk(p2_rows[vapply(p2_rows, function(x) x$scenario == "null",
                                logical(1))], "null")
  fdp_single <- mk(p2_rows[vapply(
    p2_rows, function(x) x$scenario == "single_multi_trait", logical(1))],
    "single_multi_trait")
  ## multi reused from Stage 7.5 mixed_multisignal propagation
  mx <- read.csv("inst/validation/output/stage75/stage75_fdp_propagation.csv")
  fdp_multi <- data.frame(
    scenario = "mixed_multisignal (stage75)", n = 100,
    marker_fp_total = mx$count[mx$stage == "null_marker_fp"],
    selected_total = mx$count[mx$stage == "selected_markers"],
    marker_fdp = mx$count[mx$stage == "null_marker_fp"] /
      mx$count[mx$stage == "selected_markers"],
    reps_with_fp = NA_real_,
    novel_loci_total = mx$count[mx$stage == "novel_loci"],
    loci_total = mx$count[mx$stage == "estimated_loci"],
    locus_fdp = mx$count[mx$stage == "novel_loci"] /
      mx$count[mx$stage == "estimated_loci"],
    extra_signals_total = mx$count[mx$stage == "extra_signals"],
    signals_total = mx$count[mx$stage == "estimated_signals"],
    signal_fdp = mx$count[mx$stage == "extra_signals"] /
      mx$count[mx$stage == "estimated_signals"]
  )
  fdp <- rbind(fdp_null, fdp_single, fdp_multi)
  write.csv(fdp, p2_csv, row.names = FALSE)
  print(fdp, row.names = FALSE)
}

## ==========================================================================
## P3: LD-mixing gradient (linked_pseudo_multitrait, target_r2 grid)
## ==========================================================================
p3_rds <- file.path(out_dir, "p3_rows.rds")
if (!file.exists(p3_rds)) {
  p3_rep <- function(tr2, rep_id) {
    tryCatch({
      psid <- sprintf("stage76_ldmix|r2%.1f", tr2)
      seed <- CondPED:::.seed_for_rep(20260817L, psid, rep_id)
      set.seed(seed)
      sim <- simulate_condped_data(
        n = 1000L, m = 4L, p = 1000L,
        experiment = "end_to_end", scenario = "linked_pseudo_multitrait",
        locus_pve = 0.02, secondary_signal_pve = 0.0175,
        target_r2 = tr2, correlation = "block", tolerance = tolerance
      )
      if (!isTRUE(sim$status$ok)) {
        return(list(rep_id = rep_id, status = "sim_failed"))
      }
      G <- sim$G
      truth <- sim$truth
      causals <- truth$causal_markers
      fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
      if (!isTRUE(fit$status$ok)) {
        return(list(rep_id = rep_id, status = "fit_failed"))
      }
      trait_names <- fit$trait_names
      eff1 <- estimate_mt_effects(fit, G, targets = causals[1])
      beta1 <- stats::setNames(eff1$effects_long$beta,
                               eff1$effects_long$trait)[trait_names]
      X_C <- G[, causals[2], drop = FALSE]
      proj <- CondPED:::.build_conditional_projection(fit, X_C)
      cond1 <- CondPED:::.conditional_mt_scan(
        proj, G[, causals[1], drop = FALSE],
        marker_ids = causals[1], return_effects = TRUE
      )
      beta1c <- if (!is.null(cond1$effects)) {
        stats::setNames(cond1$effects$beta[1, ], trait_names)
      } else {
        stats::setNames(rep(NA_real_, length(trait_names)), trait_names)
      }
      ## real ASSET on the marginal lead causal
      lead_se <- stats::setNames(eff1$effects_long$se,
                                 eff1$effects_long$trait)[trait_names]
      SZ <- CondPED:::.sigma_z_from_cov(eff1$covariance[, , 1L])
      asset <- if (is.null(SZ)) {
        list(asset_best_subset = character(), status = "rank_deficient")
      } else {
        CondPED:::.run_asset_comparison(
          beta1, se = lead_se, Sigma_Z = SZ, trait_names = trait_names,
          sample_size = nrow(fit$rotation$Y_tilde), backend = NULL
        )
      }
      list(
        rep_id = rep_id, status = "ok", target_r2 = tr2,
        r2_empirical = as.numeric(stats::cor(G[, causals[1]],
                                             G[, causals[2]])^2),
        contam_marginal = abs(beta1["Trait2"]),
        contam_conditional = abs(beta1c["Trait2"]),
        asset_subset = paste(asset$asset_best_subset, collapse = "+"),
        asset_mixed = all(c("Trait1", "Trait2") %in% asset$asset_best_subset)
      )
    }, error = function(e) {
      list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
    })
  }
  t0 <- Sys.time()
  p3_rows <- list()
  for (tr2 in c(0, 0.3, 0.6, 0.9)) {
    rows <- plapply(seq_len(50L), function(i) p3_rep(tr2, i))
    ok <- vapply(rows, is_ok, logical(1))
    cat(sprintf("P3 r2=%.1f: ok %d/50\n", tr2, sum(ok)))
    p3_rows <- c(p3_rows, rows[ok])
  }
  cat(sprintf("P3 wall %.1fs\n", difftime(Sys.time(), t0, units = "secs")))
  saveRDS(p3_rows, p3_rds)
} else {
  p3_rows <- readRDS(p3_rds)
  cat("P3: cached\n")
}

p3_csv <- file.path(out_dir, "stage76_ld_mixing.csv")
if (!file.exists(p3_csv)) {
  ld <- do.call(rbind, lapply(split(
    p3_rows, vapply(p3_rows, function(x) x$target_r2, numeric(1))),
    function(rows) {
      data.frame(
        target_r2 = rows[[1]]$target_r2,
        n = length(rows),
        mean_empirical_r2 = mean(vapply(rows, function(x) x$r2_empirical,
                                        numeric(1))),
        mean_marginal_contam = mean(vapply(rows, function(x) x$contam_marginal,
                                           numeric(1))),
        mean_cond_contam = mean(vapply(rows, function(x) x$contam_conditional,
                                       numeric(1))),
        contam_reduction = mean(vapply(rows, function(x) x$contam_marginal,
                                       numeric(1))) -
          mean(vapply(rows, function(x) x$contam_conditional, numeric(1))),
        P_asset_mixed = mean(vapply(rows, function(x) x$asset_mixed,
                                    logical(1)))
      )
    }))
  rownames(ld) <- NULL
  write.csv(ld, p3_csv, row.names = FALSE)
  print(ld, row.names = FALSE)
}

## ==========================================================================
## P4: attribution rho/tau boundary pilot
## ==========================================================================
p4_rds <- file.path(out_dir, "p4_rows.rds")
if (!file.exists(p4_rds)) {
  p4_rep <- function(rho, rep_id) {
    tryCatch({
      psid <- sprintf("stage76_rho|%.3f", rho)
      seed <- CondPED:::.seed_for_rep(20260817L, psid, rep_id)
      set.seed(seed)
      sim <- simulate_condped_data(
        n = 1000L, m = 4L, p = 1000L,
        experiment = "end_to_end", scenario = "single_highly_representable",
        locus_pve = 0.02, target_loss = rho,
        correlation = "block", tolerance = tolerance
      )
      if (!isTRUE(sim$status$ok)) {
        return(list(rep_id = rep_id, status = "sim_failed"))
      }
      G <- sim$G
      truth <- sim$truth
      marker_ids <- colnames(G)
      position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
      fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
      if (!isTRUE(fit$status$ok)) {
        return(list(rep_id = rep_id, status = "fit_failed"))
      }
      scan <- scan_mt_omnibus(fit, G)
      om <- scan$omnibus
      valid <- !is.na(om$p_value)
      padj <- rep(NA_real_, nrow(om))
      padj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
      selected <- om$marker_id[valid][padj[valid] <= alpha_omnibus]
      rep_exact <- NA
      cand_exact <- NA
      if (length(selected) > 0L) {
        loci_obj <- define_associated_loci(
          scan, G, chromosome = rep("chr1", ncol(G)),
          position = unname(position), selected_markers = selected,
          marker_ids = marker_ids, method = "physical",
          window_bp = window_bp
        )
        if (nrow(loci_obj$loci) > 0L) {
          resolution <- resolve_locus_signals(
            fit, G, loci_obj, marker_ids = marker_ids,
            signal_adjust = "within_locus_bonferroni",
            alpha_signal = alpha_signal)
          pipe <- CondPED:::.pipeline_condped_full(
            fit, resolution, alpha_trait = alpha_trait,
            tolerance = tolerance)
          est <- resolution$signals
          if (nrow(est) > 0L) {
            est$position <- unname(position[est$representative_snp])
          }
          estimates <- list(
            signals = est,
            candidate_sets = pipe$candidate_sets,
            minimum_representative_sets =
              pipe$subset_analysis$minimum_representative_sets,
            irreducible_modules =
              pipe$subset_analysis$irreducible_modules,
            subset_table = pipe$subset_analysis$subset_table,
            positions = position
          )
          ev <- tryCatch(
            CondPED:::.evaluate_one_replicate(
              truth, estimates,
              list(experiment = "end_to_end",
                   architecture = "single_highly_representable",
                   n = 1000, locus_pve = 0.02, correlation = "block",
                   target_loss = rho, tolerance = tolerance,
                   alpha_omnibus = alpha_omnibus),
              runtime = NA_real_, G = G),
            error = function(e) NULL)
          if (!is.null(ev)) {
            rep_exact <- ev$rep_exact_family_recovery
            cand_exact <- ev$candidate_exact_recovery
          }
        }
      }
      list(rep_id = rep_id, status = "ok", rho = rho,
           rep_exact = rep_exact, cand_exact = cand_exact)
    }, error = function(e) {
      list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
    })
  }
  t0 <- Sys.time()
  p4_rows <- list()
  ## note: rho >= tau is unreachable by design for accept-R1 scenarios
  ## (single_highly_representable requires Rep = {Trait1} within
  ## tolerance; target_loss 0.09/0.25 fails PSD acceptance ->
  ## "invalid_covariance"). The boundary pilot therefore scans rho below
  ## tau: clearly below (0.02) / mid (0.05) / near boundary (0.08).
  for (rho in c(0.02, 0.05, 0.08)) {
    rows <- plapply(seq_len(50L), function(i) p4_rep(rho, i))
    ok <- vapply(rows, is_ok, logical(1))
    cat(sprintf("P4 rho=%.2f: ok %d/50\n", rho, sum(ok)))
    p4_rows <- c(p4_rows, rows[ok])
  }
  cat(sprintf("P4 wall %.1fs\n", difftime(Sys.time(), t0, units = "secs")))
  saveRDS(p4_rows, p4_rds)
} else {
  p4_rows <- readRDS(p4_rds)
  cat("P4: cached\n")
}

p4_csv <- file.path(out_dir, "stage76_attribution_boundary.csv")
if (!file.exists(p4_csv)) {
  ab <- do.call(rbind, lapply(split(
    p4_rows, vapply(p4_rows, function(x) x$rho, numeric(1))),
    function(rows) {
      rv <- vapply(rows, function(x) x$rep_exact, numeric(1))
      cv <- vapply(rows, function(x) x$cand_exact, numeric(1))
      data.frame(
        target_rho = rows[[1]]$rho, n = length(rows),
        n_eval = sum(!is.na(rv)),
        rep_exact_rate = mean(rv[!is.na(rv)] >= 1 - 1e-8),
        cand_exact_rate = mean(cv[!is.na(cv)] >= 1 - 1e-8)
      )
    }))
  rownames(ab) <- NULL
  write.csv(ab, p4_csv, row.names = FALSE)
  print(ab, row.names = FALSE)
}

## ==========================================================================
## P5: PVE ratio descriptive + final recommendation
## ==========================================================================
p5_csv <- file.path(out_dir, "stage76_pve_sensitivity.csv")
if (!file.exists(p5_csv)) {
  cal <- read.csv(
    "inst/validation/output/stage74/oracle_calibration/calibration_summary.csv")
  cal$pve_ratio_secondary_over_primary <- cal$secondary_signal_pve / 0.02
  write.csv(cal, p5_csv, row.names = FALSE)
}

p5r_csv <- file.path(out_dir, "stage76_final_recommendation.csv")
if (!file.exists(p5r_csv)) {
  rec <- data.frame(
    conclusion = c("A_simulation_design", "B_evaluator_metric",
                   "C_method_design"),
    finding = c(
      paste("causal 49kb geometry (fixed 7.4A), secondary PVE under-power",
            "(fixed 7.4C to 0.0175, ratio 0.875 vs primary 0.02),",
            "oracle-locus benchmark with 50-marker candidate set creates",
            "M=49 Bonferroni burden mismatch vs realistic small full loci"),
      paste("signal_count_exact conflates true-signal recovery with",
            "genome-wide FDR control; locus/signal FDP inherit marker BH",
            "FPs 1:1 (no aggregation amplification); secondary recovery",
            "metric fixed in 7.4B"),
      paste("see stage76_oracle_decomposition: gap oracle-primary ->",
            "oracle-causal-set -> oracle-locus quantifies the resolver",
            "multiplicity/screening loss with frozen thresholds")
    ),
    stringsAsFactors = FALSE
  )
  write.csv(rec, p5r_csv, row.names = FALSE)
}

cat("\nSTAGE76 ANALYSIS DONE\n")
