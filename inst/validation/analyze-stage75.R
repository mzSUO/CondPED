## Stage 7.5: full-to-oracle failure decomposition audit.
## Pure diagnostic: NO CondPED statistic, threshold, locus definition,
## matching, or simulation parameter is modified. All replicates reuse the
## deterministic Stage 7.4 seeds, so the simulated datasets are identical
## to Stage 7.4 and the three analysis layers are strictly paired by
## (scenario, spve, rep_id).
##
## Part 1: two_linked_trait_specific @ spve 0.0175 and 0.0075
##   (master_seed 20260815, same as Stage 7.4 D) - layer states,
##   transitions, failure classes A-E, oracle-primary sanity.
## Part 2: mixed_multisignal @ spve 0.0175 (master_seed 20260816, same as
##   Stage 7.4 G) - extra-signal provenance + FDP propagation.
## Part 3: ldmix (reads Stage 7.4 E/F RDS) - LD-mixing evidence.
devtools::load_all(quiet = TRUE)

workers <- as.integer(Sys.getenv("STAGE75_WORKERS", "16"))
reps <- 100L
out_dir <- "inst/validation/output/stage75"
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

## ==========================================================================
## PART 1: two_linked paired three-layer audit
## ==========================================================================

two_linked_rep <- function(rep_id, spve) {
  tryCatch({
    scenario <- "two_linked_trait_specific"
    psid <- sprintf("%s|lp%.3f_sp%.4f", scenario, 0.02, spve)
    seed <- CondPED:::.seed_for_rep(20260815L, psid, rep_id)
    set.seed(seed)
    sim <- simulate_condped_data(
      n = 1000L, m = 4L, p = 1000L,
      experiment = "signal_resolution", scenario = scenario,
      locus_pve = 0.02, secondary_signal_pve = spve,
      correlation = "block", tolerance = tolerance
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
    pri_id <- truth_signals$signal_id[1]
    sec_id <- truth_signals$signal_id[2]
    truth_sets <- truth$candidate_traits

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
    in_region <- position >= tl$start & position <= tl$end
    region_members <- marker_ids[in_region]

    ## ---- Layer 0: discovery -------------------------------------------
    D1 <- any(selected %in% region_members)

    ## ---- Layer 1: locus construction -----------------------------------
    loci_obj <- NULL
    if (length(selected) > 0L) {
      loci_obj <- define_associated_loci(
        scan, G, chromosome = rep("chr1", ncol(G)),
        position = unname(position), selected_markers = selected,
        marker_ids = marker_ids, method = "physical", window_bp = window_bp
      )
    }
    causal_in_locus <- function(loci_df, causal) {
      if (is.null(loci_df) || nrow(loci_df) == 0L) return(character())
      pos <- position[causal]
      loci_df$locus_id[loci_df$start <= pos & loci_df$end >= pos]
    }
    loci_df <- if (!is.null(loci_obj)) loci_obj$loci else NULL
    lids_pri <- causal_in_locus(loci_df, causals[1])
    lids_sec <- causal_in_locus(loci_df, causals[2])
    shared_lid <- intersect(lids_pri, lids_sec)
    L_state <- if (!D1) {
      "no_discovery"
    } else if (length(shared_lid) > 0L) {
      "joint_locus"                    # both causals in one est locus
    } else if (length(lids_pri) > 0L && length(lids_sec) > 0L) {
      "split"                          # causals in different est loci
    } else {
      "partial"                        # only one causal covered
    }

    ## ---- resolution + matching helper ----------------------------------
    resolve_and_match <- function(lo) {
      res <- resolve_locus_signals(fit, G, lo, marker_ids = marker_ids,
                                   signal_adjust = "within_locus_bonferroni",
                                   alpha_signal = alpha_signal)
      est <- res$signals
      mm <- CondPED:::.match_signals(
        truth, list(signals = est, positions = position), G = G)
      matched_truth <- unique(mm$matches$truth_signal_id)
      list(
        res = res, est = est, mm = mm,
        n_signals = nrow(est),
        n_matched = nrow(mm$matches),
        n_extra = length(mm$extra),
        P1 = pri_id %in% matched_truth,
        S1 = nrow(mm$matches) >= n_truth &&
          all(truth_signals$signal_id %in% matched_truth),
        R1 = nrow(est) == n_truth && nrow(mm$matches) == n_truth &&
          length(mm$extra) == 0L,
        count_exact = nrow(est) == n_truth
      )
    }

    ## ---- attribution helper: candidate exactness on matched pairs -----
    attribution_audit <- function(resolution, mm) {
      pipe <- CondPED:::.pipeline_condped_full(fit, resolution,
                                               alpha_trait = alpha_trait,
                                               tolerance = tolerance)
      cand <- pipe$candidate_sets
      if (nrow(mm$matches) == 0L) {
        return(list(A1 = NA, n_cand_exact = 0L, cand = cand,
                    min_rep = pipe$subset_analysis$minimum_representative_sets,
                    irr = pipe$subset_analysis$irreducible_modules))
      }
      exact <- vapply(seq_len(nrow(mm$matches)), function(i) {
        t <- mm$matches$truth_signal_id[i]
        e <- mm$matches$est_signal_id[i]
        isTRUE(setequal(cand[[e]], truth_sets[[t]]))
      }, logical(1))
      list(A1 = all(exact) && nrow(mm$matches) == n_truth,
           n_cand_exact = sum(exact), cand = cand,
           min_rep = pipe$subset_analysis$minimum_representative_sets,
           irr = pipe$subset_analysis$irreducible_modules)
    }

    ## ---- full layer -----------------------------------------------------
    full <- list(n_signals = 0L, n_matched = 0L, n_extra = 0L,
                 P1 = FALSE, S1 = FALSE, R1 = FALSE, count_exact = FALSE,
                 A1 = NA, n_cand_exact = 0L)
    full_res <- NULL
    if (!is.null(loci_obj) && nrow(loci_obj$loci) > 0L) {
      rm_full <- resolve_and_match(loci_obj)
      full_res <- rm_full$res
      att <- attribution_audit(rm_full$res, rm_full$mm)
      full <- modifyList(rm_full[c("n_signals", "n_matched", "n_extra",
                                   "P1", "S1", "R1", "count_exact")],
                         list(A1 = att$A1, n_cand_exact = att$n_cand_exact))
    }

    ## ---- oracle-locus layer ---------------------------------------------
    oracle_lead <- region_members[which.min(
      om$p_value[match(region_members, om$marker_id)])]
    oracle_locus <- list(
      loci = data.frame(
        locus_id = tl$locus_id, chromosome = tl$chromosome,
        start = tl$start, end = tl$end, lead_snp = oracle_lead,
        lead_p = om$p_value[match(oracle_lead, om$marker_id)],
        n_significant_markers = length(region_members),
        n_region_markers = length(region_members), status = "ok",
        stringsAsFactors = FALSE
      ),
      membership = data.frame(
        locus_id = tl$locus_id, marker_id = region_members,
        significant_in_marginal_scan = region_members %in% selected,
        r2_to_lead = as.numeric(stats::cor(G[, region_members],
                                           G[, oracle_lead])^2),
        position = unname(position[region_members]),
        stringsAsFactors = FALSE
      )
    )
    rm_oracle <- resolve_and_match(oracle_locus)
    att_o <- attribution_audit(rm_oracle$res, rm_oracle$mm)
    oracle <- c(rm_oracle[c("n_signals", "n_matched", "n_extra",
                            "P1", "S1", "R1", "count_exact")],
                list(A1 = att_o$A1, n_cand_exact = att_o$n_cand_exact))

    ## ---- oracle-primary layer -------------------------------------------
    X_C <- G[, causals[1], drop = FALSE]
    proj <- CondPED:::.build_conditional_projection(fit, X_C)
    cs <- CondPED:::.conditional_mt_scan(
      proj, G[, causals[2], drop = FALSE],
      marker_ids = causals[2], return_effects = FALSE
    )
    ctab <- cs$conditional
    M_remaining <- length(region_members) - 1L
    thr <- alpha_signal / M_remaining
    marg_p <- om$p_value[match(causals[2], om$marker_id)]
    op <- list(
      marginal_p = marg_p,
      cond_Q = ctab$Q, cond_p = ctab$p_value, cond_df = ctab$df,
      cond_rank = ctab$rank_J, cond_status = ctab$status,
      threshold = thr,
      pass = isTRUE(ctab$p_value <= thr),
      rank_deficient = identical(ctab$status, "rank_deficient"),
      numerical_failure = !ctab$status %in% c("ok", "rank_deficient")
    )

    ## ---- failure classification (full layer, priority A > B > C > D) ---
    cls <- "success"
    sub <- ""
    if (!full$R1 || !isTRUE(full$A1)) {
      if (!D1) {
        cls <- "A"; sub <- "A_discovery"
      } else if (L_state != "joint_locus") {
        cls <- "B"
        sub <- switch(L_state,
                      split = "B1_physical_window_split",
                      partial = "B2_association_region_partial",
                      "B4_other")
      } else if (!full$S1) {
        cls <- "C"
        sub <- if (!op$pass) {
          "C1_oracle_primary_power"
        } else if (!full$P1) {
          "C2_primary_not_matched"
        } else if (full$n_signals < n_truth) {
          "C3_secondary_screening"
        } else {
          "C5_signal_matching"
        }
      } else {
        cls <- "D"
        sub <- if (isTRUE(full$A1)) {
          "D6_other"
        } else if (full$n_cand_exact < n_truth) {
          "D1_candidate_attribution"
        } else {
          "D6_other"
        }
      }
    }

    list(
      rep_id = rep_id, status = "ok", spve = spve,
      D1 = D1, L_state = L_state,
      n_selected = length(selected),
      n_selected_in_region = sum(selected %in% region_members),
      full = full, oracle = oracle, oracle_primary = op,
      failure_class = cls, failure_subclass = sub,
      E_extra_novel = full$n_extra
    )
  }, error = function(e) {
    list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
  })
}

run_condition <- function(spve) {
  tag <- sprintf("sp%.4f", spve)
  rds <- file.path(out_dir, paste0("two_linked_rows_", tag, ".rds"))
  if (file.exists(rds)) {
    out <- readRDS(rds)
    if (length(out) == reps && all(vapply(out, is_ok, logical(1)))) {
      cat("SKIP (done):", tag, "\n")
      return(out)
    }
  }
  t0 <- Sys.time()
  if (.Platform$OS.type == "unix") {
    out <- parallel::mclapply(seq_len(reps), two_linked_rep, spve = spve,
                              mc.cores = workers)
  } else {
    out <- lapply(seq_len(reps), two_linked_rep, spve = spve)
  }
  ok <- vapply(out, is_ok, logical(1))
  cat(sprintf("two_linked %s: ok %d/%d, wall %.1fs\n", tag, sum(ok),
              length(out), difftime(Sys.time(), t0, units = "secs")))
  if (any(!ok)) {
    print(table(vapply(out[!ok], function(x) as.character(x$status),
                       character(1))))
  }
  saveRDS(out, rds)
  out
}

rows_0175 <- run_condition(0.0175)
rows_0075 <- run_condition(0.0075)

## ---- per-replicate transition table --------------------------------------
trans_rows <- function(rows, spve) {
  ok <- vapply(rows, is_ok, logical(1))
  rows <- rows[ok]
  data.frame(
    spve = spve,
    rep_id = vapply(rows, function(x) x$rep_id, integer(1)),
    D1 = vapply(rows, function(x) x$D1, logical(1)),
    L_state = vapply(rows, function(x) x$L_state, character(1)),
    full_S1 = vapply(rows, function(x) x$full$S1, logical(1)),
    full_R1 = vapply(rows, function(x) x$full$R1, logical(1)),
    full_A1 = vapply(rows, function(x) isTRUE(x$full$A1), logical(1)),
    full_count_exact = vapply(rows, function(x) x$full$count_exact, logical(1)),
    oracle_S1 = vapply(rows, function(x) x$oracle$S1, logical(1)),
    oracle_R1 = vapply(rows, function(x) x$oracle$R1, logical(1)),
    oracle_A1 = vapply(rows, function(x) isTRUE(x$oracle$A1), logical(1)),
    oracle_count_exact = vapply(rows, function(x) x$oracle$count_exact,
                                logical(1)),
    opass = vapply(rows, function(x) x$oracle_primary$pass, logical(1)),
    failure_class = vapply(rows, function(x) x$failure_class, character(1)),
    failure_subclass = vapply(rows, function(x) x$failure_subclass,
                              character(1)),
    stringsAsFactors = FALSE
  )
}
trans <- rbind(trans_rows(rows_0175, 0.0175), trans_rows(rows_0075, 0.0075))
write.csv(trans, file.path(out_dir, "stage75_per_replicate_transition.csv"),
          row.names = FALSE)

## ---- failure decomposition ------------------------------------------------
fd <- do.call(rbind, lapply(split(trans, trans$spve), function(d) {
  tb <- as.data.frame(table(class = d$failure_class,
                            subclass = d$failure_subclass))
  tb$spve <- d$spve[1]
  tb
}))
fd <- fd[, c("spve", "class", "subclass", "Freq")]
write.csv(fd, file.path(out_dir, "stage75_failure_decomposition.csv"),
          row.names = FALSE)

## ---- oracle-primary table --------------------------------------------------
op_rows <- function(rows, spve) {
  ok <- vapply(rows, is_ok, logical(1))
  rows <- rows[ok]
  data.frame(
    spve = spve,
    rep_id = vapply(rows, function(x) x$rep_id, integer(1)),
    marginal_p = vapply(rows, function(x) x$oracle_primary$marginal_p,
                        numeric(1)),
    cond_Q = vapply(rows, function(x) x$oracle_primary$cond_Q, numeric(1)),
    cond_p = vapply(rows, function(x) x$oracle_primary$cond_p, numeric(1)),
    cond_df = vapply(rows, function(x) x$oracle_primary$cond_df, numeric(1)),
    cond_rank = vapply(rows, function(x) x$oracle_primary$cond_rank,
                       numeric(1)),
    cond_status = vapply(rows, function(x) x$oracle_primary$cond_status,
                         character(1)),
    threshold = vapply(rows, function(x) x$oracle_primary$threshold,
                       numeric(1)),
    pass = vapply(rows, function(x) x$oracle_primary$pass, logical(1)),
    full_S1 = vapply(rows, function(x) x$full$S1, logical(1)),
    oracle_S1 = vapply(rows, function(x) x$oracle$S1, logical(1)),
    stringsAsFactors = FALSE
  )
}
op_tab <- rbind(op_rows(rows_0175, 0.0175), op_rows(rows_0075, 0.0075))
write.csv(op_tab, file.path(out_dir, "stage75_oracle_primary.csv"),
          row.names = FALSE)

## ==========================================================================
## PART 2: mixed_multisignal extra-signal provenance + FDP propagation
## ==========================================================================

mixed_rep <- function(rep_id) {
  tryCatch({
    scenario <- "mixed_multisignal"
    psid <- sprintf("%s|lp%.3f_sp%.4f", scenario, 0.02, 0.0175)
    seed <- CondPED:::.seed_for_rep(20260816L, psid, rep_id)
    set.seed(seed)
    sim <- simulate_condped_data(
      n = 1000L, m = 4L, p = 1000L,
      experiment = "end_to_end", scenario = scenario,
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
    causals <- truth$causal_markers

    fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
    if (!isTRUE(fit$status$ok)) {
      return(list(rep_id = rep_id, status = "fit_failed"))
    }
    scan <- scan_mt_omnibus(fit, G)
    om <- scan$omnibus
    valid <- !is.na(om$p_value)
    padj <- rep(NA_real_, nrow(om))
    padj[valid] <- stats::p.adjust(om$p_value[valid], method = "BH")
    om$padj <- padj
    selected <- om$marker_id[valid][padj[valid] <= alpha_omnibus]

    tl <- as.data.frame(truth$loci)[1, ]
    in_region <- position >= tl$start & position <= tl$end
    r2_causal <- function(mid) max(stats::cor(G[, mid], G[, causals])^2)

    marker_cat <- vapply(selected, function(mid) {
      if (in_region[match(mid, marker_ids)]) "truth_region"
      else if (r2_causal(mid) >= 0.5) "causal_proxy"
      else "null_marker"
    }, character(1))

    out <- list(rep_id = rep_id, status = "ok",
                n_selected = length(selected),
                n_null_marker = sum(marker_cat == "null_marker"),
                n_truth_marker = sum(marker_cat == "truth_region"),
                extra = NULL, n_novel_locus = 0L, n_loci = 0L,
                novel_locus_widths = numeric(), novel_locus_nsig = integer(),
                truth_locus_nsig = 0L, n_signals = 0L)
    if (length(selected) == 0L) return(out)

    loci_obj <- define_associated_loci(
      scan, G, chromosome = rep("chr1", ncol(G)),
      position = unname(position), selected_markers = selected,
      marker_ids = marker_ids, method = "physical", window_bp = window_bp
    )
    if (nrow(loci_obj$loci) == 0L) return(out)
    out$n_loci <- nrow(loci_obj$loci)

    overlaps_truth <- loci_obj$loci$start <= tl$end &
      loci_obj$loci$end >= tl$start
    out$n_novel_locus <- sum(!overlaps_truth)
    novel_lids <- loci_obj$loci$locus_id[!overlaps_truth]
    memb <- loci_obj$membership
    out$novel_locus_widths <- (loci_obj$loci$end - loci_obj$loci$start)[
      !overlaps_truth]
    out$novel_locus_nsig <- vapply(novel_lids, function(lid) {
      sum(memb$locus_id == lid & memb$significant_in_marginal_scan)
    }, integer(1))
    truth_lids <- loci_obj$loci$locus_id[overlaps_truth]
    out$truth_locus_nsig <- sum(memb$locus_id %in% truth_lids &
                                  memb$significant_in_marginal_scan)

    resolution <- resolve_locus_signals(fit, G, loci_obj,
                                        marker_ids = marker_ids,
                                        signal_adjust = "within_locus_bonferroni",
                                        alpha_signal = alpha_signal)
    est <- resolution$signals
    out$n_signals <- nrow(est)
    mm <- CondPED:::.match_signals(
      truth, list(signals = est, positions = position), G = G)
    if (length(mm$extra) == 0L) return(out)

    pipe <- CondPED:::.pipeline_condped_full(fit, resolution,
                                             alpha_trait = alpha_trait,
                                             tolerance = tolerance)
    cand <- pipe$candidate_sets
    min_rep <- pipe$subset_analysis$minimum_representative_sets
    irr <- pipe$subset_analysis$irreducible_modules
    rep_sids <- if (nrow(min_rep) > 0) unique(min_rep$signal_id) else character()
    irr_sids <- if (nrow(irr) > 0) unique(irr$signal_id) else character()

    matched_est_loci <- unique(mm$matches$est_locus_id)
    matched_truth_loci <- unique(mm$matches$truth_locus_id)
    extra <- lapply(mm$extra, function(eid) {
      es <- est[est$signal_id == eid, ]
      el <- loci_obj$loci[loci_obj$loci$locus_id == es$locus_id, ]
      r2m <- r2_causal(es$representative_snp)
      ov <- el$start <= tl$end && el$end >= tl$start
      catg <- if (es$locus_id %in% matched_est_loci) {
        "within_locus_fp"
      } else if (ov) {
        if (tl$locus_id %in% matched_truth_loci && r2m >= 0.5) {
          "split_duplicate_proxy"
        } else {
          "split_region_fp"
        }
      } else {
        "novel_locus_fp"
      }
      mi <- match(es$representative_snp, om$marker_id)
      data.frame(
        rep_id = rep_id,
        est_signal_id = eid,
        marker_id = es$representative_snp,
        locus_id = es$locus_id,
        category = catg,
        novel_locus = catg == "novel_locus_fp",
        locus_n_sig_markers = sum(memb$locus_id == es$locus_id &
                                    memb$significant_in_marginal_scan),
        locus_width = el$end - el$start,
        marker_raw_p = om$p_value[mi],
        marker_adjusted_p = om$padj[mi],
        candidate_status = length(cand[[eid]]) > 0L,
        n_candidate_traits = length(cand[[eid]]),
        in_min_rep = eid %in% rep_sids,
        in_irr_module = eid %in% irr_sids,
        stringsAsFactors = FALSE
      )
    })
    out$extra <- do.call(rbind, extra)
    out
  }, error = function(e) {
    list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
  })
}

rds_m <- file.path(out_dir, "mixed_rows.rds")
if (file.exists(rds_m)) {
  mixed_rows <- readRDS(rds_m)
  if (!(length(mixed_rows) == reps &&
        all(vapply(mixed_rows, is_ok, logical(1))))) {
    mixed_rows <- NULL
  }
}
if (!exists("mixed_rows") || is.null(mixed_rows)) {
  t0 <- Sys.time()
  if (.Platform$OS.type == "unix") {
    mixed_rows <- parallel::mclapply(seq_len(reps), mixed_rep,
                                     mc.cores = workers)
  } else {
    mixed_rows <- lapply(seq_len(reps), mixed_rep)
  }
  okm <- vapply(mixed_rows, is_ok, logical(1))
  cat(sprintf("mixed: ok %d/%d, wall %.1fs\n", sum(okm), length(mixed_rows),
              difftime(Sys.time(), t0, units = "secs")))
  if (any(!okm)) {
    print(table(vapply(mixed_rows[!okm], function(x) as.character(x$status),
                       character(1))))
  }
  saveRDS(mixed_rows, rds_m)
}
mixed_rows <- mixed_rows[vapply(mixed_rows, is_ok, logical(1))]

prov <- do.call(rbind, lapply(mixed_rows, function(x) x$extra))
if (!is.null(prov)) {
  write.csv(prov, file.path(out_dir, "stage75_extra_signal_provenance.csv"),
            row.names = FALSE)
}

## ---- FDP propagation ------------------------------------------------------
n_sel <- vapply(mixed_rows, function(x) x$n_selected, numeric(1))
n_null <- vapply(mixed_rows, function(x) x$n_null_marker, numeric(1))
n_novel <- vapply(mixed_rows, function(x) x$n_novel_locus, numeric(1))
n_loci <- vapply(mixed_rows, function(x) x$n_loci, numeric(1))
n_sig <- vapply(mixed_rows, function(x) x$n_signals, numeric(1))
n_extra <- if (is.null(prov)) 0 else nrow(prov)
prop_tab <- data.frame(
  stage = c("selected_markers", "null_marker_fp", "estimated_loci",
            "novel_loci", "estimated_signals", "extra_signals",
            "extra_with_candidate", "extra_in_min_rep", "extra_in_irr"),
  count = c(sum(n_sel), sum(n_null), sum(n_loci), sum(n_novel), sum(n_sig),
            n_extra,
            if (is.null(prov)) 0 else sum(prov$candidate_status),
            if (is.null(prov)) 0 else sum(prov$in_min_rep),
            if (is.null(prov)) 0 else sum(prov$in_irr))
)
prop_tab$per_replicate <- prop_tab$count / length(mixed_rows)
prop_tab$retention_vs_previous <- c(NA_real_,
  prop_tab$count[-1] / prop_tab$count[-nrow(prop_tab)])
write.csv(prop_tab, file.path(out_dir, "stage75_fdp_propagation.csv"),
          row.names = FALSE)

## ==========================================================================
## PART 3: ldmix (reuse Stage 7.4 E/F rows)
## ==========================================================================
ld <- readRDS("inst/validation/output/stage74/ldmix/ldmix_rows.rds")
ld <- ld[vapply(ld, is_ok, logical(1))]
ldmix_tab <- data.frame(
  metric = c("n", "P_valid_mixing", "P_lead_has_T1", "P_lead_has_T2",
             "mean_marginal_contam", "mean_cond_contam", "mean_joint_contam",
             "lead_breadth_error", "resolved_breadth_error",
             "P_signal_count_recovered"),
  value = c(length(ld),
            mean(vapply(ld, function(x) x$valid_mixing, logical(1))),
            mean(vapply(ld, function(x) x$lead_has_T1, logical(1))),
            mean(vapply(ld, function(x) x$lead_has_T2, logical(1))),
            mean(vapply(ld, function(x) x$contam_marginal_snp1, numeric(1))),
            mean(vapply(ld, function(x) x$contam_cond_snp1, numeric(1))),
            mean(vapply(ld, function(x) x$contam_joint_snp1, numeric(1))),
            mean(vapply(ld, function(x) x$lead_breadth_err, numeric(1))),
            mean(unlist(lapply(ld, function(x) x$res_breadth_errs))),
            mean(vapply(ld, function(x) x$signal_count_recovered, logical(1))))
)

## ==========================================================================
## SUMMARY
## ==========================================================================
smry <- list()
add <- function(question, metric, count, denom) {
  smry[[length(smry) + 1L]] <<- data.frame(
    question = question, metric = metric, count = count, denom = denom,
    proportion = count / denom, stringsAsFactors = FALSE
  )
}
for (sp in c(0.0175, 0.0075)) {
  d <- trans[trans$spve == sp, ]
  n <- nrow(d)
  tag <- paste0("two_linked sp", sp)
  add(tag, "full count_exact", sum(d$full_count_exact), n)
  add(tag, "oracle-locus count_exact", sum(d$oracle_count_exact), n)
  add(tag, "full S1 (secondary, corrected)", sum(d$full_S1), n)
  add(tag, "oracle S1", sum(d$oracle_S1), n)
  add(tag, "oracle-primary pass", sum(d$opass), n)
  add(tag, "P(full S1 | oracle S1) count", sum(d$full_S1[d$oracle_S1]),
      sum(d$oracle_S1))
  add(tag, "P(oracle S1 | full S1) count", sum(d$oracle_S1[d$full_S1]),
      sum(d$full_S1))
  add(tag, "P(full S1 | oracle-primary pass)",
      sum(d$full_S1[d$opass]), sum(d$opass))
  add(tag, "P(oracle S1 | oracle-primary pass)",
      sum(d$oracle_S1[d$opass]), sum(d$opass))
  add(tag, "trans full fail -> oracle success",
      sum(!d$full_S1 & d$oracle_S1), n)
  add(tag, "trans full success -> oracle fail",
      sum(d$full_S1 & !d$oracle_S1), n)
  add(tag, "trans oracle fail -> oracle-primary success",
      sum(!d$oracle_S1 & d$opass), n)
  add(tag, "all success (full+oracle+op)", sum(d$full_S1 & d$oracle_S1 & d$opass), n)
  add(tag, "all fail", sum(!d$full_S1 & !d$oracle_S1 & !d$opass), n)
  add(tag, "P(final exact R1 | attribution exact A1) full",
      sum(d$full_R1[d$full_A1]), sum(d$full_A1))
  add(tag, "discovery D1", sum(d$D1), n)
  add(tag, "joint_locus L0", sum(d$L_state == "joint_locus"), n)
  op <- op_tab[op_tab$spve == sp, ]
  add(tag, "oracle-primary rank deficient", sum(op$cond_status == "rank_deficient"), nrow(op))
  add(tag, "oracle-primary non-ok status", sum(op$cond_status != "ok"), nrow(op))
}
smry <- do.call(rbind, smry)
rownames(smry) <- NULL

ldmix_print <- ldmix_tab
print(smry, row.names = FALSE)
write.csv(smry, file.path(out_dir, "stage75_summary.csv"), row.names = FALSE)
write.csv(ldmix_tab, file.path(out_dir, "stage75_ldmix_evidence.csv"),
          row.names = FALSE)

cat("\n== failure decomposition ==\n")
print(fd, row.names = FALSE)
cat("\n== FDP propagation ==\n")
print(prop_tab, row.names = FALSE)
cat("\n== extra-signal provenance ==\n")
if (!is.null(prov)) {
  print(table(prov$category))
  cat("P(extra enters candidate):", round(mean(prov$candidate_status), 3), "\n")
  cat("P(extra in min Rep set):", round(mean(prov$in_min_rep), 3), "\n")
  cat("P(extra in Irr module):", round(mean(prov$in_irr), 3), "\n")
  cat("median markers per novel-locus FP:",
      stats::median(prov$locus_n_sig_markers[prov$novel_locus]), "\n")
  cat("single-marker novel-locus FPs:",
      sum(prov$novel_locus & prov$locus_n_sig_markers == 1L), "/",
      sum(prov$novel_locus), "\n")
}
cat("\n== ldmix ==\n")
print(ldmix_print, row.names = FALSE)
cat("\nSTAGE75 ANALYSIS DONE\n")
