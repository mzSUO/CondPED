## Stage 7.7: simulation regime / evaluation design audit.
## Calibration pilots only (50 reps/cell, 4 workers). No production
## statistical method, threshold, or frozen definition is modified.
## Every cell uses deterministic seeds; all four analysis layers
## (full / oracle-locus / oracle-causal-set / oracle-primary) are computed
## on the SAME replicate and are therefore strictly paired.
##
##   P1 secondary PVE grid      (two_linked, 8 cells)
##   P2 sample-size grid        (two_linked, 6 cells)
##   P3 LD grid + contamination (two_linked, 5 cells)
##   P4 rho/tau boundary        (single_highly_representable, 4 cells)
##   P5 q=2 vs q=3 scalability  (two_linked vs three_linked)
devtools::load_all(quiet = TRUE)

workers <- as.integer(Sys.getenv("STAGE77_WORKERS", "4"))
reps <- as.integer(Sys.getenv("STAGE77_REPS", "50"))
master_seed <- 20260818L
out_dir <- "inst/validation/output/stage77"
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

## ---- unified per-replicate engine (two_linked-family) --------------------
engine_rep <- function(rep_id, cell_id, scenario, n, spve, tr2) {
  tryCatch({
    t0 <- proc.time()[["elapsed"]]
    seed <- CondPED:::.seed_for_rep(master_seed, cell_id, rep_id)
    set.seed(seed)
    sim <- simulate_condped_data(
      n = n, m = 4L, p = 1000L,
      experiment = "signal_resolution", scenario = scenario,
      locus_pve = 0.02, secondary_signal_pve = spve,
      target_r2 = tr2, correlation = "block", tolerance = tolerance
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
    truth_sets <- truth$candidate_traits
    truth_beta <- truth$beta

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

    mk_locus <- function(members) {
      pos_v <- unname(position[members])
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

    ## resolution + matching + per-pair attribution/direction
    resolve_match <- function(lo, with_attr = FALSE) {
      res <- resolve_locus_signals(fit, G, lo, marker_ids = marker_ids,
                                   signal_adjust = "within_locus_bonferroni",
                                   alpha_signal = alpha_signal)
      est <- res$signals
      mm <- CondPED:::.match_signals(
        truth, list(signals = est, positions = position), G = G)
      mt <- unique(mm$matches$truth_signal_id)
      out <- list(
        n_signals = nrow(est), n_matched = nrow(mm$matches),
        n_extra = length(mm$extra),
        recovery = nrow(mm$matches) >= n_truth &&
          all(truth_signals$signal_id %in% mt),
        count_exact = nrow(est) == n_truth,
        cand_exact = NA, dir_exact = NA, mm = mm
      )
      if (with_attr && nrow(mm$matches) > 0L) {
        pipe <- CondPED:::.pipeline_condped_full(
          fit, res, alpha_trait = alpha_trait, tolerance = tolerance)
        cand <- pipe$candidate_sets
        eff <- res$beta
        ce <- logical(nrow(mm$matches))
        de <- logical(nrow(mm$matches))
        for (i in seq_len(nrow(mm$matches))) {
          ti <- mm$matches$truth_signal_id[i]
          ei <- mm$matches$est_signal_id[i]
          ce[i] <- isTRUE(setequal(cand[[ei]], truth_sets[[ti]]))
          tb_all <- truth_beta$beta[truth_beta$signal_id == ti]
          tt_all <- truth_beta$trait[truth_beta$signal_id == ti]
          nz <- tb_all != 0                 # truth candidate (nonzero) traits
          tb <- tb_all[nz]
          tt <- tt_all[nz]
          rows_e <- eff[eff$signal_id == ei, ]
          eb <- rows_e$beta[match(tt, rows_e$trait)]
          de[i] <- length(eb) > 0L && !anyNA(eb) &&
            all(sign(eb) == sign(tb) | eb == 0)
        }
        out$cand_exact <- mean(ce)
        out$dir_exact <- mean(de)
      }
      out
    }

    ## full
    full <- list(n_signals = 0L, n_matched = 0L, n_extra = 0L,
                 recovery = FALSE, count_exact = FALSE,
                 cand_exact = NA, dir_exact = NA,
                 n_novel_extra = 0L, true_locus_candidates = 0L)
    if (length(selected) > 0L) {
      loci_obj <- define_associated_loci(
        scan, G, chromosome = rep("chr1", ncol(G)),
        position = unname(position), selected_markers = selected,
        marker_ids = marker_ids, method = "physical",
        window_bp = window_bp
      )
      if (nrow(loci_obj$loci) > 0L) {
        full <- resolve_match(loci_obj, with_attr = TRUE)
        full$n_novel_extra <- 0L
        memb <- loci_obj$membership
        tl_lids <- loci_obj$loci$locus_id[
          loci_obj$loci$start <= tl$end & loci_obj$loci$end >= tl$start]
        full$true_locus_candidates <- sum(memb$locus_id %in% tl_lids)
        ## novel-locus extras: extra signals whose est locus does not
        ## overlap the truth locus
        mm <- full$mm
        if (length(mm$extra) > 0L) {
          est_sig <- mm$est_signals
          extra_lids <- est_sig$locus_id[est_sig$signal_id %in% mm$extra]
          nov <- loci_obj$loci[match(extra_lids, loci_obj$loci$locus_id), ]
          full$n_novel_extra <- sum(!(nov$start <= tl$end &
                                        nov$end >= tl$start))
        }
        full$mm <- NULL
      }
    }

    ## oracle-locus / oracle-causal-set
    ol <- resolve_match(mk_locus(region_members))
    oc <- resolve_match(mk_locus(causals))
    ol$mm <- NULL
    oc$mm <- NULL

    ## oracle-primary: condition on true primary, test each remaining causal
    X_C <- G[, causals[1], drop = FALSE]
    proj <- CondPED:::.build_conditional_projection(fit, X_C)
    cs <- CondPED:::.conditional_mt_scan(
      proj, G[, causals[-1], drop = FALSE],
      marker_ids = causals[-1], return_effects = FALSE
    )
    ctab <- cs$conditional
    M_region <- length(region_members) - 1L
    op <- list(
      median_Q = stats::median(ctab$Q),
      pass_Mregion = all(ctab$p_value <= alpha_signal / M_region),
      pass_M1 = all(ctab$p_value <= alpha_signal),
      rank_deficient = any(ctab$status == "rank_deficient"),
      numerical_failure = any(!ctab$status %in% c("ok", "rank_deficient"))
    )

    ## cross-trait contamination (Trait1-causal's effect on Trait2 etc.)
    eff1 <- estimate_mt_effects(fit, G, targets = causals[1])
    tn <- fit$trait_names
    b1 <- stats::setNames(eff1$effects_long$beta,
                          eff1$effects_long$trait)[tn]
    contam_marg <- abs(unname(b1["Trait2"]))
    X_C2 <- G[, causals[2], drop = FALSE]
    proj2 <- CondPED:::.build_conditional_projection(fit, X_C2)
    cond1 <- CondPED:::.conditional_mt_scan(
      proj2, G[, causals[1], drop = FALSE],
      marker_ids = causals[1], return_effects = TRUE
    )
    contam_cond <- if (!is.null(cond1$effects)) {
      abs(unname(cond1$effects$beta[1, "Trait2"]))
    } else {
      NA_real_
    }

    list(
      rep_id = rep_id, status = "ok",
      runtime = proc.time()[["elapsed"]] - t0,
      r2_empirical = as.numeric(stats::cor(G[, causals[1]],
                                           G[, causals[2]])^2),
      n_selected = length(selected),
      marker_fp = sum(selected %in% marker_ids[!in_region]),
      full = full, oracle_locus = ol, oracle_causal_set = oc,
      oracle_primary = op,
      contam_marginal = contam_marg, contam_conditional = contam_cond
    )
  }, error = function(e) {
    list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
  })
}

summarize_cell <- function(rows, cell, extras = list()) {
  rows <- rows[vapply(rows, is_ok, logical(1))]
  n <- length(rows)
  g <- function(f) vapply(rows, function(x) {
    v <- f(x)
    if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v)
  }, numeric(1))
  h <- function(f) vapply(rows, function(x) {
    v <- f(x)
    if (is.null(v) || length(v) == 0L) NA else as.logical(v)
  }, logical(1))
  data.frame(
    extras,
    cell = cell, n_ok = n,
    full_recovery = mean(h(function(x) x$full$recovery)),
    full_count_exact = mean(h(function(x) x$full$count_exact)),
    oracle_locus_recovery = mean(h(function(x) x$oracle_locus$recovery)),
    oracle_causal_set_recovery = mean(h(function(x) x$oracle_causal_set$recovery)),
    oracle_primary_Mregion = mean(h(function(x) x$oracle_primary$pass_Mregion)),
    oracle_primary_M1 = mean(h(function(x) x$oracle_primary$pass_M1)),
    attribution_cand_exact = mean(g(function(x) x$full$cand_exact),
                                  na.rm = TRUE),
    direction_exact = mean(g(function(x) x$full$dir_exact), na.rm = TRUE),
    extra_signal_rate = mean(g(function(x) x$full$n_extra) > 0),
    novel_extra_rate = mean(g(function(x) x$full$n_novel_extra) > 0),
    marker_fdp = sum(g(function(x) x$marker_fp)) /
      max(sum(g(function(x) x$n_selected)), 1),
    median_cond_Q = stats::median(g(function(x) x$oracle_primary$median_Q)),
    rank_deficient = sum(h(function(x) x$oracle_primary$rank_deficient)),
    numerical_failure = sum(h(function(x) x$oracle_primary$numerical_failure)),
    contam_marginal = mean(g(function(x) x$contam_marginal)),
    contam_conditional = mean(g(function(x) x$contam_conditional),
                              na.rm = TRUE),
    median_candidates_full = stats::median(
      g(function(x) x$full$true_locus_candidates)),
    mean_empirical_r2 = mean(g(function(x) x$r2_empirical)),
    mean_runtime = mean(g(function(x) x$runtime)),
    stringsAsFactors = FALSE
  )
}

run_cell <- function(cell_id, scenario, n, spve, tr2) {
  rds <- file.path(out_dir, paste0("cell_", cell_id, ".rds"))
  if (file.exists(rds)) {
    rows <- readRDS(rds)
    if (length(rows) == reps && all(vapply(rows, is_ok, logical(1)))) {
      cat("SKIP:", cell_id, "\n")
      return(rows)
    }
  }
  t0 <- Sys.time()
  rows <- plapply(seq_len(reps), engine_rep, cell_id = cell_id,
                  scenario = scenario, n = n, spve = spve, tr2 = tr2)
  ok <- vapply(rows, is_ok, logical(1))
  cat(sprintf("CELL %s ok %d/%d wall %.0fs %s\n", cell_id, sum(ok),
              length(rows), difftime(Sys.time(), t0, units = "secs"),
              if (any(!ok)) paste(unique(vapply(
                rows[!ok], function(x) x$status, character(1))),
                collapse = ",") else ""))
  saveRDS(rows, rds)
  rows
}

## =========================== P1: secondary PVE grid =======================
p1_csv <- file.path(out_dir, "stage77_p1_pve.csv")
if (!file.exists(p1_csv)) {
  cells <- sprintf("p1_sp%.4f",
                   c(0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.025))
  p1 <- do.call(rbind, lapply(cells, function(cid) {
    sp <- as.numeric(sub("p1_sp", "", cid))
    rows <- run_cell(cid, "two_linked_trait_specific", 1000L, sp, 0.3)
    summarize_cell(rows, cid, extras = list(secondary_pve = sp))
  }))
  rownames(p1) <- NULL
  write.csv(p1, p1_csv, row.names = FALSE)
  print(p1[, c("secondary_pve", "n_ok", "full_recovery",
               "oracle_locus_recovery", "oracle_causal_set_recovery",
               "oracle_primary_Mregion", "attribution_cand_exact",
               "extra_signal_rate", "median_cond_Q")], row.names = FALSE)
}

## =========================== P2: sample size grid =========================
p2_csv <- file.path(out_dir, "stage77_p2_n.csv")
if (!file.exists(p2_csv)) {
  p2 <- do.call(rbind, lapply(c(500, 750, 1000, 1250, 1500, 2000), function(nn) {
    cid <- sprintf("p2_n%d", nn)
    rows <- run_cell(cid, "two_linked_trait_specific", nn, 0.0175, 0.3)
    summarize_cell(rows, cid, extras = list(n = nn))
  }))
  rownames(p2) <- NULL
  write.csv(p2, p2_csv, row.names = FALSE)
  print(p2[, c("n", "n_ok", "full_recovery", "oracle_causal_set_recovery",
               "oracle_primary_Mregion", "attribution_cand_exact",
               "extra_signal_rate", "mean_runtime")], row.names = FALSE)
}

## =========================== P3: LD grid ==================================
p3_csv <- file.path(out_dir, "stage77_p3_ld.csv")
if (!file.exists(p3_csv)) {
  p3 <- do.call(rbind, lapply(c(0.1, 0.3, 0.5, 0.7, 0.9), function(r2) {
    cid <- sprintf("p3_r2%.1f", r2)
    rows <- run_cell(cid, "two_linked_trait_specific", 1000L, 0.0175, r2)
    summarize_cell(rows, cid, extras = list(target_r2 = r2))
  }))
  rownames(p3) <- NULL
  write.csv(p3, p3_csv, row.names = FALSE)
  print(p3[, c("target_r2", "n_ok", "mean_empirical_r2", "contam_marginal",
               "contam_conditional", "full_recovery",
               "oracle_primary_Mregion")], row.names = FALSE)
}

## =========================== P4: rho/tau boundary =========================
p4_csv <- file.path(out_dir, "stage77_p4_rho.csv")
if (!file.exists(p4_csv)) {
  p4_rep <- function(rep_id, rho) {
    tryCatch({
      seed <- CondPED:::.seed_for_rep(master_seed,
                                      sprintf("p4_rho%.2f", rho), rep_id)
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
      rep_exact <- NA_real_
      cand_exact <- NA_real_
      irr_rate <- NA_real_
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
            signals = est, candidate_sets = pipe$candidate_sets,
            minimum_representative_sets =
              pipe$subset_analysis$minimum_representative_sets,
            irreducible_modules =
              pipe$subset_analysis$irreducible_modules,
            subset_table = pipe$subset_analysis$subset_table,
            positions = position
          )
          ev <- tryCatch(
            CondPED:::.evaluate_one_replicate(
              sim$truth, estimates,
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
            irr_rate <- ev$irr_family_recovery
          }
        }
      }
      list(rep_id = rep_id, status = "ok", rho = rho,
           rep_exact = rep_exact, cand_exact = cand_exact,
           irr_rate = irr_rate)
    }, error = function(e) {
      list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
    })
  }
  p4 <- do.call(rbind, lapply(c(0.00, 0.02, 0.05, 0.08), function(rho) {
    cid <- sprintf("p4_rho%.2f", rho)
    rds <- file.path(out_dir, paste0("cell_", cid, ".rds"))
    if (file.exists(rds)) {
      rows <- readRDS(rds)
    } else {
      t0 <- Sys.time()
      rows <- plapply(seq_len(reps), p4_rep, rho = rho)
      ok <- vapply(rows, is_ok, logical(1))
      cat(sprintf("CELL %s ok %d/%d wall %.0fs\n", cid, sum(ok),
                  length(rows), difftime(Sys.time(), t0, units = "secs")))
      saveRDS(rows, rds)
    }
    rows <- rows[vapply(rows, is_ok, logical(1))]
    rv <- vapply(rows, function(x) x$rep_exact, numeric(1))
    cv <- vapply(rows, function(x) x$cand_exact, numeric(1))
    iv <- vapply(rows, function(x) x$irr_rate, numeric(1))
    data.frame(rho = rho, n_ok = length(rows), n_eval = sum(!is.na(rv)),
               unresolved_rate = mean(is.na(rv)),
               rep_exact = mean(rv[!is.na(rv)] >= 1 - 1e-8),
               cand_exact = mean(cv[!is.na(cv)] >= 1 - 1e-8),
               irr_exact = mean(iv[!is.na(iv)] >= 1 - 1e-8))
  }))
  rownames(p4) <- NULL
  write.csv(p4, p4_csv, row.names = FALSE)
  print(p4, row.names = FALSE)
}

## =========================== P5: q=2 vs q=3 ===============================
p5_csv <- file.path(out_dir, "stage77_p5_q3.csv")
if (!file.exists(p5_csv)) {
  p5 <- do.call(rbind, lapply(
    c("two_linked_trait_specific", "three_linked_trait_specific"),
    function(sc) {
      cid <- sprintf("p5_%s", sc)
      rows <- run_cell(cid, sc, 1000L, 0.0175, 0.3)
      summarize_cell(rows, cid, extras = list(scenario = sc))
    }))
  rownames(p5) <- NULL
  write.csv(p5, p5_csv, row.names = FALSE)
  print(p5[, c("scenario", "n_ok", "full_recovery", "full_count_exact",
               "oracle_causal_set_recovery", "oracle_primary_Mregion",
               "attribution_cand_exact", "extra_signal_rate",
               "mean_runtime")], row.names = FALSE)
}

cat("\nSTAGE77 CALIBRATION DONE\n")
