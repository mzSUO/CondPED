## Stage 7.4 C: oracle-primary PVE x LD calibration for
## two_linked_trait_specific (post-geometry-fix).
## For each (secondary_signal_pve, target_r2) cell, condition on the TRUE
## primary causal and test the TRUE secondary causal with the conditional
## scan; within-locus Bonferroni uses M = (#region markers - 1), matching
## the resolver's frozen convention. No CondPED statistic is modified.
## Resumable: cells with an existing ok RDS are skipped.
devtools::load_all(quiet = TRUE)

master_seed <- 20260816L
reps <- 100L
workers <- as.integer(Sys.getenv("STAGE74_WORKERS", "16"))
out_dir <- "inst/validation/output/stage74/oracle_calibration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenario <- "two_linked_trait_specific"
experiment <- "signal_resolution"
locus_pve <- 0.02
alpha_signal <- 0.05

grid <- expand.grid(
  secondary_signal_pve = c(0.005, 0.0075, 0.01, 0.015, 0.0175, 0.02),
  target_r2 = c(0, 0.3, 0.6)
)

one_rep <- function(rep_id, spve, tr2) {
  tryCatch({
    psid <- sprintf("%s|sp%.4f_r2%.1f", scenario, spve, tr2)
    seed <- CondPED:::.seed_for_rep(master_seed, psid, rep_id)
    set.seed(seed)
    sim <- simulate_condped_data(
      n = 1000L, m = 4L, p = 1000L,
      experiment = experiment, scenario = scenario,
      locus_pve = locus_pve, secondary_signal_pve = spve,
      target_r2 = tr2, correlation = "block", tolerance = 0.10
    )
    if (!isTRUE(sim$status$ok)) {
      return(list(rep_id = rep_id, status = "sim_failed"))
    }
    G <- sim$G
    causals <- sim$truth$causal_markers
    fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
    if (!isTRUE(fit$status$ok)) {
      return(list(rep_id = rep_id, status = "fit_failed"))
    }
    X_C <- G[, causals[1], drop = FALSE]
    proj <- CondPED:::.build_conditional_projection(fit, X_C)
    cs <- CondPED:::.conditional_mt_scan(
      proj, G[, causals[2], drop = FALSE],
      marker_ids = causals[2], return_effects = FALSE
    )
    ctab <- cs$conditional
    tl <- as.data.frame(sim$truth$loci)[1, ]
    marker_ids <- colnames(G)
    position <- stats::setNames(seq_along(marker_ids) * 1000, marker_ids)
    M_remaining <- sum(position >= tl$start & position <= tl$end) - 1L
    p_adj <- min(1, M_remaining * ctab$p_value)
    list(
      rep_id = rep_id, status = "ok",
      r2_causal = as.numeric(stats::cor(G[, causals[1]], G[, causals[2]])^2),
      causal_dist = abs(diff(position[causals])),
      Q = ctab$Q, df = ctab$df, p_raw = ctab$p_value, p_adj = p_adj,
      scan_status = ctab$status,
      pass = isTRUE(p_adj <= alpha_signal)
    )
  }, error = function(e) {
    list(rep_id = rep_id, status = paste("error:", conditionMessage(e)))
  })
}

is_ok <- function(x) {
  is.list(x) && !is.null(x$status) && is.character(x$status) &&
    length(x$status) == 1L && x$status == "ok"
}

summ_rows <- list()
t0_all <- Sys.time()
for (i in seq_len(nrow(grid))) {
  spve <- grid$secondary_signal_pve[i]
  tr2 <- grid$target_r2[i]
  tag <- sprintf("sp%.4f_r2%.1f", spve, tr2)
  rds <- file.path(out_dir, paste0("cell_", tag, ".rds"))
  if (file.exists(rds)) {
    prev <- readRDS(rds)
    if (length(prev) == reps && all(vapply(prev, is_ok, logical(1)))) {
      cat("SKIP (done):", tag, "\n")
      out <- prev
      summ_rows[[tag]] <- data.frame(
        secondary_signal_pve = spve, target_r2 = tr2,
        n_ok = sum(vapply(out, is_ok, logical(1))),
        oracle_primary_power = mean(vapply(out, function(x) x$pass, logical(1))),
        mean_r2_causal = mean(vapply(out, function(x) x$r2_causal, numeric(1))),
        median_causal_dist = stats::median(vapply(out, function(x) x$causal_dist, numeric(1)))
      )
      next
    }
  }
  t0 <- Sys.time()
  if (.Platform$OS.type == "unix") {
    out <- parallel::mclapply(seq_len(reps), one_rep, spve = spve, tr2 = tr2,
                              mc.cores = workers)
  } else {
    out <- lapply(seq_len(reps), one_rep, spve = spve, tr2 = tr2)
  }
  ok <- vapply(out, is_ok, logical(1))
  cat(sprintf("CELL %s  ok %d/%d  wall %.1fs  status: %s\n", tag,
              sum(ok), length(out),
              difftime(Sys.time(), t0, units = "secs"),
              paste(names(table(vapply(out[!ok], function(x) x$status, character(1)))),
                    collapse=",")))
  saveRDS(out, rds)
  rows <- out[ok]
  summ_rows[[tag]] <- data.frame(
    secondary_signal_pve = spve, target_r2 = tr2,
    n_ok = length(rows),
    oracle_primary_power = mean(vapply(rows, function(x) x$pass, logical(1))),
    mean_r2_causal = mean(vapply(rows, function(x) x$r2_causal, numeric(1))),
    median_causal_dist = stats::median(vapply(rows, function(x) x$causal_dist, numeric(1)))
  )
}

summ <- do.call(rbind, summ_rows)
rownames(summ) <- NULL
print(summ, row.names = FALSE)
write.csv(summ, file.path(out_dir, "calibration_summary.csv"),
          row.names = FALSE)
cat("WALL_SECONDS:", round(difftime(Sys.time(), t0_all, units = "secs"), 1), "\n")
cat("STAGE74 ORACLE CALIBRATION DONE\n")
