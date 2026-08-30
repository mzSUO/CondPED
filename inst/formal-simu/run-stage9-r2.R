## Stage 9 batch 3: LD sensitivity (freeze §1.57) with E3 conditional
## effect contamination gradient (freeze §1.27): contamination MUST be
## monotone in r2. linked_pseudo_multitrait at target_r2 in {0.1, 0.5, 0.7}
## (main r2=0.3 covered by Stage 8.4), R=200/cell, paired, workers=4.
## Contamination requires a per-rep refit (deterministic regen + fit).
## Launch: bash run_sim.sh stage9-r2 inst/formal-simu/run-stage9-r2.R
devtools::load_all(quiet = TRUE)
source("inst/formal-simu/lib-stage9.R")

out_dir <- "inst/formal-simu/output/stage9-r2"
grid <- do.call(rbind, lapply(c(0.1, 0.5, 0.7), function(r2) {
  data.frame(
    experiment = "end_to_end", scenario = "linked_pseudo_multitrait",
    analysis_mode = "full", comparison_pipeline = "paired",
    n = 1000L, m = 4L, p = 1000L,
    locus_pve = 0.02, secondary_signal_pve = 0.02,
    local_ld = NA_character_, target_r2 = r2,
    correlation = "block", target_loss = NA_real_,
    tolerance = 0.10, alpha_omnibus = 0.05,
    stringsAsFactors = FALSE
  )
}))

final <- stage9_run_grid(grid, out_dir, reps_total = 200L, workers = 4L)

## ---- contamination per rep (regen + refit) ------------------------------------
trait_names <- paste0("Trait", 1:4)
grid$key <- vapply(seq_len(nrow(grid)), function(i) stage9_key(grid[i, ]),
                   character(1))
final$dir <- basename(dirname(final$file))
final$target_r2 <- grid$target_r2[match(final$dir, grid$key)]

contam_rep <- function(f) {
  o <- readRDS(f)
  if (is.null(o) || !isTRUE(o$status$ok)) return(NULL)
  s <- o$settings
  args <- list(n = s$n, m = s$m, p = s$p, experiment = s$experiment,
               scenario = s$architecture, locus_pve = s$locus_pve,
               secondary_signal_pve = s$secondary_signal_pve,
               correlation = s$correlation, tolerance = s$tolerance,
               target_r2 = s$target_r2)
  set.seed(o$seed)
  sim <- do.call(simulate_condped_data, args)
  G <- sim$G
  causals <- sim$truth$causal_markers
  fit <- fit_mt_null(sim$Y, K = sim$K_bg, control = list(maxit = 300L))
  if (!isTRUE(fit$status$ok)) return(NULL)
  eff1 <- estimate_mt_effects(fit, G, targets = causals[1])
  b1 <- stats::setNames(eff1$effects_long$beta,
                        eff1$effects_long$trait)[trait_names]
  X_C <- G[, causals[2], drop = FALSE]
  proj <- CondPED:::.build_conditional_projection(fit, X_C)
  cs <- CondPED:::.conditional_mt_scan(
    proj, G[, causals[1], drop = FALSE], marker_ids = causals[1],
    return_effects = TRUE)
  ccond <- if (!is.null(cs$effects)) {
    unname(cs$effects$beta[1, "Trait2"])
  } else NA_real_
  data.frame(
    target_r2 = s$target_r2, rep_id = o$rep_id,
    r2_empirical = as.numeric(stats::cor(G[, causals[1]], G[, causals[2]])^2),
    contam_marginal = unname(b1["Trait2"]),
    contam_conditional = ccond,
    stringsAsFactors = FALSE
  )
}
files <- final$file[final$status == "ok"]
rows <- parallel::mclapply(files, contam_rep, mc.cores = 4L)
rows <- rows[!vapply(rows, is.null, logical(1))]
d <- do.call(rbind, rows)
write.csv(d, file.path(out_dir, "contamination.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(split(d, d$target_r2), function(x) {
  data.frame(
    target_r2 = x$target_r2[1], n = nrow(x),
    mean_empirical_r2 = mean(x$r2_empirical),
    contam_marginal = mean(abs(x$contam_marginal)),
    contam_conditional = mean(abs(x$contam_conditional), na.rm = TRUE),
    reduction = mean(abs(x$contam_marginal) - abs(x$contam_conditional),
                     na.rm = TRUE),
    stringsAsFactors = FALSE)
}))
## append the Stage 8.4 main point (r2 = 0.3) for the gradient check
main_pt <- data.frame(target_r2 = 0.3, n = 500L, mean_empirical_r2 = 0.3015,
                      contam_marginal = 0.1369, contam_conditional = 0.0455,
                      reduction = 0.0913)
summ_all <- rbind(summ, main_pt)
summ_all <- summ_all[order(summ_all$target_r2), ]
mono <- all(diff(summ_all$contam_marginal) >= 0)
write.csv(summ_all, file.path(out_dir, "summary.csv"), row.names = FALSE)

lines <- c(
  "# Stage 9 batch 3: LD sensitivity (freeze §1.57) + E3 contamination gradient (§1.27)",
  "",
  sprintf("- target_r2 in {0.1, 0.5, 0.7} (R=200/cell) + main 0.3 (Stage 8.4, R=500); %d/%d ok",
          sum(final$status == "ok"), nrow(final)),
  "",
  "| target r2 | empirical r2 | |contam| marginal | |contam| conditional | reduction |",
  "|---|---|---|---|---|",
  apply(summ_all, 1, function(r) sprintf("| %.1f | %.3f | %.4f | %.4f | %.4f |",
    as.numeric(r[["target_r2"]]), as.numeric(r[["mean_empirical_r2"]]),
    as.numeric(r[["contam_marginal"]]), as.numeric(r[["contam_conditional"]]),
    as.numeric(r[["reduction"]]))),
  "",
  sprintf("- **monotone in r2: %s**", if (mono) "YES" else "NO — INVESTIGATE"),
  "- 核对点（§1.27）：marginal contamination 随 r2 单调上升，conditional 显著更低"
)
writeLines(lines, file.path(out_dir, "report.md"))
print(summ_all, row.names = FALSE)
cat("BATCH3 DONE  monotone:", mono, "\n")
