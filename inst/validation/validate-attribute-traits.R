# Heavy statistical validation for attribute_traits().
#
# NOT part of testthat. Heavy acceptance per contract section 4.5 and
# Methods sections 3.7 / 5.7:
#   1. null and projection_induced architectures: false-attribution rates;
#   2. effect architectures: trait-level TPR, FDP, exact set recovery
#      (ESR) and Jaccard of the attributed sets A;
#   3. empirical FDR of the BB mode;
#   4. empirical FWER of the Holm mode;
#   5. MT-Posthoc (mode = "none") differs from BB only through the
#      correction: identical selected sets, BB attributions are a subset;
#   6. failed replicates are kept in the results table / RDS.
#
# Pilot scale (30 replicates per configuration). The formal 300-500
# replicate run uses the same script with n_rep raised; it is not
# executed by default.
#
# Run from the package root:
#   Rscript inst/validation/validate-attribute-traits.R

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(CondPED)
}

failures <- character()
report <- function(ok, label) {
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) failures <<- c(failures, label)
}

n_rep <- 30L
n_ind <- 300L
m_tr  <- 3L
p_snp <- 400L
locus_pve <- 0.03
q_target <- 0.05
alpha_total <- 0.05
out_dir <- file.path("inst", "validation", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

trait_names <- paste0("Trait", seq_len(m_tr))

# Per-replicate pipeline: fit -> scan -> estimate at QTL -> attribute.
run_one <- function(arch, seed) {
  tryCatch({
    sim <- simulate_condped_data(n = n_ind, m = m_tr, p = p_snp,
                                 architecture = arch, n_qtl = 1L,
                                 locus_pve = locus_pve, seed = seed)
    fit <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 2L,
                       control = list(maxit = 300))
    if (!isTRUE(fit$status$ok)) return(list(ok = FALSE, status = "fit_failed"))
    scan <- scan_mt_omnibus(fit, sim$G)
    est <- estimate_mt_effects(fit, sim$G, loci = sim$qtl_index)
    att_bb <- attribute_traits(scan, est, omnibus_method = "BH",
                               attribution_mode = "bb_fdr",
                               q_target = q_target)
    att_holm <- attribute_traits(scan, est, omnibus_method = "bonferroni",
                                 attribution_mode = "holm_fwer",
                                 alpha_total = alpha_total,
                                 alpha_split = c(0.025, 0.025))
    att_none <- attribute_traits(scan, est, omnibus_method = "BH",
                                 attribution_mode = "none",
                                 q_target = q_target)
    qtl_marker <- colnames(sim$G)[sim$qtl_index]
    true_traits <- trait_names[sim$truth$A[[1L]]]
    list(ok = TRUE, qtl_marker = qtl_marker, true_traits = true_traits,
         bb = att_bb, holm = att_holm, none = att_none)
  }, error = function(e) list(ok = FALSE, status = conditionMessage(e)))
}

# Set metrics of an attributed set against the truth at the QTL.
set_metrics <- function(att_traits, true_traits) {
  att_traits <- intersect(att_traits, trait_names)
  tp <- length(intersect(att_traits, true_traits))
  fp <- length(setdiff(att_traits, true_traits))
  list(
    tpr = if (length(true_traits) > 0L) tp / length(true_traits) else NA_real_,
    fdp = if (length(att_traits) > 0L) fp / length(att_traits) else 0,
    esr = identical(sort(att_traits), sort(true_traits)),
    jaccard = {
      u <- length(union(att_traits, true_traits))
      if (u > 0L) tp / u else NA_real_
    },
    n_att = length(att_traits),
    target_attributed = "Trait1" %in% att_traits
  )
}

get_att <- function(att, marker) {
  if (marker %in% names(att$A)) att$A[[marker]] else character()
}

configs <- c("null", "projection_induced", "single_trait",
             "shared_same", "shared_opposite", "dense")
all_out <- list()
metrics_pool <- list()

for (ci in seq_along(configs)) {
  cfg <- configs[ci]
  cat(sprintf("\n== %s ==\n", cfg))
  tab <- data.frame(
    replicate = seq_len(n_rep), ok = NA_integer_, status = NA_character_,
    qtl_selected = NA_integer_,
    tpr_bb = NA_real_, fdp_bb = NA_real_, esr_bb = NA_integer_,
    jac_bb = NA_real_,
    fdp_holm = NA_real_, fw_holm = NA_integer_,
    fdp_none = NA_real_,
    target_att_bb = NA_integer_,
    subset_bb_none = NA_integer_, same_S_bb_none = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (r in seq_len(n_rep)) {
    out <- run_one(cfg, seed = 31000 + 1000 * ci + r)
    tab$ok[r] <- as.integer(isTRUE(out$ok))
    tab$status[r] <- if (isTRUE(out$ok)) "ok" else out$status
    if (!isTRUE(out$ok)) next
    bb <- out$bb; holm <- out$holm; none <- out$none
    tab$qtl_selected[r] <- as.integer(out$qtl_marker %in% bb$selected_loci)
    m_bb <- set_metrics(get_att(bb, out$qtl_marker), out$true_traits)
    m_holm <- set_metrics(get_att(holm, out$qtl_marker), out$true_traits)
    m_none <- set_metrics(get_att(none, out$qtl_marker), out$true_traits)
    tab$tpr_bb[r] <- m_bb$tpr
    tab$fdp_bb[r] <- m_bb$fdp
    tab$esr_bb[r] <- as.integer(m_bb$esr)
    tab$jac_bb[r] <- m_bb$jaccard
    tab$fdp_holm[r] <- m_holm$fdp
    tab$fw_holm[r] <- as.integer(m_holm$n_att > length(out$true_traits) ||
      any(!get_att(holm, out$qtl_marker) %in% out$true_traits))
    tab$fdp_none[r] <- m_none$fdp
    tab$target_att_bb[r] <- as.integer(m_bb$target_attributed)
    # MT-Posthoc must differ from BB only through the correction:
    # identical first-layer selection, BB attributions a subset.
    tab$same_S_bb_none[r] <-
      identical(sort(bb$selected_loci), sort(none$selected_loci))
    att_all <- function(att) {
      paste(rep(names(att$A), lengths(att$A)), unlist(att$A))
    }
    tab$subset_bb_none[r] <- all(att_all(bb) %in% att_all(none))
  }
  all_out[[cfg]] <- tab
  ok_rows <- which(tab$ok == 1L)
  cat(sprintf("replicates ok: %d/%d; QTL layer-1 selection rate: %.3f\n",
              length(ok_rows), n_rep,
              mean(tab$qtl_selected[ok_rows])))
  cat(sprintf(paste0(
    "BB    : TPR %.3f | FDP %.3f | ESR %.3f | Jaccard %.3f | ",
    "target-trait attribution %.3f\n"),
    mean(tab$tpr_bb[ok_rows], na.rm = TRUE),
    mean(tab$fdp_bb[ok_rows]),
    mean(tab$esr_bb[ok_rows]),
    mean(tab$jac_bb[ok_rows], na.rm = TRUE),
    mean(tab$target_att_bb[ok_rows])))
  cat(sprintf("Holm  : FDP %.3f | family-wise false attribution %.3f\n",
              mean(tab$fdp_holm[ok_rows]), mean(tab$fw_holm[ok_rows])))
  cat(sprintf("Posthoc (none): FDP %.3f\n", mean(tab$fdp_none[ok_rows])))
  cat(sprintf("identical S (BB vs Posthoc): %d/%d; BB subset of Posthoc: %d/%d\n",
              sum(tab$same_S_bb_none[ok_rows]), length(ok_rows),
              sum(tab$subset_bb_none[ok_rows]), length(ok_rows)))
  report(length(ok_rows) / n_rep >= 0.90,
         sprintf("%s: at least 90%% replicates completed", cfg))
  report(all(tab$same_S_bb_none[ok_rows] == 1L),
         sprintf("%s: BB and MT-Posthoc select identical loci", cfg))
  report(all(tab$subset_bb_none[ok_rows] == 1L),
         sprintf("%s: BB attributions are a subset of MT-Posthoc", cfg))
  metrics_pool[[cfg]] <- tab
}

# ---- global error-rate summaries ----------------------------------------------
cat("\n== global error rates ==\n")

# 1. null: any attribution anywhere is false. Count attributed rows.
null_tab <- metrics_pool[["null"]]
null_ok <- which(null_tab$ok == 1L)
null_fw <- mean(null_tab$fdp_bb[null_ok] > 0)  # any attribution => FDP 1
cat(sprintf("null: family-wise false attribution rate (BB) = %.3f\n", null_fw))
report(null_fw <= 0.15,
       "null: BB family-wise false attribution rate <= 0.15 (nominal 0.05)")

# 2. projection_induced: target trait has beta = 0.
pi_tab <- metrics_pool[["projection_induced"]]
pi_ok <- which(pi_tab$ok == 1L)
pi_rate <- mean(pi_tab$target_att_bb[pi_ok])
cat(sprintf(paste0(
  "projection_induced: target-trait (beta = 0) attribution rate = %.3f ",
  "(nominal ~ %.2f)\n"), pi_rate, q_target))
report(pi_rate <= 0.20,
       "projection_induced: target-trait attribution rate <= 0.20")

# 3. BB empirical FDR over effect architectures.
eff_cfgs <- c("single_trait", "shared_same", "shared_opposite", "dense")
fdp_all <- unlist(lapply(eff_cfgs, function(cfg) {
  tab <- metrics_pool[[cfg]]
  tab$fdp_bb[tab$ok == 1L]
}))
emp_fdr <- mean(fdp_all)
cat(sprintf("BB empirical FDR over effect architectures = %.3f (target %.2f)\n",
            emp_fdr, q_target))
report(emp_fdr <= q_target + 0.05,
       "BB empirical FDR <= q_target + 0.05 (pilot tolerance)")

# 4. Holm empirical FWER over effect architectures.
fw_all <- unlist(lapply(eff_cfgs, function(cfg) {
  tab <- metrics_pool[[cfg]]
  tab$fw_holm[tab$ok == 1L]
}))
emp_fwer <- mean(fw_all)
cat(sprintf("Holm empirical FWER over effect architectures = %.3f (target %.2f)\n",
            emp_fwer, alpha_total))
report(emp_fwer <= alpha_total + 0.07,
       "Holm empirical FWER <= alpha_total + 0.07 (pilot tolerance)")

rds_file <- file.path(out_dir, "attribute-traits-validation.rds")
saveRDS(list(config = list(n_rep = n_rep, n_ind = n_ind, m = m_tr, p = p_snp,
                           locus_pve = locus_pve, q_target = q_target,
                           alpha_total = alpha_total),
             configurations = all_out), rds_file)
cat(sprintf("\nAll replicates (failures included) saved to %s\n", rds_file))

cat("\n")
if (length(failures) > 0L) {
  cat("VALIDATION FAILURES:\n")
  for (f in failures) cat(" -", f, "\n")
  quit(status = 1L)
} else {
  cat("All attribute_traits validation checks passed.\n")
}
