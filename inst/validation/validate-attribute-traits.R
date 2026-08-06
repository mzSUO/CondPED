# Heavy statistical validation for attribute_traits() (v1.0 Stage 3).
#
# NOT part of testthat. Pilot scale (30 replicates per configuration):
#   1. null architecture: false-candidate rates;
#   2. effect architectures: candidate-set TPR / FDP / exact recovery /
#      Jaccard against the truth;
#   3. empirical per-locus FWER of the within-locus Holm procedure;
#   4. all_traits oracle mode: candidate set == all traits at selected
#      loci, identical selected sets;
#   5. failed replicates are kept in the results table / RDS.
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
alpha_trait <- 0.05
out_dir <- file.path("inst", "validation", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
trait_names <- paste0("Trait", seq_len(m_tr))

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
    att <- attribute_traits(scan, est, omnibus_method = "BH",
                            candidate_mode = "holm_fwer",
                            alpha_trait = alpha_trait)
    att_all <- attribute_traits(scan, est, omnibus_method = "BH",
                                candidate_mode = "all_traits")
    qtl_marker <- colnames(sim$G)[sim$qtl_index]
    true_traits <- trait_names[sim$truth$A[[1L]]]
    list(ok = TRUE, qtl_marker = qtl_marker, true_traits = true_traits,
         att = att, att_all = att_all)
  }, error = function(e) list(ok = FALSE, status = conditionMessage(e)))
}

set_metrics <- function(cand, true_traits) {
  tp <- length(intersect(cand, true_traits))
  fp <- length(setdiff(cand, true_traits))
  list(
    tpr = if (length(true_traits) > 0L) tp / length(true_traits) else NA_real_,
    fdp = if (length(cand) > 0L) fp / length(cand) else 0,
    esr = identical(sort(cand), sort(true_traits)),
    jaccard = {
      u <- length(union(cand, true_traits))
      if (u > 0L) tp / u else NA_real_
    },
    n_cand = length(cand),
    n_false = fp
  )
}

configs <- c("null", "single_trait", "shared_same", "shared_opposite", "dense")
all_out <- list()

for (ci in seq_along(configs)) {
  cfg <- configs[ci]
  cat(sprintf("\n== %s ==\n", cfg))
  tab <- data.frame(
    replicate = seq_len(n_rep), ok = NA_integer_, status = NA_character_,
    qtl_selected = NA_integer_, tpr = NA_real_, fdp = NA_real_,
    esr = NA_integer_, jaccard = NA_real_, n_false = NA_integer_,
    same_selection = NA_integer_, oracle_complete = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (r in seq_len(n_rep)) {
    out <- run_one(cfg, seed = 41000 + 1000 * ci + r)
    tab$ok[r] <- as.integer(isTRUE(out$ok))
    tab$status[r] <- if (isTRUE(out$ok)) "ok" else out$status
    if (!isTRUE(out$ok)) next
    att <- out$att
    cand <- if (out$qtl_marker %in% names(att$candidate_sets)) {
      att$candidate_sets[[out$qtl_marker]]
    } else {
      character()
    }
    m <- set_metrics(cand, out$true_traits)
    tab$qtl_selected[r] <- as.integer(out$qtl_marker %in% att$selected_loci)
    tab$tpr[r] <- m$tpr; tab$fdp[r] <- m$fdp
    tab$esr[r] <- as.integer(m$esr); tab$jaccard[r] <- m$jaccard
    tab$n_false[r] <- m$n_false
    # oracle mode: same layer-1 selection, candidates == all traits
    tab$same_selection[r] <- identical(sort(att$selected_loci),
                                       sort(out$att_all$selected_loci))
    ca <- out$att_all$candidate_sets[[out$qtl_marker]]
    tab$oracle_complete[r] <-
      isTRUE(setequal(ca, trait_names)) ||
      !out$qtl_marker %in% names(out$att_all$candidate_sets) &&
      !out$qtl_marker %in% out$att_all$selected_loci
  }
  all_out[[cfg]] <- tab
  okr <- which(tab$ok == 1L)
  cat(sprintf("replicates ok: %d/%d; QTL selected: %.2f\n",
              length(okr), n_rep, mean(tab$qtl_selected[okr])))
  cat(sprintf("Holm: TPR %.3f | FDP %.3f | ESR %.3f | Jaccard %.3f\n",
              mean(tab$tpr[okr], na.rm = TRUE), mean(tab$fdp[okr]),
              mean(tab$esr[okr]), mean(tab$jaccard[okr], na.rm = TRUE)))
  cat(sprintf("per-locus FWER (any false candidate at QTL): %.3f\n",
              mean(tab$n_false[okr] > 0)))
  report(length(okr) / n_rep >= 0.90,
         sprintf("%s: at least 90%% replicates completed", cfg))
  report(all(tab$same_selection[okr] == 1L),
         sprintf("%s: all_traits selects identical loci", cfg))
  report(all(tab$oracle_complete[okr] == 1L),
         sprintf("%s: all_traits candidate set == all traits", cfg))
}

# ---- global error rates -------------------------------------------------------
cat("\n== global error rates ==\n")
null_tab <- all_out[["null"]]
null_ok <- which(null_tab$ok == 1L)
null_fw <- mean(null_tab$n_false[null_ok] > 0)
cat(sprintf("null: false-candidate rate = %.3f (nominal <= 0.05)\n", null_fw))
report(null_fw <= 0.15, "null: false-candidate rate <= 0.15 (pilot)")

eff_cfgs <- c("single_trait", "shared_same", "shared_opposite", "dense")
fwer <- unlist(lapply(eff_cfgs, function(cfg) {
  tab <- all_out[[cfg]]
  (tab$n_false[tab$ok == 1L] > 0)
}))
emp_fwer <- mean(fwer)
cat(sprintf("within-locus Holm empirical per-locus FWER = %.3f (target %.2f)\n",
            emp_fwer, alpha_trait))
report(emp_fwer <= alpha_trait + 0.08,
       "Holm empirical FWER <= alpha_trait + 0.08 (pilot tolerance)")

rds_file <- file.path(out_dir, "attribute-traits-validation.rds")
saveRDS(list(config = list(n_rep = n_rep, n_ind = n_ind, m = m_tr,
                           p = p_snp, locus_pve = locus_pve,
                           alpha_trait = alpha_trait),
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
