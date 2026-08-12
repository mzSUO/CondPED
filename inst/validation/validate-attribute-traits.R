# Heavy validation for attribute_traits() (Stage 6B-5 signal level).
#
# NOT part of testthat. Runs the full condped(signal_mode="resolve")
# pipeline per replicate and compares the causal signal's candidate
# set against the simulation truth. Pilot scale: 20 replicates per
# configuration.
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

n_rep <- 20L
n_ind <- 300L
m_tr  <- 3L
p_snp <- 300L
locus_pve <- 0.04
alpha_trait <- 0.05
out_dir <- file.path("inst", "validation", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
trait_names <- paste0("Trait", seq_len(m_tr))

run_one <- function(arch, seed) {
  tryCatch({
    sim <- simulate_condped_data(n = n_ind, m = m_tr, p = p_snp,
                                 architecture = arch,
                                 locus_pve = locus_pve,
                                 correlation = "block", seed = seed)
    pos <- seq_len(p_snp) * 1000
    f <- condped(
      sim$Y, G = sim$G, K = sim$K_bg,
      position = pos,
      control = list(null_control = list(maxit = 300L),
                     window_bp = 5000)
    )
    causal_marker <- colnames(sim$G)[sim$causal_index]
    sig <- f$signals
    if (is.null(sig) || nrow(sig) == 0L) {
      return(list(ok = TRUE, found = FALSE, cand = character(),
                  true = sim$truth$candidate_traits,
                  n_signals = 0L))
    }
    hit <- which(sig$representative_snp == causal_marker)
    cand <- if (length(hit) > 0L) {
      sid <- sig$signal_id[hit[1L]]
      cs <- f$candidate_traits$candidate_sets[[sid]]
      if (is.null(cs)) character() else cs
    } else {
      character()
    }
    list(ok = TRUE, found = length(hit) > 0L, cand = cand,
         true = sim$truth$candidate_traits,
         n_signals = nrow(sig))
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
    n_false = fp
  )
}

configs <- c("null", "candidate_single", "candidate_pair",
             "candidate_dense")
all_out <- list()

for (ci in seq_along(configs)) {
  cfg <- configs[ci]
  cat(sprintf("\n== %s ==\n", cfg))
  tab <- data.frame(
    replicate = seq_len(n_rep), ok = NA_integer_, status = NA_character_,
    found = NA_integer_, tpr = NA_real_, fdp = NA_real_,
    esr = NA_integer_, jaccard = NA_real_, n_false = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (r in seq_len(n_rep)) {
    out <- run_one(cfg, seed = 61000 + 1000 * ci + r)
    tab$ok[r] <- as.integer(isTRUE(out$ok))
    tab$status[r] <- if (isTRUE(out$ok)) "ok" else out$status
    if (!isTRUE(out$ok)) next
    m <- set_metrics(out$cand, out$true)
    tab$found[r] <- as.integer(out$found)
    tab$tpr[r] <- m$tpr; tab$fdp[r] <- m$fdp
    tab$esr[r] <- as.integer(m$esr); tab$jaccard[r] <- m$jaccard
    tab$n_false[r] <- m$n_false
  }
  all_out[[cfg]] <- tab
  okr <- which(tab$ok == 1L)
  cat(sprintf(paste0(
    "replicates ok: %d/%d; signal found: %.2f | TPR %.3f | FDP %.3f | ",
    "ESR %.3f | Jaccard %.3f | false-candidate rate %.3f\n"),
    length(okr), n_rep, mean(tab$found[okr]),
    mean(tab$tpr[okr], na.rm = TRUE), mean(tab$fdp[okr]),
    mean(tab$esr[okr]), mean(tab$jaccard[okr], na.rm = TRUE),
    mean(tab$n_false[okr] > 0)))
  report(length(okr) / n_rep >= 0.90,
         sprintf("%s: at least 90%% replicates completed", cfg))
}

cat("\n== global error rates ==\n")
null_tab <- all_out[["null"]]
null_ok <- which(null_tab$ok == 1L)
null_fw <- mean(null_tab$n_false[null_ok] > 0)
cat(sprintf("null: false-candidate rate = %.3f (nominal <= 0.05)\n", null_fw))
report(null_fw <= 0.20, "null: false-candidate rate <= 0.20 (pilot)")

eff_cfgs <- c("candidate_single", "candidate_pair", "candidate_dense")
fwer <- unlist(lapply(eff_cfgs, function(cfg) {
  tab <- all_out[[cfg]]
  (tab$n_false[tab$ok == 1L] > 0)
}))
emp_fwer <- mean(fwer)
cat(sprintf("within-signal Holm empirical false-attribution rate = %.3f (target %.2f)\n",
            emp_fwer, alpha_trait))
report(emp_fwer <= alpha_trait + 0.10,
       "Holm empirical false-attribution rate <= alpha + 0.10 (pilot)")

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
