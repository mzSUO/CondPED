# Validation for simulate_condped_data() (v1.0 Stage 4).
#
# NOT part of testthat. Per Runbook 6.8: each architecture is
# generated 20 times and ONLY the population truth is checked (no
# CondPED fitting). Reports per architecture:
#   attempts, PSD failures, PVE error, truth margin,
#   candidate set, minimum representative sets, irreducible modules.
#
# Run from the package root:
#   Rscript inst/validation/validate-simulate-condped-data.R

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
tol <- 0.10
architectures <- c("null", "candidate_single", "candidate_pair",
                   "candidate_dense", "representative_singleton",
                   "representative_pair", "representative_full",
                   "irreducible_singleton", "irreducible_pair",
                   "multiple_modules")

loss_of <- function(tab, key) {
  i <- match(key, tab$representing_key)
  if (is.na(i)) NA_real_ else tab$representation_loss[i]
}
keys_at <- function(df) sort(df$trait_key[df$tolerance == tol])

summ <- data.frame(
  architecture = character(), n_accept = integer(),
  median_attempts = numeric(), psd_failures = integer(),
  max_pve_error = numeric(), min_margin = numeric(),
  stringsAsFactors = FALSE
)

for (ai in seq_along(architectures)) {
  arch <- architectures[ai]
  cat(sprintf("\n== %s ==\n", arch))
  n_accept <- 0L
  atts <- numeric(n_rep)
  psd <- 0L
  max_pve_err <- 0
  min_margin <- Inf
  struct_ok <- 0L
  for (r in seq_len(n_rep)) {
    s <- simulate_condped_data(
      n = 200, m = 4L, p = 200L, architecture = arch,
      correlation = "block", seed = 51000 + 1000 * ai + r
    )
    atts[r] <- s$diagnostics$attempts
    psd <- psd + s$diagnostics$psd_failures
    if (!isTRUE(s$status$ok)) next
    n_accept <- n_accept + 1L
    cand <- s$truth$candidate_traits
    tab <- s$truth$subset_table
    reps <- keys_at(s$truth$minimum_representative_sets)
    mods <- keys_at(s$truth$irreducible_modules)
    if (length(cand) > 0L) {
      max_pve_err <- max(max_pve_err,
                         abs(mean(s$truth$locus_pve[cand]) - 0.01))
      min_margin <- min(min_margin,
                        min(abs(tab$representation_loss - tol)))
    }
    ok <- switch(arch,
      null = length(cand) == 0L && nrow(tab) == 0L,
      candidate_single = identical(cand, "Trait1"),
      candidate_pair = identical(cand, c("Trait1", "Trait2")),
      candidate_dense = length(cand) == 4L,
      representative_singleton = reps[1] == "Trait1" &&
        s$truth$minimum_representative_sets$set_size[
          s$truth$minimum_representative_sets$tolerance == tol][1] == 1L,
      representative_pair = "Trait1|Trait2" %in% reps &&
        all(tab$representation_loss[tab$set_size == 1L] > tol),
      representative_full = all(
        tab$representation_loss[tab$set_size < 4L] > tol),
      irreducible_singleton = "Trait1" %in% mods,
      irreducible_pair = "Trait1|Trait2" %in% mods &&
        loss_of(tab, "Trait2|Trait3|Trait4") <= tol &&
        loss_of(tab, "Trait1|Trait3|Trait4") <= tol,
      multiple_modules = identical(mods, c("Trait1", "Trait3|Trait4")),
      FALSE
    )
    struct_ok <- struct_ok + isTRUE(ok)
    if (r <= 2L) {
      cat(sprintf(paste0("  rep %2d: attempts=%d psd_fail=%d pve_err=%.2e ",
                         "margin=%.3f cand=[%s] reps=[%s] mods=[%s]\n"),
                  r, s$diagnostics$attempts, s$diagnostics$psd_failures,
                  if (length(cand)) abs(mean(s$truth$locus_pve[cand]) - 0.01) else 0,
                  if (nrow(tab)) min(abs(tab$representation_loss - tol)) else NA,
                  paste(cand, collapse = ","),
                  paste(reps, collapse = "/"), paste(mods, collapse = "/")))
    }
  }
  summ <- rbind(summ, data.frame(
    architecture = arch, n_accept = n_accept,
    median_attempts = median(atts), psd_failures = psd,
    max_pve_error = max_pve_err,
    min_margin = if (is.finite(min_margin)) min_margin else NA_real_
  ))
  cat(sprintf("accepted %d/%d; structure rule ok %d/%d\n",
              n_accept, n_rep, struct_ok, n_accept))
  report(n_accept == n_rep,
         sprintf("%s: all %d generations accepted", arch, n_rep))
  report(struct_ok == n_accept,
         sprintf("%s: population truth satisfies its definition", arch))
}

cat("\n== summary ==\n")
print(summ, right = FALSE)
report(all(summ$max_pve_error[!is.na(summ$max_pve_error)] < 1e-8 |
             is.na(summ$max_pve_error)),
       "PVE realised within 1e-8 of target for all accepted generations")
struct_rows <- summ$architecture %in% architectures[5:10]
report(all(summ$min_margin[struct_rows] >= 0.02),
       "truth_margin >= 0.02 achieved in all structural scenes")

cat("\n")
if (length(failures) > 0L) {
  cat("VALIDATION FAILURES:\n")
  for (f in failures) cat(" -", f, "\n")
  quit(status = 1L)
} else {
  cat("All simulate_condped_data validation checks passed.\n")
}
