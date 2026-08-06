# Heavy statistical validation for scan_mt_omnibus().
#
# NOT part of testthat. Heavy acceptance per contract section 4.3:
#   1. Type I error under the null architecture;
#   2. QQ data and genomic-control lambda;
#   3. power under candidate_single / candidate_pair / candidate_dense;
#   4. rank-deficient and failure proportions;
#   5. all replicates (failures included) kept in the results table / RDS.
#
# Run from the package root:
#   Rscript inst/validation/validate-scan-mt-omnibus.R

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

## ---------------------------------------------------------------------------
## Configuration
## ---------------------------------------------------------------------------
n_rep   <- 30L
n_ind   <- 300L
m_tr    <- 2L
p_snp   <- 400L
alpha   <- 0.05
out_dir <- file.path("inst", "validation", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

run_one <- function(arch, seed) {
  tryCatch({
    sim <- simulate_condped_data(n = n_ind, m = m_tr, p = p_snp,
                                 architecture = arch,
                                 locus_pve = 0.03, seed = seed)
    fit <- fit_mt_null(sim$Y, K = sim$K_bg, n_starts = 2L,
                       control = list(maxit = 300))
    if (!isTRUE(fit$status$ok)) {
      return(list(ok = FALSE, status = "fit_failed"))
    }
    scan <- scan_mt_omnibus(fit, sim$G)
    keep <- scan$omnibus$status %in% c("ok", "rank_deficient")
    list(
      ok = TRUE,
      p_values = scan$omnibus$p_value[keep],
      Q = scan$omnibus$Q[keep],
      df = scan$omnibus$df[keep],
      qtl_p = scan$omnibus$p_value[sim$causal_index],
      n_rank_deficient = scan$diagnostics$n_rank_deficient,
      n_filtered = scan$diagnostics$n_filtered,
      converged = fit$convergence$code == 0
    )
  }, error = function(e) list(ok = FALSE, status = conditionMessage(e)))
}

collect <- function(arch, seed0) {
  res <- vector("list", n_rep)
  tab <- data.frame(
    replicate = seq_len(n_rep), ok = NA_integer_, converged = NA_integer_,
    n_rank_deficient = NA_integer_, n_filtered = NA_integer_,
    status = NA_character_, stringsAsFactors = FALSE
  )
  for (r in seq_len(n_rep)) {
    out <- run_one(arch, seed0 + r)
    res[[r]] <- out
    tab$ok[r] <- as.integer(isTRUE(out$ok))
    tab$status[r] <- if (isTRUE(out$ok)) "ok" else out$status
    if (isTRUE(out$ok)) {
      tab$converged[r] <- as.integer(out$converged)
      tab$n_rank_deficient[r] <- out$n_rank_deficient
      tab$n_filtered[r] <- out$n_filtered
    }
  }
  list(res = res, tab = tab)
}

## ---------------------------------------------------------------------------
## 1. Null architecture: type I error and lambda GC
## ---------------------------------------------------------------------------
cat("\n== 1. Null architecture: type I error / QQ / lambda GC ==\n")
null_out <- collect("null", 11000)
ok_res <- Filter(function(x) isTRUE(x$ok), null_out$res)
p_all <- unlist(lapply(ok_res, `[[`, "p_values"))
Q_all <- unlist(lapply(ok_res, `[[`, "Q"))
df_all <- unlist(lapply(ok_res, `[[`, "df"))
type1 <- mean(p_all < alpha)
se_t1 <- sqrt(type1 * (1 - type1) / length(p_all))
lambda_gc <- stats::median(Q_all) / stats::qchisq(0.5, df = stats::median(df_all))
cat(sprintf("replicates ok: %d/%d, markers tested: %d\n",
            length(ok_res), n_rep, length(p_all)))
cat(sprintf("empirical type I @ %.2f: %.4f (MC se %.4f)\n", alpha, type1, se_t1))
cat(sprintf("lambda GC (median): %.4f\n", lambda_gc))
cat(sprintf("rank-deficient markers: %d, filtered: %d\n",
            sum(null_out$tab$n_rank_deficient, na.rm = TRUE),
            sum(null_out$tab$n_filtered, na.rm = TRUE)))
# Monte Carlo tolerance for type I at 0.05 with ~30*400 obs: 0.05 +/- ~0.02
report(abs(type1 - alpha) < 3 * se_t1 + 0.01,
       "type I error within Monte Carlo tolerance")
report(lambda_gc > 0.90 && lambda_gc < 1.10,
       "genomic-control lambda in [0.90, 1.10]")
report(length(ok_res) / n_rep >= 0.90,
       "at least 90% of null replicates completed")

# QQ data (saved, not plotted here)
qq <- data.frame(
  expected = -log10(stats::ppoints(length(p_all))),
  observed = sort(-log10(p_all))
)

## ---------------------------------------------------------------------------
## 2. Power architectures
## ---------------------------------------------------------------------------
cat("\n== 2. Power under effect architectures ==\n")
power_tab <- data.frame(architecture = character(), power = numeric(),
                        n_ok = integer(), stringsAsFactors = FALSE)
for (arch in c("candidate_single", "candidate_pair", "candidate_dense")) {
  out <- collect(arch, 12000 + match(arch, c("candidate_single", "candidate_pair",
                                             "candidate_dense")) * 1000)
  ok_res <- Filter(function(x) isTRUE(x$ok), out$res)
  qtl_p <- vapply(ok_res, `[[`, numeric(1), "qtl_p")
  power <- mean(qtl_p < alpha, na.rm = TRUE)
  cat(sprintf("%-15s: power @ %.2f = %.3f (n_ok = %d/%d, qtl p non-NA = %d)\n",
              arch, alpha, power, length(ok_res), n_rep,
              sum(!is.na(qtl_p))))
  power_tab <- rbind(power_tab, data.frame(architecture = arch, power = power,
                                           n_ok = length(ok_res)))
  report(length(ok_res) / n_rep >= 0.90,
         sprintf("%s: at least 90%% replicates completed", arch))
}
report(all(power_tab$power[power_tab$architecture != "null"] > 0.5),
       "effect architectures have power > 0.5 at alpha = 0.05 (n = 300, PVE 3%)")

## ---------------------------------------------------------------------------
## 3. Persist and summarise
## ---------------------------------------------------------------------------
rds_file <- file.path(out_dir, "scan-mt-omnibus-validation.rds")
saveRDS(list(config = list(n_rep = n_rep, n_ind = n_ind, m = m_tr, p = p_snp,
                           alpha = alpha),
             null = null_out, power = power_tab, qq = qq), rds_file)
cat(sprintf("\nAll replicates (failures included) saved to %s\n", rds_file))

cat("\n")
if (length(failures) > 0L) {
  cat("VALIDATION FAILURES:\n")
  for (f in failures) cat(" -", f, "\n")
  quit(status = 1L)
} else {
  cat("All scan_mt_omnibus validation checks passed.\n")
}
