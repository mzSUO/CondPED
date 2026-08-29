## Stage 8.3 II-B rep_exact diagnostic (read-only; no code/result changes).
## Decomposes R1 rep_exact failures into threshold misclassification
## (singleton rho_hat >= tau) vs cardinality/other errors; reports the
## singleton rho_hat distribution vs tau = 0.10; same breakdown for R3.
devtools::load_all(quiet = TRUE)

out_dir <- "inst/formal-simu/output/stage83"
tau <- 0.10

rep_diag <- function(f) {
  o <- readRDS(f)
  if (is.null(o) || !isTRUE(o$status$ok)) return(NULL)
  arch <- if (!is.null(o$settings$architecture) &&
              nzchar(o$settings$architecture)) {
    o$settings$architecture
  } else {
    o$settings$scenario
  }
  truth_mrs <- o$truth$minimum_representative_sets
  ## truth-side reference: the truth min-Rep set (R1: singleton; R3: full)
  truth_key <- truth_mrs$trait_key[which.min(truth_mrs$set_size)[1]]
  tmap <- as.data.frame(o$truth$representation_map)
  truth_rho <- tmap$representation_loss[tmap$representing_key == truth_key][1]

  st <- o$subset_table
  est_rho <- st$representation_loss[st$representing_key == truth_key][1]
  tp <- o$tolerance_path[o$tolerance_path$tolerance == tau, ]
  est_min_size <- tp$minimum_set_size
  est_keys <- unique(o$minimum_representative_sets$trait_key)
  rep_exact <- o$evaluation$rep_exact_family_recovery
  rep_card <- o$evaluation$rep_min_cardinality_recovery

  cls <- "ok"
  if (isTRUE(rep_exact < 1)) {
    if (is.na(est_rho)) {
      cls <- "singleton_never_evaluated"
    } else if (est_rho >= tau && truth_rho < tau) {
      cls <- "threshold_cross"          # rho_hat crossed tau upwards
    } else if (est_min_size == 1L && !(truth_key %in% est_keys)) {
      cls <- "wrong_singleton_trait"
    } else if (est_min_size != 1L && arch == "highly_representable") {
      cls <- "cardinality_error_other"
    } else {
      cls <- "family_mismatch_other"
    }
  }
  data.frame(arch = arch, rep_id = o$rep_id, truth_key = truth_key,
             truth_rho = truth_rho, est_rho = est_rho,
             est_min_size = est_min_size,
             n_est_rep_sets = length(est_keys),
             rep_exact = rep_exact, rep_card = rep_card, cls = cls,
             stringsAsFactors = FALSE)
}

dirs <- c("3324782832", "1172990067")
files <- unlist(lapply(dirs, function(d) {
  list.files(file.path(out_dir, d), pattern = "rep_.*rds", full.names = TRUE)
}))
workers <- 4L
rows <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(files, rep_diag, mc.cores = workers)
} else {
  lapply(files, rep_diag)
}
rows <- rows[!vapply(rows, is.null, logical(1))]
d <- do.call(rbind, rows)
write.csv(d, file.path(out_dir, "repflip_diagnostics.csv"), row.names = FALSE)

r1 <- d[d$arch == "highly_representable", ]
r3 <- d[d$arch == "strongly_nonredundant", ]

cat("== R1 failure classes ==\n")
print(table(r1$cls))
cat(sprintf("R1 singleton rho_hat: median %.4f, q90 %.4f, q99 %.4f, P(>=tau) %.3f\n",
            stats::median(r1$est_rho), as.numeric(quantile(r1$est_rho, 0.9)),
            as.numeric(quantile(r1$est_rho, 0.99)), mean(r1$est_rho >= tau)))
cat(sprintf("R1 failures with est_rho >= tau: %d / %d\n",
            sum(r1$cls != "ok" & r1$est_rho >= tau), sum(r1$cls != "ok")))
cat(sprintf("R1 threshold_cross share of failures: %.3f\n",
            mean(r1$cls[r1$cls != "ok"] == "threshold_cross")))
cat("== R3 failure classes ==\n")
print(table(r3$cls))

png(file.path(out_dir, "r1_singleton_rho_hat.png"), width = 900, height = 560)
hist(r1$est_rho, breaks = 60, col = "grey80", border = "grey40",
     main = "R1 singleton rho_hat distribution (n = 500)",
     xlab = expression(hat(rho)[{Trait1}]))
abline(v = tau, col = "red", lwd = 2, lty = 2)
abline(v = stats::median(r1$est_rho), col = "blue", lwd = 2, lty = 3)
legend("topright",
       legend = c(sprintf("tau = %.2f", tau),
                  sprintf("median = %.3f", stats::median(r1$est_rho)),
                  sprintf("P(rho_hat >= tau) = %.3f", mean(r1$est_rho >= tau))),
       col = c("red", "blue", NA), lty = c(2, 3, NA), lwd = c(2, 2, NA),
       bty = "n")
dev.off()

## truth-side singleton rho (fixed at target 0.05 for R1) for reference
cat(sprintf("R1 truth singleton rho: unique values %s\n",
            paste(round(unique(r1$truth_rho), 4), collapse = ", ")))
cat("REFFLIP DIAGNOSTIC DONE\n")
