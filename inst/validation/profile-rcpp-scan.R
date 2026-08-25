## profvis hotspot check after the Rcpp migration (p = 5k subset suffices).
## Expect: scan time dominated by crossprod / gls_blocks_cpp, with no R-level
## per-marker loop hotspot.
devtools::load_all(quiet = TRUE)

set.seed(20260826)
n <- 800L; m <- 6L; p <- 5000L
Y <- matrix(rnorm(n * m), n, m, dimnames = list(NULL, paste0("T", 1:m)))
K <- CondPED:::.make_grm(matrix(rnorm(n * 100), n, 100))
fit <- fit_mt_null(Y, K = K, control = list(maxit = 200))
G <- matrix(rbinom(n * p, 2, 0.3), n, p)
colnames(G) <- paste0("M", seq_len(p))

prof <- profvis::profvis(scan_mt_omnibus(fit, G), interval = 0.01)
out <- "inst/validation/output/rcpp-profvis.html"
tryCatch(
  htmlwidgets::saveWidget(prof, out, selfcontained = FALSE),
  error = function(e) cat("html save skipped:", conditionMessage(e), "\n")
)

## text summary: self-time by function
tab <- prof$x$message$prof
if (!is.null(tab)) {
  agg <- aggregate(time ~ label, data = tab, sum)
  agg <- agg[order(-agg$time), ]
  cat("top self-time labels:\n")
  print(head(agg, 10), row.names = FALSE)
}
cat("PROFVIS DONE\n")
