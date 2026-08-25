## Rcpp migration benchmark: scan_mt_omnibus() wall time, n=800, p=100k.
## Compares the C++ chunk kernel against the R reference path (per-marker
## .gls_block_components_r) on the same fitted null model and genotype
## matrix (with a small number of NA dosages).
devtools::load_all(quiet = TRUE)
options(digits = 10)

set.seed(20260826)
n <- 800L; m <- 6L; p <- 100000L
Y <- matrix(rnorm(n * m), n, m, dimnames = list(NULL, paste0("T", 1:m)))
K <- CondPED:::.make_grm(matrix(rnorm(n * 100), n, 100))
fit <- fit_mt_null(Y, K = K, control = list(maxit = 200))
stopifnot(isTRUE(fit$status$ok))

G <- matrix(rbinom(n * p, 2, 0.3), n, p)
colnames(G) <- paste0("M", seq_len(p))
na_idx <- sample(length(G), 200)
G[na_idx] <- NA_integer_

out_file <- "inst/validation/output/rcpp-benchmark.rds"

## ---- C++ path (current scan_mt_omnibus) -----------------------------------
t_cpp <- system.time(scan_cpp <- scan_mt_omnibus(fit, G))[["elapsed"]]

## ---- R reference path: same scan flow, per-marker R blocks -----------------
## Replicates scan_mt_omnibus exactly, swapping only the block kernel.
scan_r_reference <- function(null_fit, G, chunk_size = 2000L,
                             rank_tol = sqrt(.Machine$double.eps)) {
  rot <- null_fit$rotation
  n <- nrow(rot$Y_tilde); m <- ncol(rot$Y_tilde)
  A_arr <- rot$Vinv
  AM_arr <- CondPED:::.precompute_AM(rot)
  Ar <- CondPED:::.precompute_Ar(rot)
  G_inv <- rot$XtVinvX_inv
  p <- ncol(G)
  Qv <- numeric(p); dfv <- integer(p); pv <- numeric(p)
  chunk_starts <- seq(1L, p, by = chunk_size)
  for (cs in chunk_starts) {
    idx <- cs:min(p, cs + chunk_size - 1L)
    G_chunk <- G[, idx, drop = FALSE]
    G_imp <- G_chunk
    if (anyNA(G_imp)) {
      cmu <- colMeans(G_imp, na.rm = TRUE)
      for (s in seq_len(ncol(G_imp))) {
        if (anyNA(G_imp[, s]) && is.finite(cmu[s])) {
          G_imp[is.na(G_imp[, s]), s] <- cmu[s]
        }
      }
    }
    X_tilde <- crossprod(rot$U, G_imp)
    for (s in seq_along(idx)) {
      block <- CondPED:::.gls_block_components_r(X_tilde[, s], AM_arr, Ar,
                                                 A_arr, G_inv, rank_tol)
      qs <- CondPED:::.q_from_block(block, m)
      Qv[idx[s]] <- qs$Q
      dfv[idx[s]] <- block$rank
      pv[idx[s]] <- qs$p_value
    }
  }
  list(Q = Qv, df = dfv, p = pv)
}

## R reference on a 5k-marker subset (extrapolated); full R run is infeasible
sub <- seq_len(5000L)
t_r_sub <- system.time(ref <- scan_r_reference(fit, G[, sub]))[["elapsed"]]
t_r_est <- t_r_sub / length(sub) * p

om <- scan_cpp$omnibus
max_dQ <- max(abs(om$Q[sub] - ref$Q), na.rm = TRUE)
max_dp <- max(abs(om$p_value[sub] - ref$p), na.rm = TRUE)
same_rank <- identical(om$rank_J[sub], as.integer(ref$df))

res <- list(
  n = n, m = m, p = p,
  wall_cpp_seconds = t_cpp,
  wall_r_subset_seconds = t_r_sub,
  wall_r_estimated_seconds = t_r_est,
  speedup_estimated = t_r_est / t_cpp,
  max_abs_dQ_subset = max_dQ,
  max_abs_dp_subset = max_dp,
  rank_identical_subset = same_rank
)
saveRDS(res, out_file)
cat(sprintf("C++ wall: %.1fs | R subset (5k): %.1fs -> est full %.0fs | speedup x%.0f\n",
            t_cpp, t_r_sub, t_r_est, t_r_est / t_cpp))
cat(sprintf("numerics: max|dQ|=%.2e max|dp|=%.2e rank identical: %s\n",
            max_dQ, max_dp, same_rank))
cat("BENCHMARK DONE\n")
