## Stage 8.3 supplementary probe: is strict ++++ feasible for BOTH R1 and
## R3 under Sigma_P = I4 (correlation = "independent")? 200 draws each.
## Result decides whether the II-B-I4 strict-sign supplementary runs.
devtools::load_all(quiet = TRUE)

canon <- function(beta) {
  s <- sign(beta)
  nz <- which(s != 0)
  if (length(nz) == 0L) return(NA_character_)
  if (s[nz[1]] < 0) s <- -s
  paste(ifelse(s > 0, "+", ifelse(s < 0, "-", "0")), collapse = "")
}

one_cell <- function(scenario, target_loss = NULL, n_draws = 200L) {
  ok <- 0L
  concordant <- 0L
  for (i in seq_len(n_draws)) {
    args <- list(n = 200L, m = 4L, p = 300L, experiment = "trait_representation",
                 scenario = scenario, locus_pve = 0.02,
                 correlation = "independent", tolerance = 0.10,
                 seed = 400000 + i)
    if (!is.null(target_loss)) args$target_loss <- target_loss
    s <- do.call(simulate_condped_data, args)
    if (!isTRUE(s$status$ok)) next
    ok <- ok + 1L
    if (identical(canon(s$truth$B_Q[1L, ]), "++++")) concordant <- concordant + 1L
  }
  data.frame(scenario = scenario, correlation = "independent",
             draws = n_draws, accepted = ok, accept_rate = ok / n_draws,
             pp_pp = concordant, pp_rate = concordant / n_draws)
}

cells <- list(
  list("highly_representable", 0.05),
  list("strongly_nonredundant", NULL)
)
out <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(cells, function(c) one_cell(c[[1]], c[[2]]),
                     mc.cores = 2L)
} else {
  lapply(cells, function(c) one_cell(c[[1]], c[[2]]))
}
res <- do.call(rbind, out)
print(res, row.names = FALSE)
write.csv(res, "inst/formal-simu/output/stage83/i4_sign_probe.csv",
          row.names = FALSE)
cat("I4 PROBE DONE\n")
