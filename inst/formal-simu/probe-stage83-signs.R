## Stage 8.3 sign-pattern feasibility probe (freeze revision 2026-08-27,
## per-rep sign matching). For each candidate sign pattern (canonical
## form: first trait positive), draw 200 attempts for R1
## (highly_representable, rho = 0.05) and R3 (strongly_nonredundant) and
## report the acceptance rate of draws whose truth is conforming with the
## architecture AND whose canonical sign pattern equals the target.
## Acceptance depends only on population-level truth quantities, so the
## probe uses n = 200 / p = 300 for speed (documented).
devtools::load_all(quiet = TRUE)

patterns <- list(
  "++++" = c(1, 1, 1, 1),
  "+++-" = c(1, 1, 1, -1),
  "++-+" = c(1, 1, -1, 1),
  "++--" = c(1, 1, -1, -1),
  "+-++" = c(1, -1, 1, 1),
  "+-+-" = c(1, -1, 1, -1),
  "+--+" = c(1, -1, -1, 1),
  "+---" = c(1, -1, -1, -1)
)
canon <- function(beta) {
  s <- sign(beta)
  nz <- which(s != 0)
  if (length(nz) == 0L) return(NA_character_)
  if (s[nz[1]] < 0) s <- -s
  paste(ifelse(s > 0, "+", ifelse(s < 0, "-", "0")), collapse = "")
}

archs <- list(
  R1 = list(scenario = "highly_representable", target_loss = 0.05),
  R3 = list(scenario = "strongly_nonredundant", target_loss = NULL)
)

one_cell <- function(arch, pat_name, pat, n_draws = 200L) {
  ok <- 0L
  for (i in seq_len(n_draws)) {
    args <- list(n = 200L, m = 4L, p = 300L, experiment = "trait_representation",
                 scenario = arch$scenario, locus_pve = 0.02,
                 correlation = "block", tolerance = 0.10,
                 seed = 300000 + i)
    if (!is.null(arch$target_loss)) args$target_loss <- arch$target_loss
    s <- do.call(simulate_condped_data, args)
    if (!isTRUE(s$status$ok)) next
    if (identical(canon(s$truth$B_Q[1L, ]), pat_name)) ok <- ok + 1L
  }
  data.frame(arch = arch$scenario, pattern = pat_name, draws = n_draws,
             accepted = ok, rate = ok / n_draws)
}

cells <- expand.grid(arch = names(archs), pattern = names(patterns),
                     stringsAsFactors = FALSE)
workers <- 4L
t0 <- Sys.time()
if (.Platform$OS.type == "unix") {
  out <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
    one_cell(archs[[cells$arch[i]]], cells$pattern[i], patterns[[cells$pattern[i]]])
  }, mc.cores = workers)
} else {
  out <- lapply(seq_len(nrow(cells)), function(i) {
    one_cell(archs[[cells$arch[i]]], cells$pattern[i], patterns[[cells$pattern[i]]])
  })
}
res <- do.call(rbind, out)
rownames(res) <- NULL
print(res, row.names = FALSE)
dir.create("inst/formal-simu/output/stage83", showWarnings = FALSE,
           recursive = TRUE)
write.csv(res, "inst/formal-simu/output/stage83/sign_pattern_probe.csv",
          row.names = FALSE)
cat(sprintf("PROBE DONE wall %.0fs\n", difftime(Sys.time(), t0, units = "secs")))
