## Stage 7.4 H: failed-replicate audit across all Stage 7.4 outputs.
## Scans every per-replicate RDS under inst/validation/output/stage74/
## and tabulates non-ok statuses with their messages.
base <- "inst/validation/output/stage74"
rds_files <- list.files(base, pattern = "\\.rds$", recursive = TRUE,
                        full.names = TRUE)
rds_files <- rds_files[!grepl("manifest", rds_files)]

audit <- list()
for (f in rds_files) {
  out <- readRDS(f)
  if (!is.list(out) || length(out) == 0L || !is.list(out[[1L]]) ||
      is.null(out[[1L]]$status)) next
  st <- vapply(out, function(x) {
    if (is.null(x$status)) "missing_status" else as.character(x$status)
  }, character(1))
  bad <- st != "ok"
  audit[[f]] <- data.frame(
    file = f,
    n_total = length(st),
    n_ok = sum(!bad),
    n_failed = sum(bad),
    status_breakdown = paste(sprintf("%s=%d", names(table(st[bad])),
                                     as.integer(table(st[bad]))),
                             collapse = "; "),
    stringsAsFactors = FALSE
  )
}
tab <- do.call(rbind, audit)
if (is.null(tab)) {
  cat("no replicates found\n")
} else {
  rownames(tab) <- NULL
  print(tab, row.names = FALSE)
  write.csv(tab, file.path(base, "failed_replicate_audit.csv"),
            row.names = FALSE)
}
cat("STAGE74 FAILED-REPLICATE AUDIT DONE\n")
