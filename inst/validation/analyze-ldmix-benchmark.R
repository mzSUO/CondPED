## Stage 7.3 B analysis: summarise the controlled LD-mixing benchmark.
res <- readRDS("inst/validation/output/ldmix-benchmark/benchmark_rows.rds")
res <- res[res$status == "ok", ]

summ <- do.call(rbind, lapply(split(res, res$psid), function(x) {
  data.frame(
    psid = x$psid[1],
    n = nrow(x),
    empirical_r2_median = round(stats::median(x$empirical_r2_causal,
                                              na.rm = TRUE), 3),
    P_lead_has_T1_T2 = round(mean(x$lead_has_both_causal_traits), 3),
    lead_breadth_error = round(mean(x$lead_breadth_error), 3),
    P_signal_count_recovered = round(mean(x$signal_count_recovered), 3),
    res_exact_trait_rate = round(mean(x$res_exact_trait_rate,
                                      na.rm = TRUE), 3),
    res_breadth_error = round(mean(x$res_mean_breadth_error,
                                   na.rm = TRUE), 3),
    P_paired_improvement = round(mean(x$paired_improvement > 0,
                                      na.rm = TRUE), 3),
    median_rep_to_causal_r2 = round(stats::median(
      x$median_rep_to_causal_r2, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
}))
rownames(summ) <- NULL
# order by r2 then scale then ratio
summ <- summ[order(summ$psid), ]
options(width = 200)
print(summ, row.names = FALSE)
write.csv(summ, "inst/validation/output/ldmix-benchmark/benchmark_summary.csv",
          row.names = FALSE)
cat("\nBENCHMARK SUMMARY DONE\n")
