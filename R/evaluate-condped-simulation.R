#' Summarise CondPED simulation results
#'
#' Aggregates per-replicate simulation outputs into evaluation metrics
#' (type-I error, power, set-recovery metrics, bias/RMSE/coverage of
#' effects, conditional effects and PVE, and runtime), with Monte Carlo
#' standard errors. Failed or unstable replicates are never silently
#' dropped.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param files Character vector of RDS files produced by
#'   [run_condped_simulation()].
#' @param metrics Character vector of metrics to compute; any subset of
#'   `"type1_omnibus"`, `"power_omnibus"`, `"tpr_A"`, `"fdp_A"`, `"esr_A"`,
#'   `"jaccard_A"`, `"tpr_D"`, `"fdp_D"`, `"esr_D"`, `"bias_beta"`,
#'   `"rmse_beta"`, `"coverage_beta"`, `"bias_eta"`, `"rmse_eta"`,
#'   `"coverage_eta"`, `"bias_pve"`, `"rmse_pve"`, `"coverage_pve"`,
#'   `"runtime"`.
#' @param include_unstable Logical; include unstable replicates in the
#'   summary (they are still reported separately, never silently dropped).
#' @param group_by Character vector of grouping variables for the summary.
#' @param conf_level Confidence level for Monte Carlo intervals.
#'
#' @return A list with components `summary` (data.frame),
#'   `replicate_level` (data.frame), `failures` (data.frame),
#'   `monte_carlo_se` (data.frame) and `status`, as specified in the
#'   interface contract.
#' @export
evaluate_condped_simulation <- function(
  files,
  metrics = c(
    "type1_omnibus",
    "power_omnibus",
    "tpr_A",
    "fdp_A",
    "esr_A",
    "jaccard_A",
    "tpr_D",
    "fdp_D",
    "esr_D",
    "bias_beta",
    "rmse_beta",
    "coverage_beta",
    "bias_eta",
    "rmse_eta",
    "coverage_eta",
    "bias_pve",
    "rmse_pve",
    "coverage_pve",
    "runtime"
  ),
  include_unstable = TRUE,
  group_by = c("experiment", "n", "architecture", "pve", "correlation"),
  conf_level = 0.95
) {
  .not_implemented("evaluate_condped_simulation")
}
