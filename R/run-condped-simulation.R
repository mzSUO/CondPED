#' Run CondPED simulation scenarios
#'
#' Executes a grid of simulation scenarios with reproducible random seeds.
#' Each replicate is uniquely determined by `(master_seed, scenario_id,
#' rep_id)`; parallel scheduling never changes the random stream. Failed
#' replicates are saved with their status instead of being silently
#' dropped.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param grid Data.frame describing the scenario grid (one row per
#'   scenario).
#' @param reps Integer; number of replicates per scenario.
#' @param out_dir Output directory for per-replicate RDS files.
#' @param master_seed Master random seed; replicates use
#'   `RNGkind("L'Ecuyer-CMRG")` streams derived from it.
#' @param workers Integer; number of parallel workers.
#' @param backend Execution backend: `"sequential"` or `"parallel"`.
#' @param resume Logical; skip replicates whose output files already exist.
#' @param overwrite Logical; overwrite existing output files.
#' @param save_data Logical; also save the simulated data.
#' @param save_fit Logical; save the fitted objects.
#' @param fail_policy `"save_unstable"` saves failed replicates with status
#'   information; `"stop"` aborts on the first failure.
#' @param progress Logical; show a progress indicator.
#'
#' @return A list with components `manifest` (data.frame), `completed`,
#'   `failed`, `skipped`, `elapsed` and `status`, as specified in the
#'   interface contract.
#' @export
run_condped_simulation <- function(
  grid,
  reps,
  out_dir,
  master_seed = 20260727L,
  workers = 1L,
  backend = c("sequential", "parallel"),
  resume = TRUE,
  overwrite = FALSE,
  save_data = FALSE,
  save_fit = TRUE,
  fail_policy = c("save_unstable", "stop"),
  progress = interactive()
) {
  .not_implemented("run_condped_simulation")
}
