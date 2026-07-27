#' Decompose marker effects into conditional deviations
#'
#' Computes the conditional effects \eqn{\widehat\eta = c^\top
#' \widehat\beta} for attributed locus-trait pairs and tests the conditional
#' deviation from the projection-implied value.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param effects Effect estimates returned by [estimate_mt_effects()].
#' @param attribution Attribution result returned by [attribute_traits()].
#' @param contrasts Contrast object returned by
#'   [derive_conditional_contrasts()].
#' @param gate Logical; when `TRUE` conditional tests are computed and
#'   reported only for traits inside the attributed set of each locus.
#' @param alpha Significance level for the conditional tests.
#' @param p_adjust Adjustment method for the conditional p-values:
#'   `"holm"`, `"BH"` or `"none"`.
#' @param se_method Standard error method: `"plugin"` (uses
#'   \eqn{c^\top S c}) or `"bootstrap"`.
#' @param bootstrap_result Optional bootstrap result returned by
#'   [bootstrap_condped()], required when `se_method = "bootstrap"`.
#'
#' @return A list with components `conditional_table` (a data.frame with
#'   columns `marker_id`, `trait`, `eta`, `se`, `p_raw`, `p_adjusted`,
#'   `tested`, `deviating`), `D` (per-locus deviation sets), `status` and
#'   `diagnostics`, as specified in the interface contract.
#' @export
decompose_conditional_effects <- function(
  effects,
  attribution,
  contrasts,
  gate = TRUE,
  alpha = 0.05,
  p_adjust = c("holm", "BH", "none"),
  se_method = c("plugin", "bootstrap"),
  bootstrap_result = NULL
) {
  .not_implemented("decompose_conditional_effects")
}
