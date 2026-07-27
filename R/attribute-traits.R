#' Trait attribution with hierarchical error control
#'
#' Attributes the omnibus-selected loci to individual traits, with
#' hierarchical multiplicity correction (Benjamini-Bogomolov style FDR or
#' Holm FWER) layered on top of the first-layer omnibus selection.
#'
#' @section Status:
#' S0 skeleton: this function is an API stub and currently throws an error of
#' class `condped_not_implemented`.
#'
#' @param omnibus Omnibus scan result returned by [scan_mt_omnibus()].
#' @param effects Effect estimates returned by [estimate_mt_effects()].
#' @param omnibus_method Multiple-testing adjustment for the first-layer
#'   omnibus p-values: `"BH"`, `"bonferroni"` or `"none"`.
#' @param alpha_omnibus First-layer omnibus significance level.
#' @param attribution_mode Second-layer correction mode: `"bb_fdr"`
#'   (Benjamini-Bogomolov style FDR on the selected family), `"holm_fwer"`
#'   (split error budget with Holm on all subsequent hypotheses) or `"none"`.
#' @param q_target Target FDR level for `"bb_fdr"`.
#' @param alpha_total Total error budget for `"holm_fwer"`.
#' @param alpha_split Numeric vector of length 2 giving the split of
#'   `alpha_total` between the first and subsequent layers.
#' @param dependence Dependence assumption for the BB step: `"assumed"` or
#'   `"empirical"`.
#' @param return_all Logical; return the full per-trait table, including
#'   non-attributed rows.
#'
#' @return A list with components `selected_loci`, `trait_table` (a
#'   data.frame with columns `marker_id`, `trait`, `beta`, `se`, `p_raw`,
#'   `p_adjusted`, `attributed`), `A` (per-locus associated-trait sets),
#'   `mode`, `assumptions`, `status` and `diagnostics`, as specified in the
#'   interface contract.
#' @export
attribute_traits <- function(
  omnibus,
  effects,
  omnibus_method = c("BH", "bonferroni", "none"),
  alpha_omnibus = 0.05,
  attribution_mode = c("bb_fdr", "holm_fwer", "none"),
  q_target = 0.05,
  alpha_total = 0.05,
  alpha_split = c(0.025, 0.025),
  dependence = c("assumed", "empirical"),
  return_all = TRUE
) {
  .not_implemented("attribute_traits")
}
