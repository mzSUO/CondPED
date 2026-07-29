# Contract tests for the CondPED v0.3 S0 API skeleton.
# Verifies, against docs/design/CondPED_v0.3_Interface_Contract.md:
#   1. all 13 public functions are exported;
#   2. every function's formals (names, order, defaults) match the contract;
#   3. every stub throws an error of class condped_not_implemented;
#   4. .safe_inverse() returns the contracted fields on full-rank and
#      rank-deficient matrices.

# Expected signatures, copied verbatim from the interface contract (section 4).
api_formals <- list(
  simulate_condped_data = alist(
    n = 1000L,
    m = 4L,
    p = 2000L,
    maf_range = c(0.05, 0.50),
    h2 = rep(0.50, m),
    R_G = NULL,
    R_E = NULL,
    architecture = c(
      "null",
      "single_trait",
      "shared_same",
      "shared_opposite",
      "dense",
      "covariance_aligned",
      "conditional_deviation",
      "projection_induced"
    ),
    locus_pve = 0.01,
    delta = 0,
    target_trait = 1L,
    n_groups = 8L,
    fst = 0.05,
    structured = TRUE,
    n_qtl = 1L,
    max_attempts = 100L,
    seed = NULL,
    return_latent = TRUE
  ),
  fit_mt_null = alist(
    Y = ,
    W = NULL,
    K = ,
    add_intercept = TRUE,
    optimizer = c("BFGS", "nlminb"),
    n_starts = 5L,
    h2_start = 0.50,
    start_values = NULL,
    eig_tol = 1e-8,
    inverse_tol = sqrt(.Machine$double.eps),
    chol_floor = 1e-8,
    control = list(maxit = 500L, reltol = 1e-8),
    compute_gamma = TRUE,
    return_rotation = TRUE,
    verbose = FALSE
  ),
  scan_mt_omnibus = alist(
    null_fit = ,
    G = ,
    marker_ids = colnames(G),
    maf_min = 0.05,
    chunk_size = 2000L,
    rank_tol = sqrt(.Machine$double.eps),
    return_score = FALSE,
    return_effects = FALSE,
    bootstrap_p = FALSE,
    n_boot = 0L,
    seed = NULL
  ),
  estimate_mt_effects = alist(
    null_fit = ,
    G = ,
    loci = ,
    marker_ids = colnames(G),
    rank_tol = sqrt(.Machine$double.eps),
    return_covariance = TRUE
  ),
  attribute_traits = alist(
    omnibus = ,
    effects = ,
    omnibus_method = c("BH", "bonferroni", "none"),
    alpha_omnibus = 0.05,
    attribution_mode = c("bb_fdr", "holm_fwer", "none"),
    q_target = 0.05,
    alpha_total = 0.05,
    alpha_split = c(0.025, 0.025),
    dependence = c("assumed", "empirical"),
    return_all = TRUE
  ),
  derive_conditional_contrasts = alist(
    Sigma_P = ,
    trait_names = colnames(Sigma_P),
    inverse_tol = sqrt(.Machine$double.eps),
    ridge = 0,
    standardize = FALSE
  ),
  decompose_conditional_effects = alist(
    effects = ,
    attribution = ,
    contrasts = ,
    gate = TRUE,
    alpha = 0.05,
    p_adjust = c("holm", "BH", "none"),
    se_method = c("plugin", "bootstrap"),
    bootstrap_result = NULL
  ),
  estimate_locus_pve = alist(
    effects = ,
    genotype_variance = ,
    Sigma_P = ,
    conditional = NULL,
    correction = c("raw", "noise_corrected", "both"),
    truncate_zero = TRUE,
    gate = TRUE,
    attribution = NULL
  ),
  crossfit_effect_pve = alist(
    Y = ,
    W = ,
    G = ,
    fold_id = NULL,
    n_folds = 2L,
    n_repeats = 10L,
    analysis_args = list(),
    gamma_mode = c("global_fixed", "fold_specific"),
    seed = 1L,
    workers = 1L,
    keep_logs = TRUE
  ),
  bootstrap_condped = alist(
    object = ,
    Y = NULL,
    W = NULL,
    G = NULL,
    loci = NULL,
    B = 200L,
    type = c("parametric", "residual"),
    refit_null = TRUE,
    reselect = FALSE,
    interval = c("percentile", "bca"),
    seed = 1L,
    workers = 1L,
    min_valid = 0.80
  ),
  condped = alist(
    Y = ,
    W = NULL,
    G = ,
    K = NULL,
    chromosome = NULL,
    alpha_omnibus = 0.05,
    omnibus_adjust = "BH",
    attribution_mode = "bb_fdr",
    q_target = 0.05,
    conditional_adjust = "holm",
    pve_correction = "both",
    crossfit = FALSE,
    bootstrap = FALSE,
    control = list(),
    seed = 1L
  ),
  evaluate_condped_simulation = alist(
    files = ,
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
  ),
  run_condped_simulation = alist(
    grid = ,
    reps = ,
    out_dir = ,
    master_seed = 20260727L,
    workers = 1L,
    backend = c("sequential", "parallel"),
    resume = TRUE,
    overwrite = FALSE,
    save_data = FALSE,
    save_fit = TRUE,
    fail_policy = c("save_unstable", "stop"),
    progress = interactive()
  )
)

deparse_default <- function(x) paste(deparse(x), collapse = " ")

test_that("all 13 public API functions are exported", {
  expect_length(api_formals, 13L)
  exports <- getNamespaceExports("CondPED")
  missing <- setdiff(names(api_formals), exports)
  expect_identical(missing, character())
})

test_that("every public function's formals match the interface contract", {
  ns <- asNamespace("CondPED")
  for (fn_name in names(api_formals)) {
    fn <- get(fn_name, envir = ns)
    expect_true(is.function(fn), label = fn_name)
    actual <- formals(fn)
    expected <- api_formals[[fn_name]]
    expect_identical(
      names(actual), names(expected),
      label = paste0(fn_name, " argument names/order")
    )
    for (arg in names(expected)) {
      expect_identical(
        deparse_default(actual[[arg]]),
        deparse_default(expected[[arg]]),
        label = paste0(fn_name, " default for `", arg, "`")
      )
    }
  }
})

test_that("every stub throws condped_not_implemented", {
  # Functions that have graduated from the S0 skeleton to a real
  # implementation no longer throw the stub error.
  implemented <- c("simulate_condped_data", "fit_mt_null",
                   "scan_mt_omnibus", "estimate_mt_effects",
                   "attribute_traits")
  ns <- asNamespace("CondPED")
  for (fn_name in setdiff(names(api_formals), implemented)) {
    fn <- get(fn_name, envir = ns)
    expect_error(
      do.call(fn, list()),
      class = "condped_not_implemented",
      label = fn_name
    )
  }
})

test_that(".safe_inverse() returns the contracted fields on a full-rank matrix", {
  A <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  res <- CondPED:::.safe_inverse(A)
  expect_identical(
    names(res),
    c("inverse", "rank", "eigenvalues", "condition_number", "method", "status")
  )
  expect_identical(res$status, "ok")
  expect_identical(res$method, "chol")
  expect_identical(res$rank, 2L)
  expect_length(res$eigenvalues, 2L)
  expect_true(is.finite(res$condition_number))
  expect_equal(A %*% res$inverse, diag(2), tolerance = 1e-10)
})

test_that(".safe_inverse() degrades gracefully on a rank-deficient matrix", {
  B <- matrix(c(1, 2, 2, 4), 2, 2)
  res <- CondPED:::.safe_inverse(B)
  expect_identical(
    names(res),
    c("inverse", "rank", "eigenvalues", "condition_number", "method", "status")
  )
  expect_identical(res$status, "rank_deficient")
  expect_identical(res$method, "eigen_pinv")
  expect_identical(res$rank, 1L)
  expect_identical(res$condition_number, Inf)
  # Moore-Penrose defining property: B %*% B+ %*% B == B
  expect_equal(B %*% res$inverse %*% B, B, tolerance = 1e-8)
})

test_that(".safe_inverse() treats a numerically near-singular PD matrix as rank deficient", {
  A <- diag(c(1, 1e-20))
  out <- CondPED:::.safe_inverse(A, tol = 1e-8)
  expect_identical(out$rank, 1L)
  expect_identical(out$status, "rank_deficient")
  expect_identical(out$method, "eigen_pinv")
})

test_that(".safe_inverse() returns the zero matrix for a zero matrix (symmetric)", {
  Z <- matrix(0, 2, 2)
  out <- CondPED:::.safe_inverse(Z)
  expect_identical(out$rank, 0L)
  expect_equal(out$inverse, matrix(0, 2, 2))
  expect_identical(out$status, "rank_deficient")
  expect_true(is.infinite(out$condition_number))
})

test_that(".safe_inverse() returns the zero matrix for a zero matrix (non-symmetric)", {
  Z <- matrix(0, 2, 2)
  out <- CondPED:::.safe_inverse(Z, symmetric = FALSE)
  expect_identical(out$rank, 0L)
  expect_equal(out$inverse, matrix(0, 2, 2))
  expect_identical(out$status, "rank_deficient")
  expect_true(is.infinite(out$condition_number))
})
