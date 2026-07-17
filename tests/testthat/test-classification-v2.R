library(CondPED)


# ------------------------------------------------------------------------------
# 测试辅助函数：构造符合 bidirectional_mr() 接口的单方向 MR 结果
# ------------------------------------------------------------------------------

make_mr_direction <- function(
    sig = FALSE,
    pval = 0.50,
    n_iv = 10L,
    status = "ok") {
    list(
        gamma = if (sig) 0.30 else 0.02,
        se = 0.05,
        z = if (sig) 6 else 0.40,
        pval = pval,
        sig = sig,
        n_iv = as.integer(n_iv),
        status = status,
        method = "test",
        ivs = paste0("SNP", seq_len(n_iv))
    )
}


make_bidirectional_mr <- function(
    ab_sig = FALSE,
    ba_sig = FALSE,
    ab_pval = if (ab_sig) 0.001 else 0.50,
    ba_pval = if (ba_sig) 0.001 else 0.50,
    ab_n_iv = 10L,
    ba_n_iv = 10L,
    ab_status = "ok",
    ba_status = "ok") {
    list(
        AB = make_mr_direction(
            sig = ab_sig,
            pval = ab_pval,
            n_iv = ab_n_iv,
            status = ab_status
        ),
        BA = make_mr_direction(
            sig = ba_sig,
            pval = ba_pval,
            n_iv = ba_n_iv,
            status = ba_status
        )
    )
}


# ------------------------------------------------------------------------------
# Pattern 1
# ------------------------------------------------------------------------------

test_that("pattern1: only one trait is marginally significant", {
    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = FALSE,
        cond_sig_trait1 = TRUE,
        cond_sig_trait2 = FALSE,
        traits = c("Trait1", "Trait2"),
        return_details = TRUE
    )

    expect_s3_class(res, "data.frame")
    expect_identical(res$pattern, "pattern1")
    expect_identical(res$direction, "none")
    expect_identical(res$confidence, "descriptive")
})


# ------------------------------------------------------------------------------
# No MR provided
# ------------------------------------------------------------------------------

test_that("multi-trait candidate is unresolved when MR was not provided", {
    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = TRUE,
        cond_sig_trait1 = TRUE,
        cond_sig_trait2 = TRUE,
        mr_results = NULL,
        traits = c("Trait1", "Trait2"),
        return_details = TRUE
    )

    expect_identical(res$pattern, "unresolved")
    expect_identical(res$unresolved_reason, "mr_not_provided")
})


# ------------------------------------------------------------------------------
# Pattern 2
# ------------------------------------------------------------------------------

test_that("pattern2: both conditional associations retained and both MR directions non-significant", {
    mr_results <- make_bidirectional_mr(
        ab_sig = FALSE,
        ba_sig = FALSE
    )

    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = TRUE,
        cond_sig_trait1 = TRUE,
        cond_sig_trait2 = TRUE,
        mr_results = mr_results,
        traits = c("Trait1", "Trait2"),
        return_details = TRUE
    )

    expect_identical(res$pattern, "pattern2")
    expect_identical(res$direction, "none")
    expect_identical(res$confidence, "moderate")
    expect_true(is.na(res$unresolved_reason))
})


# ------------------------------------------------------------------------------
# Pattern 3
# ------------------------------------------------------------------------------

test_that("pattern3: forward MR supported and outcome conditional association attenuated", {
    mr_results <- make_bidirectional_mr(
        ab_sig = TRUE,
        ba_sig = FALSE
    )

    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = TRUE,
        cond_sig_trait1 = TRUE,
        cond_sig_trait2 = FALSE,
        mr_results = mr_results,
        traits = c("Trait1", "Trait2"),

        # 本单元测试只验证核心分类逻辑；
        # Steiger、Q 与稳健估计另设诊断测试。
        require_steiger = FALSE,
        reject_heterogeneity = FALSE,
        require_robust_direction = FALSE,
        return_details = TRUE
    )

    expect_identical(res$pattern, "pattern3")
    expect_identical(res$direction, "Trait1->Trait2")
    expect_identical(res$confidence, "direction_supported")
})


test_that("pattern3 also works for the reverse supported direction", {
    mr_results <- make_bidirectional_mr(
        ab_sig = FALSE,
        ba_sig = TRUE
    )

    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = TRUE,

        # BA 为支持方向，因此 Trait1 是结局；
        # Trait1 条件关联未保留，对应 pattern3。
        cond_sig_trait1 = FALSE,
        cond_sig_trait2 = TRUE,
        mr_results = mr_results,
        traits = c("Trait1", "Trait2"),
        require_steiger = FALSE,
        reject_heterogeneity = FALSE,
        require_robust_direction = FALSE,
        return_details = TRUE
    )

    expect_identical(res$pattern, "pattern3")
    expect_identical(res$direction, "Trait2->Trait1")
})


# ------------------------------------------------------------------------------
# Pattern 4
# ------------------------------------------------------------------------------

test_that("pattern4: forward MR supported and outcome conditional association retained", {
    mr_results <- make_bidirectional_mr(
        ab_sig = TRUE,
        ba_sig = FALSE
    )

    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = TRUE,
        cond_sig_trait1 = TRUE,
        cond_sig_trait2 = TRUE,
        mr_results = mr_results,
        traits = c("Trait1", "Trait2"),
        require_steiger = FALSE,
        reject_heterogeneity = FALSE,
        require_robust_direction = FALSE,
        return_details = TRUE
    )

    expect_identical(res$pattern, "pattern4")
    expect_identical(res$direction, "Trait1->Trait2")
    expect_identical(res$confidence, "direction_supported")
})


# ------------------------------------------------------------------------------
# Pattern 5
# ------------------------------------------------------------------------------

test_that("pattern5: both MR directions are significant", {
    mr_results <- make_bidirectional_mr(
        ab_sig = TRUE,
        ba_sig = TRUE
    )

    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = TRUE,
        cond_sig_trait1 = TRUE,
        cond_sig_trait2 = TRUE,
        mr_results = mr_results,
        traits = c("Trait1", "Trait2"),
        return_details = TRUE
    )

    expect_identical(res$pattern, "pattern5")
    expect_identical(res$direction, "bidirectional_or_ambiguous")
    expect_identical(res$confidence, "ambiguous")
})


# ------------------------------------------------------------------------------
# Unresolved
# ------------------------------------------------------------------------------

test_that("unresolved when one or both MR directions have insufficient IVs", {
    mr_results <- make_bidirectional_mr(
        ab_sig = TRUE,
        ba_sig = FALSE,
        ab_n_iv = 1L,
        ba_n_iv = 1L,
        ab_status = "insufficient_iv",
        ba_status = "insufficient_iv"
    )

    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = TRUE,
        cond_sig_trait1 = TRUE,
        cond_sig_trait2 = FALSE,
        mr_results = mr_results,
        traits = c("Trait1", "Trait2"),
        return_details = TRUE
    )

    expect_identical(res$pattern, "unresolved")
    expect_match(
        res$unresolved_reason,
        "insufficient_iv"
    )
})


test_that("unresolved when conditional attenuation exists but neither direction is supported", {
    mr_results <- make_bidirectional_mr(
        ab_sig = FALSE,
        ba_sig = FALSE
    )

    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = TRUE,
        cond_sig_trait1 = TRUE,
        cond_sig_trait2 = FALSE,
        mr_results = mr_results,
        traits = c("Trait1", "Trait2"),
        return_details = TRUE
    )

    expect_identical(res$pattern, "unresolved")
    expect_identical(
        res$unresolved_reason,
        "no_directional_support_with_attenuated_conditional_pattern"
    )
})


# ------------------------------------------------------------------------------
# Not detected
# ------------------------------------------------------------------------------

test_that("not_detected when neither trait is marginally significant", {
    res <- classify_pair_pattern(
        marg_sig_trait1 = FALSE,
        marg_sig_trait2 = FALSE,
        cond_sig_trait1 = FALSE,
        cond_sig_trait2 = FALSE,
        traits = c("Trait1", "Trait2"),
        return_details = TRUE
    )

    expect_identical(res$pattern, "not_detected")
    expect_identical(res$confidence, "not_applicable")
})


# ------------------------------------------------------------------------------
# Scalar return
# ------------------------------------------------------------------------------

test_that("return_details FALSE returns one pattern string", {
    res <- classify_pair_pattern(
        marg_sig_trait1 = TRUE,
        marg_sig_trait2 = FALSE,
        cond_sig_trait1 = TRUE,
        cond_sig_trait2 = FALSE,
        return_details = FALSE
    )

    expect_identical(res, "pattern1")
})
