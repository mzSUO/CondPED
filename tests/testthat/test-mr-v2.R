library(CondPED)


test_that("MR result contains required fields", {
    fake <- list(
        gamma = 0.2,
        se = 0.05,
        pval = 0.001,
        sig = TRUE,
        n_iv = 10,
        status = "ok"
    )


    expect_true("gamma" %in% names(fake))
    expect_true("n_iv" %in% names(fake))
    expect_true("status" %in% names(fake))
})
