test_that("conditional projection removes correlation", {
    set.seed(1)


    Y <- matrix(
        rnorm(1000),
        500,
        2
    )


    Y[, 2] <- Y[, 1] * 0.8 + rnorm(500)


    res <- compute_conditional_phenotype(
        Y
    )


    cors <- cor(
        res$Y_cond[, 1],
        Y[, 2]
    )


    expect_lt(
        abs(cors),
        0.1
    )
})
