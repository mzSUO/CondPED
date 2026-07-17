library(CondPED)


test_that("four traits generate all directional pairs", {
    traits <- c(
        "T1",
        "T2",
        "T3",
        "T4"
    )


    pairs <- combn(
        traits,
        2
    )


    expect_equal(
        ncol(pairs),
        6
    )
})
