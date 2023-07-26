test_that("the foundation has the right dimensions and columns", {
    testStartingPathways <- createFoundation(
        mat = pathwayDistancesJaccard,
        maxDistance = 0.8
    )

    expect_equal(dim(testStartingPathways), c(6852, 5))

    expect_setequal(
        colnames(testStartingPathways),
        c(
            "pathwayName1",
            "pathwayName2",
            "distance",
            "pathway1",
            "pathway2"
        )
    )
})
