testStartingPathways <- createFoundation(
    mat = pathwayDistancesJaccard,
    maxDistance = 0.8
)

test_that("the foundation has the right dimensions", {
  expect_equal(dim(testStartingPathways), c(6852, 5))
})

test_that("the foundation has the right columns", {
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
