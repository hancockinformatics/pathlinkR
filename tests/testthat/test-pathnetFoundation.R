test_that("the foundation has the right dimensions and columns", {

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=slice_head(arrange(sigoraDatabase, pathwayId), prop = 0.25),
        distMethod="jaccard"
    )

    testStartingPathways <- pathnetFoundation(
        mat=pathwayDistancesJaccard,
        maxDistance=0.8
    )

    expect_equal(dim(testStartingPathways), c(770, 5))

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
