test_that("the foundation has the right dimensions and columns", {
    data("sigoraDatabase")

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=dplyr::slice_head(
            dplyr::arrange(sigoraDatabase, pathwayId),
            prop=0.05
        ),
        distMethod="jaccard"
    )

    testStartingPathways <- pathnetFoundation(
        mat=pathwayDistancesJaccard,
        maxDistance=0.8
    )

    expect_equal(dim(testStartingPathways), c(80, 5))

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
