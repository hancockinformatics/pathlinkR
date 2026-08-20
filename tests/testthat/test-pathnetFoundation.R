test_that("the foundation has the right dimensions and columns", {
    data("sigoraDatabaseHS")

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=dplyr::slice_head(
            dplyr::arrange(sigoraDatabaseHS, pathwayId),
            prop=0.05
        ),
        distMethod="jaccard"
    )

    testStartingPathways <- pathnetFoundation(
        mat=pathwayDistancesJaccard,
        maxDistance=0.8,
        species="human"
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

test_that("the foundation has the right dimensions and columns in mice", {
    data("sigoraDatabaseMM")

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=dplyr::slice_head(
            dplyr::arrange(sigoraDatabaseMM, pathwayId),
            prop=0.05
        ),
        distMethod="jaccard"
    )

    testStartingPathways <- pathnetFoundation(
        mat=pathwayDistancesJaccard,
        maxDistance=0.8,
        species="mouse"
    )

    expect_equal(dim(testStartingPathways), c(52, 5))

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
