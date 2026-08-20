test_that("pathVisNetwork returns the right plot", {
    data("sigoraDatabaseHS", "sigoraExamplesHS")

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=dplyr::slice_head(
            dplyr::arrange(sigoraDatabaseHS, pathwayId),
            prop=0.05
        ),
        distMethod="jaccard"
    )

    startingPathways <- pathnetFoundation(
        mat=pathwayDistancesJaccard,
        species="human",
        maxDistance=0.8
    )

    exPathnet <- pathnetCreate(
        pathwayEnrichmentResult=dplyr::filter(
            sigoraExamplesHS,
            comparison == "COVID Pos Over Time"
        ),
        species="human",
        foundation=startingPathways,
        trim=TRUE,
        trimOrder=1
    )

    expect_no_error(pathnetVisNetwork(exPathnet))
})
