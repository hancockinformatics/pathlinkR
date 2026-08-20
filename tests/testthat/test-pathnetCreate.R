test_that("the pathway network is the right size", {
    data("sigoraDatabaseHS", "sigoraExamplesHS")

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=dplyr::slice_head(
            dplyr::arrange(sigoraDatabaseHS, pathwayId),
            prop=0.05
        ),
        distMethod="jaccard"
    )

    testStartingPathways <- pathnetFoundation(
        mat=pathwayDistancesJaccard,
        species="human",
        maxDistance=0.8
    )

    testExPathwayNetworkInput <- sigoraExamplesHS %>%
        filter(comparison == "COVID Pos Over Time")

    testMyPathwayNetwork <- pathnetCreate(
        pathwayEnrichmentResult=testExPathwayNetworkInput,
        species="human",
        foundation=testStartingPathways,
        trim=TRUE,
        trimOrder=1
    )

    expect_length(testMyPathwayNetwork, 2)
})
