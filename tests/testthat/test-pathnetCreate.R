test_that("the pathway network is the right size", {
    data("sigoraDatabase", "sigoraExamples")

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

    testExPathwayNetworkInput <- sigoraExamples %>%
        filter(comparison == "COVID Pos Over Time")

    testMyPathwayNetwork <- pathnetCreate(
        pathwayEnrichmentResult=testExPathwayNetworkInput,
        foundation=testStartingPathways,
        trim=TRUE,
        trimOrder=1
    )

    expect_length(testMyPathwayNetwork, 2)
})
