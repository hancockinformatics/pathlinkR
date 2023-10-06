test_that("the pathway network is the right size", {
    pathwayDistancesJaccard <- getPathwayDistances(
      pathwayData=slice_head(arrange(sigoraDatabase, pathwayId), prop = 0.25),
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

    expect_length(testMyPathwayNetwork, 27)
})
