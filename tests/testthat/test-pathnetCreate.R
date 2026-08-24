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
    expect_equal(nrow(as_tibble(testMyPathwayNetwork)), 2)
})

test_that("the pathway network is the right size in mice", {
    data("sigoraDatabaseMM", "sigoraExamplesMM")

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=dplyr::slice_head(
            dplyr::arrange(sigoraDatabaseMM, pathwayId),
            prop=0.1
        ),
        distMethod="jaccard"
    )

    testStartingPathways <- pathnetFoundation(
        mat=pathwayDistancesJaccard,
        species="mouse",
        maxDistance=0.9
    )

    testMyPathwayNetwork <- pathnetCreate(
        pathwayEnrichmentResult=sigoraExamplesMM,
        species="mouse",
        foundation=testStartingPathways,
        trim=TRUE,
        trimOrder=1
    )

    expect_length(testMyPathwayNetwork, 4)
    expect_equal(nrow(as_tibble(testMyPathwayNetwork)), 4)
})
