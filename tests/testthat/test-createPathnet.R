test_that("the pathway network is the right size", {
    testStartingPathways <- createFoundation(
        mat = pathwayDistancesJaccard,
        maxDistance = 0.8
    )

    testExPathwayNetworkInput <- sigoraExamples %>%
        filter(comparison == "COVID Pos Over Time")

    testMyPathwayNetwork <- createPathnet(
        enrichPathwayResult = testExPathwayNetworkInput,
        foundation = testStartingPathways,
        trim = TRUE,
        trimOrder = 1
    )

    expect_length(testMyPathwayNetwork, 98)
})
