test_that("there are the right number of nodes", {
    testStartingPathways <- createFoundation(
        mat = pathwayDistancesJaccard,
        maxDistance = 0.8
    )

    testExPathwayNetworkInput <- sigoraExamples %>%
        filter(comparison == "COVID Pos Over Time")

    testMyPathwayNetwork <- createPathnet(
        sigoraResult = testExPathwayNetworkInput,
        foundation = testStartingPathways,
        trim = TRUE,
        trimOrder = 1
    )

    expect_length(testMyPathwayNetwork, 98)
})
