testStartingPathways <- createFoundation(
    mat = pathwayDistancesJaccard,
    maxDistance = 0.8
)

testExPathwayNetworkInput <- sigoraExamples %>%
    filter(comparison == "COVID Pos Over Time")

testMyPathwayNetwork <- createPathnet(
    sigora_result = testExPathwayNetworkInput,
    foundation = testStartingPathways,
    trim = TRUE,
    trimOrder = 1
)

test_that("there are the right number of nodes", {
    expect_length(testMyPathwayNetwork, 168)
})
