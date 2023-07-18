startingPathways <- createFoundation(
    mat = pathwayDistancesJaccard,
    maxDistance = 0.8
)

exPathnet <- createPathnet(
    sigoraResult = dplyr::filter(
        sigoraExamples,
        comparison == "COVID Pos Over Time"
    ),
    foundation = startingPathways,
    trim = TRUE,
    trimOrder = 1
)

set.seed(1)

exPathnetVisNetwork <- pathnetVisNetwork(exPathnet)

test_that("multiplication works", {
    vdiffr::expect_doppelganger(
        "pathnetVisNetworkExample",
        exPathnetVisNetwork
    )
})
