test_that("pathVisNetwork returns the right plot", {
    set.seed(1)

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

    vdiffr::expect_doppelganger(
        "pathnetVisNetworkExample",
        pathnetVisNetwork(exPathnet)
    )
})
