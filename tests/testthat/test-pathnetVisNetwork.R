test_that("pathVisNetwork returns the right plot", {
    set.seed(1)

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=slice_head(arrange(sigoraDatabase, pathwayId), prop = 0.25),
        distMethod="jaccard"
    )

    startingPathways <- pathnetFoundation(
        mat = pathwayDistancesJaccard,
        maxDistance = 0.8
    )

    exPathnet <- pathnetCreate(
        pathwayEnrichmentResult=dplyr::filter(
            sigoraExamples,
            comparison == "COVID Pos Over Time"
        ),
        foundation=startingPathways,
        trim=TRUE,
        trimOrder=1
    )

    vdiffr::expect_doppelganger(
        "pathnetVisNetworkExample",
        pathnetVisNetwork(exPathnet)
    )
})
