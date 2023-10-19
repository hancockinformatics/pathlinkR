test_that("pathVisNetwork returns the right plot", {
    set.seed(1)

    data("sigoraDatabase", "sigoraExamples")

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=dplyr::slice_head(
            dplyr::arrange(sigoraDatabase, pathwayId),
            prop=0.1
        ),
        distMethod="jaccard"
    )

    startingPathways <- pathnetFoundation(
        mat=pathwayDistancesJaccard,
        maxDistance=0.8
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
