test_that("we get the right plot", {
    set.seed(1)

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
        "pathnetGGraphExample",
        pathnetGGraph(
            exPathnet,
            labelProp=0.1,
            nodeLabelSize=4,
            nodeLabelOverlaps=8,
            segColour="red"
        )
    )
})
