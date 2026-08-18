test_that("we get the right plot", {
    data("sigoraDatabaseHS", "sigoraExamples")

    pathwayDistancesJaccard <- getPathwayDistances(
        pathwayData=dplyr::slice_head(
            dplyr::arrange(sigoraDatabaseHS, pathwayId),
            prop=0.05
        ),
        distMethod="jaccard"
    )

    startingPathways <- pathnetFoundation(
        mat=pathwayDistancesJaccard,
        species="human",
        maxDistance=0.8
    )

    exPathnet <- pathnetCreate(
        pathwayEnrichmentResult=dplyr::filter(
            sigoraExamples,
            comparison == "COVID Pos Over Time"
        ),
        species="human",
        foundation=startingPathways,
        trim=TRUE,
        trimOrder=1
    )

    expect_no_error(
        pathnetGGraph(
            exPathnet,
            labelProp=0.1,
            nodeLabelSize=4,
            nodeLabelOverlaps=8,
            segColour="red"
        )
    )
})
