## Note we don't test with `label=TRUE` because even with a seed, the placement
## of the labels varies enough (though slightly) between runs to consistently
## fail the test
test_that("legend toggle is working", {
    data("exampleDESeqResultsHS")

    exNetwork <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResultsHS[[1]],
        filterInput=TRUE,
        order="zero"
    )

    expect_no_error(
        ppiPlotNetwork(
            exNetwork,
            fillColumn=LogFoldChange,
            fillType="foldChange",
            legend=FALSE,
            label=FALSE
        )
    )
})


test_that("we get the right plot output", {
    data("exampleDESeqResultsHS")

    exNetwork <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResultsHS[[1]],
        filterInput=TRUE,
        order="zero"
    )

    expect_no_error(
        ppiPlotNetwork(
            exNetwork,
            fillColumn=LogFoldChange,
            fillType="foldChange",
            legend=TRUE,
            label=FALSE
        )
    )
})

test_that("plotting subnetworks works as expected", {
    data("exampleDESeqResultsHS")

    exNetwork2 <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResultsHS[[1]],
        filterInput=TRUE,
        order="zero"
    )

    exPathways <- ppiEnrichNetwork(
        network=exNetwork2,
        analysis="hallmark"
    )

    exSubnetwork <- ppiExtractSubnetwork(
        network=exNetwork2,
        pathwayEnrichmentResult=exPathways,
        pathwayToExtract="INTERFERON ALPHA RESPONSE"
    )

    expect_no_error(
        ppiPlotNetwork(
            network=exSubnetwork,
            fillColumn=degree,
            fillType="oneSided",
            legendTitle="Degree",
            label=FALSE
        )
    )
})
