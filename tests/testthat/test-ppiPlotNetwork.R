## Note we don't test with `label=TRUE` because even with a seed, the placement
## of the labels varies enough (though slightly) between runs to constantly fail
## the test
test_that("we get the right plot output", {
    set.seed(1)

    exNetwork <- ppiBuildNetwork(
        deseqResults=exampleDESeqResults[[1]],
        filterInput=TRUE,
        order="zero"
    )

    vdiffr::expect_doppelganger(
        "example-network",
        ppiPlotNetwork(
            exNetwork,
            fillColumn=log2FoldChange,
            fillType="foldChange",
            legend=TRUE,
            label=FALSE
        )
    )
})

test_that("plotting subnetworks works as expected", {
    set.seed(1)

    exNetwork2 <- ppiBuildNetwork(
        deseqResults=exampleDESeqResults[[1]],
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

    vdiffr::expect_doppelganger(
        "example-subnetwork",
        ppiPlotNetwork(
            network=exSubnetwork,
            fillType="oneSided",
            fillColumn=degree,
            label=TRUE,
            labelColumn=hgncSymbol
        )
    )
})
