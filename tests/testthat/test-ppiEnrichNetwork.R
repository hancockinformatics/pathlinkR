test_that("network enrichment works properly", {

    exNetwork <- ppiBuildNetwork(
        deseqResults=exampleDESeqResults[[1]],
        filterInput=TRUE,
        order="zero"
    )

    exOutput <- ppiEnrichNetwork(
        network=exNetwork,
        analysis="sigora",
        gpsRepo="default"
    )

    expect_equal(dim(exOutput), c(20, 10))
})
