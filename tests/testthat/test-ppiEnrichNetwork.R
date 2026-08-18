test_that("network enrichment works properly", {
    data("exampleDESeqResults")

    exNetwork <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResults[[1]],
        species="human",
        filterInput=TRUE,
        order="zero"
    )

    exOutput <- ppiEnrichNetwork(
        network=exNetwork,
        species="human",
        analysis="sigora",
        gpsRepo="reaH"
    )

    expect_equal(dim(exOutput), c(20, 10))
})
