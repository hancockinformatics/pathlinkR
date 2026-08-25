test_that("network enrichment works properly", {
    data("exampleDESeqResultsHS")
    data("reaH", package="sigora")

    exNetwork <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResultsHS[[1]],
        species="human",
        filterInput=TRUE,
        order="zero"
    )

    exOutput <- ppiEnrichNetwork(
        network=exNetwork,
        species="human",
        analysis="sigora",
        gpsRepo=reaH,
        gpsLevel=4
    )

    expect_equal(dim(exOutput), c(20, 10))
})

test_that("network enrichment works properly", {
    data("exampleDESeqResultsMM", "innateDbPPIMM")
    data("reaM", package="sigora")

    exNetwork <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResultsMM[[1]],
        species="mouse",
        filterInput=TRUE,
        order="zero",
        ppiData = innateDbPPIMM
    )

    exOutput <- ppiEnrichNetwork(
        network=exNetwork,
        species="mouse",
        analysis="sigora",
        gpsRepo=reaM,
        gpsLevel=4
    )

    expect_equal(dim(exOutput), c(4, 10))
})
