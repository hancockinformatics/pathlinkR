test_that("pathway plots are correct", {
    data("sigoraExamplesHS")
    expect_no_error(pathwayPlots(sigoraExamplesHS, columns=2))
})

test_that("pathway plots work with Sigora/KEGG", {
    data("exampleDESeqResultsHS")

    suppressMessages(
        testResultSigoraKEGG <- pathwayEnrichment(
            inputList=exampleDESeqResultsHS[1],
            species="human",
            analysis="sigora",
            gpsRepo="kegH",
            gpsLevel=2
        )
    )

    expect_no_error(pathwayPlots(testResultSigoraKEGG, columns=2))
})


test_that("pathwayPlots works with fgsea results", {
    data("exampleDESeqResultsHS")

    suppressMessages(
        testResultFgseaReactome <- pathwayEnrichment(
            inputList=exampleDESeqResultsHS[1],
            species="human",
            analysis="fgsea_reactome"
        )
    )

    expect_no_error(pathwayPlots(testResultFgseaReactome, columns=3))
})
