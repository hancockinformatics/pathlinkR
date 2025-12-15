test_that("pathway plots are correct", {
    

    expect_no_error(pathwayPlots(sigoraExamples, columns=2))
})

test_that("pathway plots work with Sigora/KEGG", {
    data("exampleDESeqResults")

    suppressMessages(
        testResultSigoraKEGG <- pathwayEnrichment(
            inputList=exampleDESeqResults[1],
            analysis="sigora",
            gpsRepo="kegH"
        )
    )

    expect_no_error(pathwayPlots(testResultSigoraKEGG, columns=2))
})


test_that("pathwayPlots works with fgsea results", {
    data("exampleDESeqResults")

    suppressMessages(
        testResultFgseaReactome <- pathwayEnrichment(
            inputList=exampleDESeqResults[1],
            analysis="fgsea_reactome"
        )
    )

    expect_no_error(pathwayPlots(testResultFgseaReactome, columns=3))
})
