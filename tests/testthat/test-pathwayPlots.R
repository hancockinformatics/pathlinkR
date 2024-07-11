test_that("pathway plots are correct", {
    set.seed(1)

    data("sigoraExamples")

    vdiffr::expect_doppelganger(
        "pathwayPlotExample",
        pathwayPlots(sigoraExamples, columns=2)
    )
})

test_that("pathway plots work with Sigora/KEGG", {
    set.seed(1)

    data("exampleDESeqResults")

    suppressMessages(
        testResultSigoraKEGG <- pathwayEnrichment(
            inputList=exampleDESeqResults[1],
            analysis="sigora",
            gpsRepo="kegH"
        )
    )

    vdiffr::expect_doppelganger(
        "pathwayPlotExampleKEGG",
        pathwayPlots(testResultSigoraKEGG, columns=2)
    )
})


test_that("pathwayPlots works with fgsea results", {
    set.seed(1)

    data("exampleDESeqResults")

    suppressMessages(
        testResultFgseaReactome <- pathwayEnrichment(
            inputList=exampleDESeqResults[1],
            analysis="fgsea_reactome"
        )
    )

    vdiffr::expect_doppelganger(
        "pathwayPlotExampleFgseaReactome",
        pathwayPlots(testResultFgseaReactome, columns=3)
    )
})
