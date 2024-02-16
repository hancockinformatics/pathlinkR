test_that("fold change heatmaps are working with some customization", {
    set.seed(1)

    data("exampleDESeqResults")

    vdiffr::expect_doppelganger(
        "plotFoldChangeExample",
        plotFoldChange(
            inputList=exampleDESeqResults,
            pathName="Interferon alpha/beta signaling",
            hideNonsigFC=FALSE,
            invert=TRUE,
            clusterColumns=TRUE
        )
    )
})
