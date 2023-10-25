test_that("fold change heatmaps are working", {
    set.seed(1)

    data("exampleDESeqResults")

    vdiffr::expect_doppelganger(
        "plotFoldChangeExample",
        plotFoldChange(
            inputList=exampleDESeqResults,
            pathName="Generation of second messenger molecules"
        )
    )
})

test_that("customization options work correctly", {
    set.seed(1)

    data("exampleDESeqResults")

    vdiffr::expect_doppelganger(
        "plotFoldChangeExample2",
        plotFoldChange(
            inputList=exampleDESeqResults,
            pathName="Interferon alpha/beta signaling",
            hideNonsigFC=FALSE,
            invert=TRUE,
            clusterColumns=TRUE
        )
    )
})
