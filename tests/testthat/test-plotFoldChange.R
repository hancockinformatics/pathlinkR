test_that("fold change heatmaps are working", {
    set.seed(1)

    vdiffr::expect_doppelganger(
        "plotFoldChangeExample",
        plotFoldChange(
            exampleDESeqResults,
            pathName="Generation of second messenger molecules"
        )
    )
})

test_that("customization options work correctly", {
    set.seed(1)

    vdiffr::expect_doppelganger(
        "plotFoldChangeExample2",
        plotFoldChange(
            exampleDESeqResults,
            pathName="Interferon alpha/beta signaling",
            hideNonsigFC=FALSE,
            invert=TRUE,
            clusterColumns=TRUE
        )
    )
})
