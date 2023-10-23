test_that("fold change heatmaps are working", {
    set.seed(1)

    data("exampleDESeqResults")

    vdiffr::expect_doppelganger(
        "plotFoldChangeExample",
        plotFoldChange(
            inputList=exampleDESeqResults,
            columnFC="log2FoldChange",
            columnP="padj",
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
            columnFC="log2FoldChange",
            columnP="padj",
            pathName="Interferon alpha/beta signaling",
            hideNonsigFC=FALSE,
            invert=TRUE,
            clusterColumns=TRUE
        )
    )
})
