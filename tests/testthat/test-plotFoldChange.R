exPlotFoldChange <- plotFoldChange(
    deseqExampleList,
    pathName = "Generation of second messenger molecules"
)

test_that("fold change heatmaps are working", {
    vdiffr::expect_doppelganger("plotFoldChangeExample", exPlotFoldChange)
})
