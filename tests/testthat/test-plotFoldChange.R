exPlotFoldChange <- plotFoldChange(
    deseqExampleList,
    pathName = "Generation of second messenger molecules"
)

set.seed(1)

test_that("fold change heatmaps are working", {
    vdiffr::expect_doppelganger("plotFoldChangeExample", exPlotFoldChange)
})
