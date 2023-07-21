exNetwork <- ppiBuildNetwork(
    deseqResults = deseqExampleList[[1]],
    filterInput = TRUE,
    order = "zero"
)

set.seed(1)
exNetworkPlot <- ppiPlotNetwork(
    exNetwork,
    fillColumn = log2FoldChange,
    fillType = "foldChange",
    layout = "lgl",
    legend = TRUE,
    label = FALSE
)

test_that("we get the right plot output", {
    vdiffr::expect_doppelganger("example-network", exNetworkPlot)
})
