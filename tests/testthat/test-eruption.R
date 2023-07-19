set.seed(1)

exVolcanoPlot <- eruption(deseqResults = deseqExampleList[[1]])

test_that("the volcano plot is correct", {
  vdiffr::expect_doppelganger("volcanoPlotExample", exVolcanoPlot)
})
