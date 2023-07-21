set.seed(1)
exPathwayPlot <- plotPathways(
    sigoraExamples,
    columns = 2
)

test_that("pathway plots are correct", {
    vdiffr::expect_doppelganger("pathwayPlotExample", exPathwayPlot)
})
