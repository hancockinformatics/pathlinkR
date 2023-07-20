exPathwayPlot <- plotPathways(
    sigoraExamples,
    columns = 2
)

set.seed(1)

test_that("pathway plots are correct", {
    vdiffr::expect_doppelganger("pathwayPlotExample", exPathwayPlot)
})
