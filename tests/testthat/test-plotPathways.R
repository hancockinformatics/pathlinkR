exPathwayPlot <- plotPathways(
    sigoraExamples,
    columns = 2
)

test_that("multiplication works", {
    vdiffr::expect_doppelganger("pathwayPlotExample", exPathwayPlot)
})
