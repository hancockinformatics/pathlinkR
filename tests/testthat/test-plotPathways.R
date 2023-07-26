test_that("pathway plots are correct", {
    set.seed(1)

    vdiffr::expect_doppelganger(
        "pathwayPlotExample",
        plotPathways(
            sigoraExamples,
            columns = 2
        )
    )
})
