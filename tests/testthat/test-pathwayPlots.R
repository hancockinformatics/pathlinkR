test_that("pathway plots are correct", {
    set.seed(1)

    vdiffr::expect_doppelganger(
        "pathwayPlotExample",
        pathwayPlots(
            sigoraExamples,
            columns=2
        )
    )
})
