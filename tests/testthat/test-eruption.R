test_that("the volcano plot is correct", {
    set.seed(1)

    vdiffr::expect_doppelganger(
        "volcanoPlotExample",
        eruption(deseqResults = deseqExampleList[[1]])
    )
})