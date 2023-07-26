test_that("fold change heatmaps are working", {
    set.seed(1)

    vdiffr::expect_doppelganger(
        "plotFoldChangeExample",
        plotFoldChange(
            deseqExampleList,
            pathName = "Generation of second messenger molecules"
        )
    )
})
