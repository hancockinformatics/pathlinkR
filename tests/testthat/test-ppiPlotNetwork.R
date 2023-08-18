test_that("we get the right plot output", {
    set.seed(1)

    exNetwork <- ppiBuildNetwork(
        deseqResults=deseqExampleList[[1]],
        filterInput=TRUE,
        order="zero"
    )

    vdiffr::expect_doppelganger(
        "example-network",
        ppiPlotNetwork(
            exNetwork,
            fillColumn=log2FoldChange,
            fillType="foldChange",
            legend=TRUE,
            label=FALSE
        )
    )
})
