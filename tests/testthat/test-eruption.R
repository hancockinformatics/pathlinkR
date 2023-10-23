test_that("a basic volcano plot is correct", {
    set.seed(1)

    data("exampleDESeqResults", package = "pathlinkR")

    vdiffr::expect_doppelganger(
        "volcanoPlotExample",
        eruption(
            rnaseqResult=exampleDESeqResults[[1]],
            columnFC="log2FoldChange",
            columnP="padj"
        )
    )
})

test_that("some of the options are working as expected", {
    set.seed(1)

    data(
        list = c("exampleDESeqResults", "sigoraDatabase"),
        package = "pathlinkR"
    )

    interferonGenes <- sigoraDatabase %>%
        filter(pathwayName == "Interferon Signaling") %>%
        pull(ensemblGeneId)

    vdiffr::expect_doppelganger(
        "volcanoPlotExample2",
        eruption(
            rnaseqResult=exampleDESeqResults[[1]],
            columnFC="log2FoldChange",
            columnP="padj",
            xaxis=c(-4, 4),
            yaxis=c(0, 8),
            highlightGenes=interferonGenes
        )
    )
})
