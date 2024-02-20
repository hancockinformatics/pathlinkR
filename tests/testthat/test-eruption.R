test_that("a volcano plot works with some of customizations", {
    set.seed(1)

    data("exampleDESeqResults", "sigoraDatabase")

    interferonGenes <- sigoraDatabase %>%
        filter(pathwayName == "Interferon Signaling") %>%
        pull(ensemblGeneId)

    vdiffr::expect_doppelganger(
        "volcanoPlotExample",
        eruption(
            rnaseqResult=exampleDESeqResults[[1]],
            xaxis=c(-4, 4),
            yaxis=c(0, 8),
            highlightGenes=interferonGenes
        )
    )
})
