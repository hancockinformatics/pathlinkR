test_that("a basic volcano plot is correct", {
    set.seed(1)

    vdiffr::expect_doppelganger(
        "volcanoPlotExample",
        eruption(deseqResult=deseqExampleList[[1]])
    )
})

test_that("some of the options are working as expected", {
    set.seed(1)

    interferonGenes <- sigoraDatabase %>%
        filter(pathwayName == "Interferon Signaling") %>%
        pull(ensemblGeneId)

    vdiffr::expect_doppelganger(
        "volcanoPlotExample2",
        eruption(
            deseqResult=deseqExampleList[[1]],
            xaxis=c(-4, 4),
            yaxis=c(0, 8),
            highlightGenes=interferonGenes
        )
    )
})
