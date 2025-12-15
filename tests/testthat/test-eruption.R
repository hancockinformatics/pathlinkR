test_that("a volcano plot works with some of customizations", {
    data("exampleDESeqResults", "sigoraDatabase")

    interferonGenes <- sigoraDatabase %>%
        filter(pathwayName == "Interferon Signaling") %>%
        pull(ensemblGeneId)

    expect_no_error(
        eruption(
            rnaseqResult=exampleDESeqResults[[1]],
            xaxis=c(-4, 4),
            yaxis=c(0, 8),
            highlightGenes=interferonGenes
        )
    )
})
