test_that("pathway plots are correct", {
    set.seed(1)

    data("sigoraExamples")

    vdiffr::expect_doppelganger(
        "pathwayPlotExample",
        pathwayPlots(sigoraExamples, columns=2)
    )
})

test_that("pathwayPlots works with fgsea results", {
    set.seed(1)

    data("sigoraExamples")

    uniquePathways <- unique(as.character(sigoraExamples[["pathwayName"]]))

    fgseaInput <- data.frame(
        comparison = rep(c("A", "B"), each = length(uniquePathways) / 2),
        pathway = uniquePathways,
        pval = abs(rnorm(n = length(uniquePathways), mean = 0.0001, sd = 0.001)),
        log2err = NA,
        ES = NA,
        NES = rnorm(n = length(uniquePathways), mean = 0, sd = 1),
        size = NA,
        leadingEdge = NA
    )

    fgseaInput$padj <- p.adjust(fgseaInput$pval, method = "BH")
    fgseaInput <- fgseaInput[fgseaInput$padj < 0.01, ]

    vdiffr::expect_doppelganger(
        "pathwayPlotExampleFGSEA",
         pathwayPlots(fgseaInput, columns=2)
    )
})
