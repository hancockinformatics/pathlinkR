test_that("zero order network behave as expected", {
    exNetworkZero <- ppiBuildNetwork(
        deseqResults=exampleDESeqResults[[1]],
        filterInput=TRUE,
        order="zero"
    )

    expect_length(exNetworkZero, 497)

    expect_equal(
        nrow(as_tibble(tidygraph::activate(exNetworkZero, "edges"))),
        997
    )

    expect_equal(
        colnames(as_tibble(exNetworkZero)),
        c(
            "name",
            "degree",
            "betweenness",
            "seed",
            "hubScoreBtw",
            "hgncSymbol",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj"
        )
    )
})

test_that("simple minimum order networks behave as expected", {
    exNetworkSimple <- ppiBuildNetwork(
        deseqResults=exampleDESeqResults[[1]],
        filterInput=TRUE,
        order="minSimple"
    )

    expect_length(exNetworkSimple, 3960)

    expect_equal(
        nrow(as_tibble(tidygraph::activate(exNetworkSimple, "edges"))),
        15824
    )

    expect_equal(
        colnames(as_tibble(exNetworkSimple)),
        c(
            "name",
            "degree",
            "betweenness",
            "seed",
            "hubScoreBtw",
            "hgncSymbol",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj"
        )
    )
})

test_that("Steinder-trimmed networks behave as expected", {
    set.seed(1)

    exNetworkSteiner <- ppiBuildNetwork(
        deseqResults=exampleDESeqResults[[1]],
        filterInput=TRUE,
        order="minSteiner"
    )

    expect_length(exNetworkSteiner, 1376)

    expect_equal(
        nrow(as_tibble(tidygraph::activate(exNetworkSteiner, "edges"))),
        1375
    )

    expect_equal(
        colnames(as_tibble(exNetworkSteiner)),
        c(
            "name",
            "degree",
            "betweenness",
            "seed",
            "hubScoreBtw",
            "hgncSymbol",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj"
        )
    )
})
