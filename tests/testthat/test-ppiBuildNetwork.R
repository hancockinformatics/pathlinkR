exNetworkZero <- ppiBuildNetwork(
    deseqResults = deseqExampleList[[1]],
    filterInput = TRUE,
    order = "zero"
)

exNetworkSimple <- ppiBuildNetwork(
    deseqResults = deseqExampleList[[1]],
    filterInput = TRUE,
    order = "minSimple"
)

exNetworkSteiner <- ppiBuildNetwork(
    deseqResults = deseqExampleList[[1]],
    filterInput = TRUE,
    order = "minSteiner"
)

expectedColNames <- c(
    "name", "degree", "betweenness", "seed", "hubScoreBtw", "hgncSymbol",
    "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
)

test_that("we get the right number of nodes and edges", {
    expect_length(exNetworkZero, 497)
    expect_equal(nrow(as_tibble(tidygraph::activate(exNetworkZero, "edges"))), 997)

    expect_length(exNetworkSimple, 3960)
    expect_equal(
        nrow(as_tibble(tidygraph::activate(exNetworkSimple, "edges"))),
        15824
    )

    expect_length(exNetworkSteiner, 6939)
    expect_equal(
        nrow(as_tibble(tidygraph::activate(exNetworkSteiner, "edges"))),
        18944
    )
})

test_that("we have the right column names", {
    expect_equal(colnames(as_tibble(exNetworkZero)), expectedColNames)
    expect_equal(colnames(as_tibble(exNetworkSimple)), expectedColNames)
    expect_equal(colnames(as_tibble(exNetworkSteiner)), expectedColNames)
})
