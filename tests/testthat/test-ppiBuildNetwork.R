exNetwork <- ppiBuildNetwork(
    deseqResults = deseqExampleList[[1]],
    filterInput = TRUE,
    order = "zero"
)

exNetworkTibble <- as_tibble(exNetwork)

expectedColNames <- c(
    "name", "degree", "betweenness", "seed", "hubScoreBtw", "hgncSymbol",
    "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
)

test_that("we get the right number of nodes and edges", {
    expect_length(exNetwork, 497)
    expect_equal(nrow(as_tibble(tidygraph::activate(exNetwork, "edges"))), 997)
})

test_that("we have the right column names", {
    expect_equal(colnames(exNetworkTibble), expectedColNames)
})
