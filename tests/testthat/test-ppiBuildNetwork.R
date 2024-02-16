test_that("zero order network behave as expected", {
    data("exampleDESeqResults")

    exNetworkZero <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResults[[1]],
        filterInput=TRUE,
        order="zero"
    )

    expect_length(exNetworkZero, 497)

    expect_equal(
        nrow(as_tibble(tidygraph::activate(exNetworkZero, "edges"))),
        997
    )

    expect_contains(
        colnames(as_tibble(exNetworkZero)),
        c(
            "name",
            "degree",
            "betweenness",
            "seed",
            "hubScoreBtw",
            "hgncSymbol"
        )
    )
})

test_that("simple minimum order networks behave as expected", {
    data("exampleDESeqResults")

    exNetworkSimple <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResults[[1]],
        filterInput=TRUE,
        order="minSimple"
    )

    expect_length(exNetworkSimple, 3960)

    expect_equal(
        nrow(as_tibble(tidygraph::activate(exNetworkSimple, "edges"))),
        15822
    )

    expect_contains(
        colnames(as_tibble(exNetworkSimple)),
        c(
            "name",
            "degree",
            "betweenness",
            "seed",
            "hubScoreBtw",
            "hgncSymbol"
        )
    )
})
