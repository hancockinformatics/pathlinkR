test_that("zero order network behave as expected", {
    data("exampleDESeqResultsHS")

    exNetworkZero <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResultsHS[[1]],
        species="human",
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
    data("exampleDESeqResultsHS")

    suppressMessages(
        exNetworkSimple <- ppiBuildNetwork(
            rnaseqResult=exampleDESeqResultsHS[[1]],
            species="human",
            filterInput=TRUE,
            order="minSimple"
        )
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

test_that("zero order network behave as expected", {
    data("exampleDESeqResultsMM", "innateDbPPIMM")

    exNetworkZero <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResultsMM[[1]],
        species="mouse",
        filterInput=TRUE,
        order="zero",
        ppiData=innateDbPPIMM
    )

    expect_length(exNetworkZero, 126)

    expect_equal(
        nrow(as_tibble(tidygraph::activate(exNetworkZero, "edges"))),
        185
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
    data("exampleDESeqResultsMM", "innateDbPPIMM")

    suppressMessages(
        exNetworkSimple <- ppiBuildNetwork(
            rnaseqResult=exampleDESeqResultsMM[[1]],
            species="mouse",
            filterInput=TRUE,
            order="minSimple",
            ppiData=innateDbPPIMM
        )
    )


    expect_length(exNetworkSimple, 797)

    expect_equal(
        nrow(as_tibble(tidygraph::activate(exNetworkSimple, "edges"))),
        1799
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
