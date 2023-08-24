test_that("the output looks right", {

    tj1 <- jsonlite::read_json(
        system.file(
            "extdata/networkAnalystExample.json",
            package="pathlinkR"
        ),
        simplifyVector=TRUE
    )

    tj2 <- graph_from_data_frame(
        d=select(tj1$edges, source, target),
        directed=FALSE,
        vertices=select(
            tj1$nodes,
            id,
            label,
            x,
            y,
            "types"=molType,
            expr
        )
    )

    tj3 <- ppiCleanNetwork(as_tbl_graph(tj2))

    tj4 <- as_tibble(tj3)

    expect_length(tj3, 146)
    expect_equal(ncol(tj4), 11)

    expectedColnames <- c(
        "name",
        "degree",
        "betweenness",
        "seed",
        "hubScoreBtw",
        "geneName"
    )

    expect_true(all(expectedColnames %in% colnames(tj4)))
})
