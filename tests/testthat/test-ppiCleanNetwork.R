test_that("the output looks right", {

    tj1 <- jsonlite::read_json(
        system.file("extdata/networkAnalystExample.json", package = "pathnet"),
        simplifyVector = TRUE
    )

    tj2 <- igraph::graph_from_data_frame(
        d = dplyr::select(tj1$edges, source, target),
        directed = FALSE,
        vertices = dplyr::select(
            tj1$nodes,
            id,
            label,
            x,
            y,
            "types" = molType,
            expr
        )
    )

    tj3 <- ppiCleanNetwork(tidygraph::as_tbl_graph(tj2))

    tj4 <- dplyr::as_tibble(tj3)

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
