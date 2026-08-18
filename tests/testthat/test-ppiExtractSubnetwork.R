test_that("subnetwork extraction works with a pathway name", {
    data("exampleDESeqResults", "mappingFileHS")

    exNetwork <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResults[[1]],
        species="human",
        filterInput=TRUE,
        order="zero"
    )

    exPathways <- ppiEnrichNetwork(
        network=exNetwork,
        species="human",
        analysis="hallmark"
    )

    exSubnetwork <- ppiExtractSubnetwork(
        network=exNetwork,
        species="human",
        pathwayEnrichmentResult=exPathways,
        pathwayToExtract="INTERFERON ALPHA RESPONSE"
    )

    expect_equal(nrow(as_tibble(exSubnetwork)), 74)
})

test_that("subnetwork extraction works with a character vector of genes", {
    data("exampleDESeqResults", "mappingFileHS")

    exNetwork2 <- ppiBuildNetwork(
        rnaseqResult=exampleDESeqResults[[1]],
        species="human",
        filterInput=TRUE,
        order="zero"
    )

    exPathways2 <- ppiEnrichNetwork(
        network=exNetwork2,
        species="human",
        analysis="hallmark"
    )

    myGenes <- mappingFileHS %>%
        filter(hgncSymbol %in% unlist(strsplit(exPathways2[[2, 5]], ";"))) %>%
        pull(ensemblGeneId)

    exSubnetwork2 <- ppiExtractSubnetwork(
        network=exNetwork2,
        species="human",
        genes=myGenes
    )

    expect_equal(nrow(as_tibble(exSubnetwork2)), 74)

    expect_error(
        ppiExtractSubnetwork(
            network=exNetwork2,
            species="human",
            genes=list(myGenes)
        )
    )
})
