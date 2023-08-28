test_that("subnetwork extraction works with a pathway name", {
    exNetwork <- ppiBuildNetwork(
        deseqResults=deseqExampleList[[1]],
        filterInput=TRUE,
        order="zero"
    )

    exPathways <- ppiEnrichNetwork(
        network=exNetwork,
        analysis="hallmark"
    )

    exSubnetwork <- ppiExtractSubnetwork(
        network=exNetwork,
        pathwayEnrichmentResult=exPathways,
        pathwayToExtract="INTERFERON ALPHA RESPONSE"
    )

    expect_equal(nrow(as_tibble(exSubnetwork)), 74)
})

test_that("subnetwork extraction works with a list of genes", {
    exNetwork2 <- ppiBuildNetwork(
        deseqResults=deseqExampleList[[1]],
        filterInput=TRUE,
        order="zero"
    )

    exPathways2 <- ppiEnrichNetwork(
        network=exNetwork2,
        analysis="hallmark"
    )

    myGenes <- mappingFile %>%
        filter(hgncSymbol %in% unlist(strsplit(exPathways2[[2, 7]], ";"))) %>%
        pull(ensemblGeneId)

    exSubnetwork2 <- ppiExtractSubnetwork(
        network=exNetwork2,
        genes=myGenes
    )

    expect_equal(nrow(as_tibble(exSubnetwork2)), 74)
})
