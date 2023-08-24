test_that("subnetwork extraction works with a pathway name", {
    exNetwork <- ppiBuildNetwork(
        deseqResults=deseqExampleList[[1]],
        filterInput=TRUE,
        order="zero"
    )

    exPathways <- ppiEnrichNetwork(
        network=exNetwork,
        analysis="sigora"
    )

    exSubnetwork <- ppiExtractSubnetwork(
        network=exNetwork,
        pathwayEnrichmentResult=exPathways,
        pathwayToExtract="Interferon alpha/beta signaling"
    )

    expect_equal(nrow(as_tibble(exSubnetwork)), 27)
})

test_that("subnetwork extraction works with a list of genes", {
    exNetwork2 <- ppiBuildNetwork(
        deseqResults=deseqExampleList[[1]],
        filterInput=TRUE,
        order="zero"
    )

    exPathways2 <- ppiEnrichNetwork(
        network=exNetwork2,
        analysis="sigora"
    )

    myGenes <- mappingFile %>%
        filter(hgncSymbol %in% unlist(strsplit(exPathways2[[1, 7]], ";"))) %>%
        pull(ensemblGeneId)

    exSubnetwork2 <- ppiExtractSubnetwork(
        network=exNetwork2,
        genes=myGenes
    )

    expect_equal(nrow(as_tibble(exSubnetwork2)), 27)
})
