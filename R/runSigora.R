#' INTERNAL .runSigora
#'
#' @param enrichGenes Vector of genes to enrich
#' @param gpsRepo GPS object to use for testing pathways
#' @param pValFilter Desired threshold for filtering results
#'
#' @return Data frame of results from Sigora
#'
#' @import dplyr
#' @import sigora
#' @importFrom purrr map imap
#'
#' @description Internal wrapper function to run Sigora and return the results
#'   with desired columns
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
.runSigora <- function(
        enrichGenes,
        gpsRepo,
        pValFilter = NA
) {

    # Run SIGORA based on default settings (GPSrepo = reaH, level = 4)
    invisible(capture.output(
        sigoraResult1 <- sigora(
            GPSrepo = gpsRepo,
            level = 4,
            markers = TRUE,
            queryList = enrichGenes
        )
    ))

    sigoraResult2 <-
        if (is.na(pValFilter)) {
            sigoraResult1$summary_results
        } else {
            filter(sigoraResult1$summary_results, Bonferroni < pValFilter)
        }

    # Start by calculating n
    nGenes <- length(enrichGenes[enrichGenes %in% idmap$Ensembl.Gene.ID])

    # Get the DE genes that were enriched for each pathway
    sigoraDetailedList1 <- sigoraResult1$detailed_results %>%
        as_tibble() %>%
        select(pathway, contains("gene")) %>%
        split(x = ., f = .$pathway)

    sigoraDetailedList2 <- sigoraDetailedList1 %>% imap(
        ~tibble(
            pathwy.id = .y,
            EntrezGene.ID = unique(c(pull(.x, gene1), pull(.x, gene2)))
        ) %>%
            left_join(idmap, by = "EntrezGene.ID", multiple = "all") %>%
            filter(Ensembl.Gene.ID %in% enrichGenes, Symbol != "^$") %>%
            select(pathwy.id, Symbol) %>%
            mutate(across(everything(), as.character)) %>%
            distinct() %>%
            arrange(pathwy.id, Symbol)
    ) %>%
        bind_rows() %>%
        group_by(pathwy.id) %>%
        summarize(genes = paste0(Symbol, collapse = ";")) %>%
        mutate(numCandidateGenes = 1 + str_count(genes, ";")) %>%
        ungroup() %>%
        mutate(geneRatio = numCandidateGenes / nGenes)

    left_join(
        sigoraResult2,
        sigoraDetailedList2,
        by = "pathwy.id",
        multiple = "all"
    ) %>%
        mutate(numBgGenes = nGenes) %>%
        select(
            "pathwayId" = pathwy.id,
            "pathwayName" = description,
            "pValue" = pvalues,
            "pValueAdjusted" = Bonferroni,
            genes,
            numCandidateGenes,
            numBgGenes,
            geneRatio
        ) %>% as_tibble()
}
