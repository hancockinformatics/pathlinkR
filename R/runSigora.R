#' INTERNAL Wrapper around Sigora's enrichment function
#'
#' @param enrichGenes Vector of genes to enrich
#' @param gpsRepo GPS object to use for testing pathways
#' @param pValFilter Desired threshold for filtering results
#'
#' @return Data frame of results from Sigora
#'
#' @importFrom dplyr %>% across arrange bind_rows contains distinct everything
#'   filter group_by left_join mutate select summarise ungroup
#' @importFrom purrr imap
#' @importFrom sigora sigora
#' @importFrom stringr str_count
#' @importFrom tibble as_tibble tibble
#'
#' @description Internal wrapper function to run Sigora and return the results
#'   with desired columns
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
.runSigora <- function(
        enrichGenes,
        gpsRepo,
        pValFilter=NA
) {

    stopifnot(
        "Input must be a character vector" = { is(enrichGenes, "character") }
    )
    stopifnot(
        "Input must be a character vector" = { is(enrichGenes, "vector") }
    )

    stopifnot(
        "Your input vector doesn't look like Ensembl genes." = {
            any(grepl(pattern="^ENSG", enrichGenes))
        }
    )

    data_env <- new.env(parent=emptyenv())
    data("idmap", "reaH", envir=data_env, package="sigora")
    idmap <- data_env[["idmap"]]
    reaH <- data_env[["reaH"]]

    if (gpsRepo == "default") gpsRepo <- reaH

    ## Run SIGORA based on default settings (GPSrepo=reaH, level=4)
    invisible(capture.output(
        sigoraResult1 <- sigora(
            GPSrepo=gpsRepo,
            level=4,
            markers=TRUE,
            queryList=enrichGenes
        )
    ))

    sigoraResult2 <-
        if (is.na(pValFilter)) {
            sigoraResult1$summary_results
        } else {
            filter(sigoraResult1$summary_results, Bonferroni < pValFilter)
        }

    ## Start by calculating n
    nGenes <- length(enrichGenes[enrichGenes %in% idmap$Ensembl.Gene.ID])

    ## Get the DE genes that were enriched for each pathway
    sigoraDetailedList1 <- sigoraResult1$detailed_results %>%
        as_tibble() %>%
        select(pathway, contains("gene")) %>%
        split(x=., f=.$pathway)

    sigoraDetailedList2 <- sigoraDetailedList1 %>% imap(
        ~tibble(
            pathwy.id=.y,
            EntrezGene.ID=unique(c(pull(.x, gene1), pull(.x, gene2)))
        ) %>%
            left_join(idmap, by="EntrezGene.ID", multiple="all") %>%
            filter(Ensembl.Gene.ID %in% enrichGenes, Symbol != "^$") %>%
            select(pathwy.id, Symbol) %>%
            mutate(across(everything(), as.character)) %>%
            distinct() %>%
            arrange(pathwy.id, Symbol)
    ) %>%
        bind_rows() %>%
        group_by(pathwy.id) %>%
        summarise(genes=paste0(Symbol, collapse=";")) %>%
        mutate(numCandidateGenes=1 + str_count(genes, ";")) %>%
        ungroup() %>%
        mutate(geneRatio=numCandidateGenes / nGenes)

    left_join(
        sigoraResult2,
        sigoraDetailedList2,
        by="pathwy.id",
        multiple="all"
    ) %>%
        mutate(numBgGenes=nGenes) %>%
        select(
            "pathwayId"=pathwy.id,
            "pathwayName"=description,
            "pValue"=pvalues,
            "pValueAdjusted"=Bonferroni,
            genes,
            numCandidateGenes,
            numBgGenes,
            geneRatio
        ) %>% as_tibble()
}
