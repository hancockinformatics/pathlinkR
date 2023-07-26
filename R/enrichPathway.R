#' Test lists of genes for enriched pathways
#'
#' @param inputList list of data frames of DESeq2 results. The list names
#'   are used as the comparison for each dataframe (e.g. COVID vs Healthy). Data
#'   frames must have Ensembl gene IDs as the rownames.
#' @param filterInput If providing list of data frames containing the
#'   unfiltered output from `DESeq2::results()`, set this to TRUE to filter for
#'   DE genes using the thresholds set by the `pCutoff` and `fcCutoff`
#'   arguments. When FALSE it's assumed your passing the filtered
#'   results into `inputList` and no more filtering will be done.
#' @param pCutoff Adjusted p value cutoff, defaults to 0.05.
#' @param fcCutoff Minimum absolute fold change, defaults to 1.5.
#' @param split Boolean (TRUE); Split into up and down-regulated DEGs using the
#'   column "log2FoldChange", and do enrichment separately
#' @param analysis Default is "sigora", but can also be "reactomepa" or
#'   "hallmark"
#' @param filterResults Should the output be filtered for significance? Use
#'   `1` to return the unfiltered results, or any number less than 1 for a
#'   custom p-value cutoff. If `default`, the significance cutoff for Sigora
#'   is 0.001, and for ReactomePA or Hallmark is 0.05.
#' @param gpsRepo Only applies to `analysis = "sigora"`. Gene Pair Signature
#'   object for Sigora to use to test for enriched pathways. We recommend using
#'   the one which ships with Sigora, which is already loaded as "reaH".
#' @param geneUniverse Only applies when `analysis` is "reactomepa" or
#'   "hallmark". The set of background genes to use when testing with ReactomePA
#'   or Hallmark gene sets. For ReactomePA this must be a character vector of
#'   Entrez genes. For Hallmark, it must be Ensembl IDs.
#'
#' @return A data frame of pathway enrichment results for all input comparisons
#' @export
#'
#' @import dplyr
#' @import purrr
#' @importFrom clusterProfiler enricher
#'
#' @description This function provides a simple and consistent interface to
#'   three different pathway enrichment tools: Sigora and ReactomePA (which both
#'   test for Reactome pathways), and MSigDB Hallmark gene set enrichment. The
#'   input must be a named list of data frames, which can be pre-filtered or
#'   "raw", in which case the function can filter using user-defined cutoffs.
#'   Column names are expected to comply with those output by
#'   `DESeq2::results()` function, namely: padj and log2FoldChange. Rownames are
#'   assumed to contain the input genes to be tested.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' enrichPathway(
#'     inputList = deseqExampleList[1],
#'     filterInput = TRUE,
#'     split = TRUE,
#'     analysis = "reactomepa",
#'     filterResults = "default",
#' )
#'
enrichPathway <- function(
        inputList,
        filterInput = TRUE,
        pCutoff = 0.05,
        fcCutoff = 1.5,
        split = TRUE,
        analysis = "sigora",
        filterResults = "default",
        gpsRepo = reaH,
        geneUniverse = NULL
) {

    ## Check inputs
    stopifnot(analysis %in% c("sigora", "reactomepa", "hallmark"))

    stopifnot(
        "Provide a named list of data frames of results, with the name
        of each item in the list as the comparison name."  = {
            is.list(inputList)
            !is.null(names(inputList))
            all(unlist(lapply(inputList, is.data.frame)))
        }
    )

    ## Iterate through each element of "inputList"
    resultList <- imap(inputList, function(x, comparison) {

        ## Basic checks before continuing
        stopifnot(
            "Elements of 'inputList' should be named" = !is.null(comparison)
        )

        if (all(rownames(x) == seq_len(nrow(x)))) {
            stop(
                "The rownames of the data frame for the element with name '",
                comparison, "' don't look like gene IDs!"
            )
        }

        message("Comparison being analyzed: ", comparison)

        ## Filter the input genes if specified
        deseqResults <-
            if (filterInput) {
                stopifnot(
                    c("padj", "log2FoldChange") %in% colnames(x)
                )
                message("\tFiltering the results before testing...")
                filter(
                    x,
                    padj < pCutoff,
                    abs(log2FoldChange) > log2(fcCutoff)
                )
            } else {
                x
            }

        ## Turn the input into a list of gene IDs, split by direction or not
        if (split) {
            preppedGenes <- list(
                "Up"   = rownames(filter(deseqResults, log2FoldChange > 0)),
                "Down" = rownames(filter(deseqResults, log2FoldChange < 0))
            )
            message(
                "\tDEGs used: ",
                length(preppedGenes$Up), " Up, ",
                length(preppedGenes$Down), " Down..."
            )
        } else {
            preppedGenes <- list("All" = rownames(deseqResults))
            message("\tDEGs used: ", length(preppedGenes$All), "...")
        }

        ## Sigora
        if (analysis == "sigora") {
            message("\tRunning enrichment using Sigora...")

            runSigoraSafely <- possibly(.runSigora)

            resultFinal <- imap_dfr(
                .x  = preppedGenes,
                .id = "direction",
                function(y, direction) {
                    runSigoraSafely(
                        enrichGenes = y,
                        gpsRepo = gpsRepo,
                        pValFilter = ifelse(
                            filterResults == "default",
                            0.001,
                            filterResults
                        )
                    )
                }
            )

            resultFinal$totalGenes <- nrow(deseqResults)

            message(
                "\tDone! Found ",
                nrow(resultFinal),
                " enriched pathways.\n"
            )
            return(resultFinal)
        }


        ## ReactomePA or Hallmark
        if (analysis %in% c("reactomepa", "hallmark")) {

            ## ReactomePA
            if (analysis == "reactomepa") {
                message("\tRunning enrichment using ReactomePA")

                rpaHallResult <- imap_dfr(
                    .x  = preppedGenes,
                    .id = "direction",
                    function(y, direction) {

                        genesEntrez <- idmap %>%
                            filter(Ensembl.Gene.ID %in% y) %>%
                            pull(EntrezGene.ID) %>%
                            unique()

                        enricher(
                            genesEntrez,
                            TERM2GENE = select(
                                reactomeDatabase,
                                pathwayId,
                                entrezGeneId
                            ),
                            TERM2NAME = select(
                                reactomeDatabase,
                                pathwayId,
                                pathwayName
                            ),
                            universe = geneUniverse,
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = ifelse(
                                filterResults == "default",
                                0.05,
                                filterResults
                            )
                        ) %>% as_tibble()
                    }
                ) %>%
                    mutate(geneID = as.character(geneID)) %>%
                    separate_longer_delim(geneID, delim = "/") %>%
                    left_join(
                        mappingFile,
                        by = c("geneID" = "entrezGeneId"),
                        multiple = "all"
                    ) %>%
                    select(
                        -any_of(c("geneID", "entrezGeneId", "ensemblGeneId"))
                    ) %>%
                    group_by(ID) %>%
                    mutate(genes = paste(hgncSymbol, collapse = ";")) %>%
                    ungroup() %>%
                    select(-hgncSymbol) %>%
                    distinct()
            }

            ## Hallmark
            if (analysis == "hallmark") {
                message("\tRunning enrichment using Hallmark...")

                rpaHallResult <- imap_dfr(
                    .x  = preppedGenes,
                    .id = "direction",
                    function(y, direction) {

                        enricher(
                            y,
                            TERM2GENE = mSigDbTermToGene,
                            universe = geneUniverse,
                            pvalueCutoff = ifelse(
                                filterResults == "default",
                                0.05,
                                filterResults
                            )
                        ) %>% as_tibble()
                    }
                ) %>%
                    separate_longer_delim(geneID, delim = "/") %>%
                    left_join(
                        mappingFile,
                        by = c("geneID" = "ensemblGeneId"),
                        multiple = "all"
                    ) %>%
                    select(
                        -any_of(c("geneID", "entrezGeneId", "ensemblGeneId"))
                    ) %>%
                    group_by(ID) %>%
                    mutate(genes = paste(hgncSymbol, collapse = ";")) %>%
                    ungroup() %>%
                    select(-hgncSymbol) %>%
                    distinct()
            }

            ## ReactomePA or Hallmark
            resultFinal <- rpaHallResult %>%
                separate_wider_delim(
                    cols = GeneRatio,
                    delim = "/",
                    names = c("numCandidateGenes", "numBgGenes")
                ) %>%
                mutate(
                    across(c(numCandidateGenes, numBgGenes), as.numeric),
                    geneRatio = numCandidateGenes / numBgGenes,
                    totalGenes = nrow(deseqResults)
                ) %>%
                select(
                    direction,
                    "pathwayId" = ID,
                    "pathwayName" = Description,
                    "pValue" = pvalue,
                    "pValueAdjusted" = p.adjust,
                    genes,
                    numCandidateGenes,
                    numBgGenes,
                    geneRatio,
                    totalGenes
                )

            message(
                "\tDone! Found ",
                nrow(resultFinal),
                " enriched pathways.\n"
            )
            return(resultFinal)
        }
    })

    results_all_comparisons <- resultList %>%
        bind_rows(.id = "comparison") %>%
        left_join(
            .,
            select(topPathwaysMore, pathwayId, topPathways),
            by = "pathwayId",
            multiple = "all"
        ) %>%
        as_tibble()

    return(results_all_comparisons)
}
