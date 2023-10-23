#' Test significant DE genes from DESeq2 for enriched pathways
#'
#' @param inputList A list, with each element containing RNA-Seq results as a
#'   "DESeqResults", "TopTags", or "data.frame" object. Rownames of each table
#'   must contain Ensembl Gene IDs. The list names are used as the comparison
#'   name for each element (e.g. "COVID vs Healthy"). See Details for more
#'   information on supported input types.
#' @param columnFC Character; optional column name containing fold change
#'   values. Defaults to NA.
#' @param columnP Character; optional column name containing p values. Defaults
#'   to NA.
#' @param filterInput When providing list of data frames containing the
#'   unfiltered output from `DESeq2::results()` (default format), set this to
#'   `TRUE` to filter for significant genes using the thresholds set by the
#'   `pCutoff` and `fcCutoff`. When this argument is `FALSE` it's
#'   assumed your passing a pre-filtered data frame to `inputList`, and no
#'   more filtering will be done.
#' @param pCutoff Adjusted p value cutoff, defaults to <0.05
#' @param fcCutoff Minimum absolute fold change, defaults to >1.5
#' @param split Boolean (TRUE); Split into up and down-regulated DE genes using
#'   the requisite column `columnFC`, and do enrichment independently. Results
#'   are combined at the end, with an added "direction" column.
#' @param analysis Default is "sigora", but can also be "reactomepa" or
#'   "hallmark"
#' @param filterResults Should the output be filtered for significance? Use
#'   `1` to return the unfiltered results, or any number less than 1 for a
#'   custom p-value cutoff. If left as `default`, the significance cutoff
#'   for Sigora is 0.001, or 0.05 for ReactomePA and Hallmark.
#' @param gpsRepo Only applies to `analysis="sigora"`. Gene Pair Signature
#'   object for Sigora to use to test for enriched pathways. Leaving this set
#'   as "default" will use the "reaH" GPS object from `Sigora`, or you can
#'   provide your own custom GPS repository.
#' @param geneUniverse Only applies when `analysis` is "reactomepa" or
#'   "hallmark". The set of background genes to use when testing with ReactomePA
#'   or Hallmark gene sets. For ReactomePA this must be a character vector of
#'   Entrez genes. For Hallmark, it must be Ensembl IDs.
#'
#' @return A "data.frame" (tibble) of pathway enrichment results for all input
#'   comparisons, with the following columns:
#'   \item{comparison}{Source comparison from the names of `inputList`}
#'   \item{direction}{Whether the pathway was enriched in all genes
#'   (`split=FALSE`), or up- or down-regulated genes (`split=TRUE`)}
#'   \item{pathwayId}{Pathway identifier}
#'   \item{pathwayName}{Pathway name}
#'   \item{pValue}{Nominal p value for the pathway}
#'   \item{pValueAdjusted}{p value corrected for multiple testing}
#'   \item{genes}{Candidate genes, which were DE for the comparison and also in
#'   the pathway}
#'   \item{numCandidateGenes}{Number of candidate genes}
#'   \item{numBgGenes}{Number of background genes for the pathway}
#'   \item{geneRatio}{Ratio of candidate and background genes}
#'   \item{totalGenes}{Number of DE genes which were tested for enriched
#'   pathways}
#'   \item{topLevelPathway}{High level Reactome term which serves to group
#'   similar pathways}
#'
#' @export
#'
#' @import dplyr
#'
#' @importFrom purrr imap imap_dfr possibly
#' @importFrom tidyr separate_longer_delim separate_wider_delim
#'
#' @description This function provides a simple and consistent interface to
#'   three different pathway enrichment tools: Sigora and ReactomePA (which both
#'   test for Reactome pathways), and MSigDB Hallmark gene set enrichment.
#'
#' @details `inputList` must be a named list of RNA-Seq results, with each
#'   element being of class "DESeqResults" from `DESeq2`, "TopTags" from
#'   `edgeR`, or a simple data frame. For the first two cases, column names are
#'   expected to be the standard defined by each class ("log2FoldChange" and
#'   "padj" for "DESeqResults", and "logFC" and "FDR" for "TopTags"). Hence for
#'   these two cases the arguments `columnFC` and `columnP`
#'   can be left as `NA`.
#'
#'   In the last case (elements are "data.frame"), `columnFC` and
#'   `columnP` must be supplied when `filterInput=TRUE`, and
#'   `columnFC` must be given if `split=TRUE`.
#'
#' @references
#'   Sigora: <https://cran.r-project.org/package=sigora>
#'   ReactomePA: <https://www.bioconductor.org/packages/ReactomePA/>
#'   MSigDB/Hallmark: <https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp>
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' data("exampleDESeqResults")
#'
#' pathwayEnrichment(
#'     inputList=exampleDESeqResults[1],
#'     columnFC="log2FoldChange",
#'     columnP="padj",
#'     filterInput=TRUE,
#'     split=TRUE,
#'     analysis="hallmark",
#'     filterResults="default"
#' )
#'
pathwayEnrichment <- function(
        inputList,
        columnFC=NA,
        columnP=NA,
        filterInput=TRUE,
        pCutoff=0.05,
        fcCutoff=1.5,
        split=TRUE,
        analysis="sigora",
        filterResults="default",
        gpsRepo="default",
        geneUniverse=NULL
) {
    stopifnot(analysis %in% c("sigora", "reactomepa", "hallmark"))

    stopifnot(
        "Provide a named list of data frames of results, with the name
        of each item in the list as the comparison name."  = {
            is.list(inputList)
            !is.null(names(inputList))
        }
    )

    ## Coerce the input
    if (is(inputList[[1]], "DESeqResults")) {
        inputListCleaned <- lapply(inputList, function(x) {
            as.data.frame(x) %>%
                rename("LogFoldChange"=log2FoldChange, "PAdjusted"=padj)
        })

    } else if (is(inputList[[1]], "TopTags")) {
        inputListCleaned <- lapply(inputList, function(x) {
            as.data.frame(x) %>%
                rename("LogFoldChange"=logFC, "PAdjusted"=FDR)
        })

    } else {

        if (filterInput) {
            stopifnot(
                "Must provide 'columnFC' and 'columnP' when filtering input" = {
                    !any(is.na(columnFC), is.na(columnP))
                }
            )
            inputList <- lapply(inputList, function(x) {
                rename(
                    x,
                    "LogFoldChange"=any_of(columnFC),
                    "PAdjusted"=any_of(columnP)
                )
            })
        }
        if (split) {
            stopifnot(
                "Must provide 'columnFC' when splitting input" = {
                    !is.na(columnFC)
                }
            )
            inputList <- lapply(inputList, function(x) {
                rename(x, "LogFoldChange"=any_of(columnFC))
            })
        }
        inputListCleaned <- inputList
    }

    data_env <- new.env(parent=emptyenv())
    data("idmap", envir=data_env, package="sigora")
    idmap <- data_env[["idmap"]]

    data(
        "pathwayCategories",
        "reactomeDatabase",
        "hallmarkDatabase",
        "mappingFile",
        envir=data_env,
        package = "pathlinkR"
    )
    pathwayCategories <- data_env[["pathwayCategories"]]
    reactomeDatabase <- data_env[["reactomeDatabase"]]
    hallmarkDatabase <- data_env[["hallmarkDatabase"]]
    mappingFile <- data_env[["mappingFile"]]


    ## Iterate through each element of "inputListCleaned"
    resultList <- imap(inputListCleaned, function(x, comparison) {
        stopifnot(
            "Elements of 'inputList' should be named" = !is.null(comparison)
        )

        if (all(rownames(x) == seq_len(nrow(x)))) {
            stop(
                "The rownames of the data frame for the element with name '",
                comparison, "' don't look like Ensembl gene IDs!"
            )
        }

        message("Comparison being analyzed: ", comparison)

        ## Filter the input genes if specified
        rnaseqResults <-
            if (filterInput) {
                # stopifnot(c("padj", "log2FoldChange") %in% colnames(x))
                message("\tFiltering the results before testing...")
                filter(x, PAdjusted < pCutoff, abs(LogFoldChange) > log2(fcCutoff))
            } else {
                x
            }

        ## Turn the input into a list of gene IDs, split by direction or not
        if (split) {
            preppedGenes <- list(
                "Up"=rownames(filter(rnaseqResults, LogFoldChange > 0)),
                "Down"=rownames(filter(rnaseqResults, LogFoldChange < 0))
            )
            message(
                "\tDEGs used: ",
                length(preppedGenes$Up), " Up, ",
                length(preppedGenes$Down), " Down..."
            )
        } else {
            preppedGenes <- list("All"=rownames(rnaseqResults))
            message("\tDEGs used: ", length(preppedGenes$All), "...")
        }

        ## Sigora
        if (analysis == "sigora") {
            message("\tRunning enrichment using Sigora...")

            runSigoraSafely <- possibly(.runSigora)

            resultFinal <- imap_dfr(
                .x =preppedGenes,
                .id="direction",
                function(y, direction) {
                    runSigoraSafely(
                        enrichGenes=y,
                        gpsRepo=gpsRepo,
                        pValFilter=ifelse(
                            filterResults == "default",
                            0.001,
                            filterResults
                        )
                    )
                }
            )
            resultFinal$totalGenes <- nrow(rnaseqResults)

            message(
                "\tDone, found ", nrow(resultFinal), " enriched pathways.\n"
            )
            return(resultFinal)
        }


        ## ReactomePA or Hallmark
        if (analysis %in% c("reactomepa", "hallmark")) {

            ## ReactomePA
            if (analysis == "reactomepa") {
                message("\tRunning enrichment using ReactomePA")

                rpaHallResult <- imap_dfr(
                    .x =preppedGenes,
                    .id="direction",
                    function(y, direction) {

                        genesEntrez <- idmap %>%
                            filter(Ensembl.Gene.ID %in% y) %>%
                            pull(EntrezGene.ID) %>%
                            unique()

                        tibble::as_tibble(clusterProfiler::enricher(
                            genesEntrez,
                            TERM2GENE=select(
                                reactomeDatabase,
                                pathwayId,
                                entrezGeneId
                            ),
                            TERM2NAME=select(
                                reactomeDatabase,
                                pathwayId,
                                pathwayName
                            ),
                            universe=geneUniverse,
                            minGSSize=10,
                            maxGSSize=500,
                            pvalueCutoff=ifelse(
                                filterResults == "default",
                                0.05,
                                filterResults
                            )
                        ))
                    }
                ) %>%
                    mutate(geneID=as.character(geneID)) %>%
                    separate_longer_delim(geneID, delim="/") %>%
                    left_join(
                        mappingFile,
                        by=c("geneID" = "entrezGeneId"),
                        multiple="all",
                        relationship="many-to-many"
                    ) %>%
                    select(
                        -any_of(c("geneID", "entrezGeneId", "ensemblGeneId"))
                    ) %>%
                    group_by(ID) %>%
                    mutate(genes=paste(hgncSymbol, collapse=";")) %>%
                    ungroup() %>%
                    select(-hgncSymbol) %>%
                    distinct()
            }

            ## Hallmark
            if (analysis == "hallmark") {
                message("\tRunning enrichment using Hallmark...")

                rpaHallResult <- imap_dfr(
                    .x =preppedGenes,
                    .id="direction",
                    function(y, direction) {

                        tibble::as_tibble(clusterProfiler::enricher(
                            y,
                            TERM2GENE=hallmarkDatabase,
                            universe=geneUniverse,
                            pvalueCutoff=ifelse(
                                filterResults == "default",
                                0.05,
                                filterResults
                            )
                        ))
                    }
                ) %>%
                    separate_longer_delim(geneID, delim="/") %>%
                    left_join(
                        mappingFile,
                        by=c("geneID" = "ensemblGeneId"),
                        multiple="all",
                        relationship="many-to-many"
                    ) %>%
                    select(
                        -any_of(c("geneID", "entrezGeneId", "ensemblGeneId"))
                    ) %>%
                    group_by(ID) %>%
                    mutate(genes=paste(hgncSymbol, collapse=";")) %>%
                    ungroup() %>%
                    select(-hgncSymbol) %>%
                    distinct()
            }

            ## ReactomePA or Hallmark
            resultFinal <- rpaHallResult %>%
                separate_wider_delim(
                    cols=GeneRatio,
                    delim="/",
                    names=c("numCandidateGenes", "numBgGenes")
                ) %>%
                mutate(
                    across(c(numCandidateGenes, numBgGenes), as.numeric),
                    geneRatio=numCandidateGenes / numBgGenes,
                    totalGenes=nrow(rnaseqResults)
                ) %>%
                select(
                    direction,
                    "pathwayId"=ID,
                    "pathwayName"=Description,
                    "pValue"=pvalue,
                    "pValueAdjusted"=p.adjust,
                    genes,
                    numCandidateGenes,
                    numBgGenes,
                    geneRatio,
                    totalGenes
                )

            message(
                "\tDone, found ", nrow(resultFinal), " enriched pathways.\n"
            )
            return(resultFinal)
        }
    })

    results_all_comparisons <- resultList %>%
        bind_rows(.id="comparison") %>%
        left_join(
            select(pathwayCategories, pathwayId, topLevelPathway),
            by="pathwayId",
            multiple="all"
        ) %>%
        tibble::as_tibble()

    return(results_all_comparisons)
}
