#' Test significant DE genes for enriched pathways
#'
#' @param inputList A list, with each element containing RNA-Seq results as a
#'   "DESeqResults", "TopTags", or "data.frame" object. Rownames of each table
#'   must contain Ensembl Gene IDs. The list names are used as the comparison
#'   name for each element (e.g. "COVID vs Healthy"). See Details for more
#'   information on supported input types.
#' @param columnFC Character; Column to plot along the x-axis, typically log2
#'   fold change values. Only required when `rnaseqResult` is a simple data
#'   frame. Defaults to NA.
#' @param columnP Character; Column to plot along the y-axis, typically nominal
#'   or adjusted p values. Only required when `rnaseqResult` is a simple data
#'   frame. Defaults to NA.
#' @param filterInput When providing list of data frames containing the
#'   unfiltered RNA-Seq results (i.e. not all genes are significant), set this
#'   to `TRUE` to remove non-significant genes using the thresholds set by the
#'   `pCutoff` and `fcCutoff`. When this argument is `FALSE`
#'   its assumed your passing a pre-filtered data in `inputList`, and no
#'   more filtering will be done.
#' @param pCutoff Adjusted p value cutoff when filtering. Defaults to < 0.05.
#' @param fcCutoff Minimum absolute fold change value when filtering. Defaults
#'   to > 1.5
#' @param split Boolean (TRUE); Split into up- and down-regulated DE genes using
#'   the fold change column, and do enrichment independently on each. Results
#'   are combined at the end, with an added "direction" column.
#' @param analysis Method/database to use for enrichment analysis. The default
#'   is "sigora", but can also be "reactome"/"reactomepa", "hallmark", "kegg",
#'   "fgsea_reactome" or "fgsea_hallmark".
#' @param filterResults Should the output be filtered for significance? Use `1`
#'   to return the unfiltered results, or any number less than 1 for a custom
#'   p-value cutoff. If left as `default`, the significance cutoff for
#'   `analysis="sigora"` is 0.001, or 0.05 for "reactome", "hallmark", and
#'   "kegg".
#' @param gpsRepo Only applies to `analysis="sigora"`. Gene Pair Signature (GPS)
#'   object for Sigora to use to test for enriched pathways. "reaH" (default)
#'   will use the Reactome GPS object from `Sigora`; "kegH" will use the KEGG
#'   GPS. One can also provide their own GPS object; see Sigora's documentation
#'   for details.
#' @param gpsLevel Only applies to `analysis="sigora"`. If left as `default`,
#'   will be set to `4` for `gpsRepo="reaH"` or `2` for `gpeRepo="kegH"`. If
#'   providing your own GPS object, can be set as desired; see Sigora's
#'   documentation for details.
#' @param geneUniverse Only applies when `analysis` is "reactome"/"reactomepa",
#'   "hallmark", or "kegg". The set of background genes to use when testing with
#'   Reactome, Hallmark, or KEGG gene sets. For Reactome this must be a
#'   character vector of Entrez genes. For Hallmark or KEGG, it must be Ensembl
#'   IDs.
#' @param verbose Logical; If FALSE (the default), don't print info/progress
#'   messages.
#'
#' @return A "data.frame" (tibble) of pathway enrichment results for all input
#'   comparisons, with the following columns:
#'   \item{comparison}{Source comparison from the names of `inputList`}
#'   \item{direction}{Whether the pathway was enriched in all genes
#'   (`split=FALSE`), or up- or down-regulated genes (`split=TRUE`)}
#'   \item{pathwayId}{Pathway identifier}
#'   \item{pathwayName}{Pathway name}
#'   \item{pValue}{Nominal p value for the pathway}
#'   \item{pValueAdjusted}{p value, corrected for multiple testing}
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
#'   In the last case (elements are "data.frame"), both `columnFC` and
#'   `columnP` must be supplied when `filterInput=TRUE`,
#'   and `columnFC` must be given if `split=TRUE`.
#'
#'   Setting `analysis` to any of "reactome", "reactomepa", "hallmark", or
#'   "kegg" will execute traditional over-representation analysis, the only
#'   difference being the database used ("reactome" and "reactomepa" are treated
#'   the same). Setting `analysis="sigora"` will use a gene pair-based approach,
#'   which can be performed on either Reactome data when `gpsRepo="reaH"` or
#'   KEGG data with `gpsRepo="kegH"`.
#'
#' @references
#'   Sigora: <https://cran.r-project.org/package=sigora>
#'   ReactomePA: <https://www.bioconductor.org/packages/ReactomePA/>
#'   Reactome: <https://reactome.org/>
#'   MSigDB/Hallmark: <https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp>
#'   KEGG: <https://www.kegg.jp/>
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' data("exampleDESeqResults")
#'
#' pathwayEnrichment(
#'     inputList=exampleDESeqResults[1],
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
        gpsRepo="reaH",
        gpsLevel="default",
        geneUniverse=NULL,
        verbose=FALSE
) {
    stopifnot(
        analysis %in% c(
            "sigora",
            "reactome",
            "reactomepa",
            "hallmark",
            "kegg",
            "fgsea_reactome",
            "fgsea_hallmark"
        )
    )

    stopifnot(
        "Provide a named list of data frames of results, with the name
        of each item in the list as the comparison name." ={
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
                "Must provide 'columnFC' and 'columnP' when filtering input"={
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
                "Must provide 'columnFC' when splitting input of class
                'data.frame'"={
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
        "keggDatabase",
        "mappingFile",
        envir=data_env,
        package="pathlinkR"
    )
    pathwayCategories <- data_env[["pathwayCategories"]]
    reactomeDatabase <- data_env[["reactomeDatabase"]]
    hallmarkDatabase <- data_env[["hallmarkDatabase"]]
    keggDatabase <- data_env[["keggDatabase"]]
    mappingFile <- data_env[["mappingFile"]]


    ## Iterate through each element of "inputListCleaned"
    resultList <- imap(inputListCleaned, function(x, comparison) {
        stopifnot(
            "Elements of 'inputList' should be named"=!is.null(comparison)
        )

        if (all(rownames(x) == seq_len(nrow(x)))) {
            stop(
                "The rownames of the data frame for the element with name '",
                comparison, "' don't look like Ensembl gene IDs!"
            )
        }

        if (verbose) message("Comparison being analyzed: ", comparison)

        ## Filter the input genes if specified
        rnaseqResults <-
            if (filterInput) {
                if (verbose) {
                    message("\tFiltering the results before testing...")
                }
                filter(
                    x,
                    PAdjusted < pCutoff,
                    abs(LogFoldChange) > log2(fcCutoff)
                )
            } else {
                x
            }

        ## Turn the input into a list of gene IDs, split by direction or not
        if (split) {
            preppedGenesTable <- list(
                "Up"=filter(rnaseqResults, LogFoldChange > 0),
                "Down"=filter(rnaseqResults, LogFoldChange < 0)
            )

            if (verbose) {
                message(
                    "\tDEGs used: ",
                    nrow(preppedGenesTable$Up), " Up, ",
                    nrow(preppedGenesTable$Down), " Down..."
                )
            }
        } else {
            preppedGenesTable <- list("All"=rnaseqResults)
            if (verbose) {
                message("\tDEGs used: ", nrow(preppedGenesTable$All), "...")
            }
        }

        ## Sigora
        if (analysis == "sigora") {
            if (verbose) message("\tRunning enrichment using Sigora...")

            runSigoraSafely <- possibly(.runSigora)

            resultFinal <- imap_dfr(
                .x=preppedGenesTable,
                .id="direction",
                function(y, direction) {
                    runSigoraSafely(
                        enrichGenes=rownames(y),
                        gpsRepo=gpsRepo,
                        gpsLevel=gpsLevel,
                        pValFilter=ifelse(
                            filterResults == "default",
                            0.001,
                            filterResults
                        )
                    )
                }
            )
            resultFinal$totalGenes <- nrow(rnaseqResults)

            if (verbose) {
                message(
                    "\tDone, found ", nrow(resultFinal), " enriched pathways.\n"
                )
            }
            return(resultFinal)
        }


        ## Reactome, Hallmark, or KEGG
        if (analysis %in% c("reactome", "reactomepa", "hallmark", "kegg")) {

            ## ReactomePA
            if (analysis %in% c("reactomepa", "reactome")) {
                if (verbose) message("\tRunning enrichment using Reactome")

                oraResult <- imap_dfr(
                    .x=preppedGenesTable,
                    .id="direction",
                    function(y, direction) {

                        genesEntrez <- idmap %>%
                            filter(Ensembl.Gene.ID %in% rownames(y)) %>%
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
                if (verbose) message("\tRunning enrichment using Hallmark...")

                oraResult <- imap_dfr(
                    .x=preppedGenesTable,
                    .id="direction",
                    function(y, direction) {
                        tibble::as_tibble(clusterProfiler::enricher(
                            rownames(y),
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

            ## KEGG
            if (analysis == "kegg") {
                if (verbose) message("\tRunning enrichment using KEGG")

                oraResult <- imap_dfr(
                    .x=preppedGenesTable,
                    .id="direction",
                    function(y, direction) {
                        tibble::as_tibble(clusterProfiler::enricher(
                            rownames(y),
                            TERM2GENE=select(
                                keggDatabase,
                                pathwayId,
                                ensemblGeneId
                            ),
                            TERM2NAME=select(
                                keggDatabase,
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


            ## Reactome, Hallmark, or KEGG
            resultFinal <- oraResult %>%
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

            if (verbose) {
                message(
                    "\tDone, found ", nrow(resultFinal), " enriched pathways.\n"
                )
            }
            return(resultFinal)
        }

        if (grepl(x=analysis, pattern="fgsea")) {
            if (analysis == "fgsea_reactome") {
                reactomeGeneSets <- reactomeDatabase %>%
                    mutate(
                        geneSetName = paste0(pathwayId, ";", pathwayName)
                    ) %>%
                    distinct(geneSetName, entrezGeneId) %>%
                    split(x=.$entrezGeneId, f=.$geneSetName)

                imap_dfr(
                    .x=preppedGenesTable,
                    .id="direction",
                    function(y, direction) {
                        gseaInput <- y %>%
                            tibble::as_tibble(rownames="ensemblGeneId") %>%
                            left_join(mappingFile, by="ensemblGeneId") %>%
                            mutate(
                                geneRank=-log10(PAdjusted) * LogFoldChange
                            ) %>%
                            group_by(entrezGeneId) %>%
                            summarise(geneRank = mean(geneRank)) %>%
                            ungroup() %>%
                            distinct(entrezGeneId, geneRank) %>%
                            tibble::deframe()

                        fgsea::fgsea(
                            pathways = reactomeGeneSets,
                            stats = gseaInput,
                            scoreType="pos",
                            minSize = 10,
                            maxSize = 200
                        ) %>%
                            tidyr::separate_wider_delim(
                                pathway,
                                delim=";",
                                names=c("pathwayId", "pathwayName")
                            ) %>%
                            rename(
                                "pValue"=pval,
                                "pValueAdjusted"=padj,
                            ) %>%
                            mutate(totalGenes=nrow(rnaseqResults))
                    }
                )
            } else if (analysis == "fgsea_hallmark") {
                hallmarkGeneSets <- hallmarkDatabase %>%
                    distinct() %>%
                    split(x = .$ensemblGeneId, f = .$pathwayId)

                imap_dfr(
                    .x=preppedGenesTable,
                    .id="direction",
                    function(y, direction) {
                        gseaInput <- y %>%
                            tibble::as_tibble(rownames="ensemblGeneId") %>%
                            mutate(
                                geneRank=-log10(PAdjusted) * LogFoldChange
                            ) %>%
                            group_by(ensemblGeneId) %>%
                            summarise(geneRank = mean(geneRank)) %>%
                            ungroup() %>%
                            distinct(ensemblGeneId, geneRank) %>%
                            tibble::deframe()

                        fgsea::fgsea(
                            pathways = hallmarkGeneSets,
                            stats = gseaInput,
                            scoreType="pos",
                            minSize = 10,
                            maxSize = 200
                        ) %>%
                            rename(
                                "pathwayId"=pathway,
                                "pValue"=pval,
                                "pValueAdjusted"=padj,
                            ) %>%
                            mutate(
                                pathwayName = pathwayId, .after="pathwayId"
                            ) %>%
                            mutate(totalGenes=nrow(rnaseqResults))
                    }
                )
            } else {
                return(NULL)
            }
        }
    })

    resultsAllComparisons <- resultList %>%
        bind_rows(.id="comparison") %>%
        left_join(
            select(pathwayCategories, pathwayId, topLevelPathway),
            by="pathwayId",
            multiple="all"
        ) %>%
        tibble::as_tibble()
    return(resultsAllComparisons)
}


#' INTERNAL Wrapper around Sigora's enrichment function
#'
#' @param enrichGenes Vector of genes to enrich
#' @param gpsRepo GPS object to use for testing pathways
#' @param gpsLevel Level to use for enrichment testing
#' @param pValFilter Desired threshold for filtering results
#'
#' @return A "data.frame" (tibble) of results from Sigora
#'
#' @import dplyr
#'
#' @importFrom purrr imap
#' @importFrom stringr str_count
#' @importFrom tibble tibble
#'
#' @description Internal wrapper function to run Sigora and return the results
#'   with desired columns
#'
#' @references <https://cran.r-project.org/package=sigora>
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
.runSigora <- function(
        enrichGenes,
        gpsRepo,
        gpsLevel,
        pValFilter=NA
) {

    stopifnot(
        "Input must be a character vector"={ is(enrichGenes, "character") }
    )
    stopifnot(
        "Input must be a character vector"={ is(enrichGenes, "vector") }
    )

    stopifnot(
        "Your input vector doesn't look like Ensembl genes."={
            any(grepl(pattern="^ENSG", enrichGenes))
        }
    )

    data_env <- new.env(parent=emptyenv())
    data("idmap", "reaH", "kegH", envir=data_env, package="sigora")
    idmap <- data_env[["idmap"]]
    reaH <- data_env[["reaH"]]
    kegH <- data_env[["kegH"]]

    if (gpsRepo %in% c("default", "reaH")) {
        gpsRepo <- reaH
        gpsLevel <- 4
    } else if (gpsRepo == "kegH") {
        gpsRepo <- kegH
        gpsLevel <- 2
    }

    invisible(capture.output(
        sigoraResult1 <- sigora::sigora(
            GPSrepo=gpsRepo,
            level=gpsLevel,
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
        tibble::as_tibble() %>%
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
        ) %>% tibble::as_tibble()
}
