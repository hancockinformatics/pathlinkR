#' Create a pathway network from enrichment results and a pathway
#' interaction foundation
#'
#' @param pathwayEnrichmentResult Data frame of results from Sigora or
#'   ReactomePA (should be based on Reactome data)
#' @param columnId Character; column containing the Reactome pathway IDs.
#'   Defaults to "pathwayID".
#' @param columnP Character; column containing the adjusted p values. Defaults
#'   to "pValueAdjusted".
#' @param foundation List of pathway pairs to use in constructing a network.
#'   Typically this will be the output from `createFoundation`.
#' @param trim Remove independent subgraphs which don't contain any enriched
#'   pathways (default is `TRUE`).
#' @param trimOrder Order to use when removing subgraphs; Higher values will
#'   keep more non-enriched pathway nodes. Defaults to `1`.
#'
#' @return A pathway network as a "tidygraph" object, with the following columns
#'  for nodes:
#'  \item{pathwayId}{Reactome pathway ID}
#'  \item{pathwayName}{Reactome pathway name}
#'  \item{comparison}{Name of source comparison, if this pathway was enriched}
#'  \item{direction}{Whether an enriched pathway was found in all genes or up-
#'  or down-regulated genes}
#'  \item{pValue}{Nominal p-value from the enrichment result}
#'  \item{pValueAdjusted}{Corrected p-value from the enrichment}
#'  \item{genes}{Candidate genes for the given pathway if it was enriched}
#'  \item{numCandidateGenes}{Number of candidate genes}
#'  \item{numBgGenes}{Number of background genes}
#'  \item{geneRatio}{Ratio of candidate and background genes}
#'  \item{totalGenes}{Total number of DE genes tested, for an enriched pathway}
#'  \item{topLevelPathway}{Highest level Reactome term for a given pathway}
#'  \item{groupedPathway}{Custom pathway category used in visualizations}
#'
#'  For edges, the following information is also included:
#'  \item{from}{Starting node (row number) for the edge}
#'  \item{to}{Ending node (row number) for the edge}
#'  \item{similarity}{Similarity of two nodes/pathways}
#'  \item{distance}{Inverse of similarity}
#'
#' @export
#'
#' @import dplyr
#'
#' @importFrom igraph as.igraph neighborhood
#' @importFrom purrr map
#' @importFrom tidygraph activate tbl_graph
#'
#' @description Creates a tidygraph network object from the provided pathway
#'   information, ready to be visualized with `pathnetGGraph` or
#'   `pathnetVisNetwork`.
#'
#' @details With the "trim" option enabled, nodes (pathways) and subgraphs which
#'   are not sufficiently connected to enriched pathways will be removed. How
#'   aggressively this is done can be controlled via the `trimOrder` argument,
#'   and the optimal value will depend on the number of enriched pathways and
#'   the number of interacting pathways (i.e. number of rows in "foundation").
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' data("sigoraDatabase", "sigoraExamples")
#'
#' pathwayDistancesJaccard <- getPathwayDistances(
#'     pathwayData=dplyr::slice_head(
#'         dplyr::arrange(sigoraDatabase, pathwayId),
#'         prop=0.1
#'     ),
#'     distMethod="jaccard"
#' )
#'
#' startingPathways <- pathnetFoundation(
#'     mat=pathwayDistancesJaccard,
#'     maxDistance=0.8
#' )
#'
#' pathnetCreate(
#'     pathwayEnrichmentResult=sigoraExamples[grepl(
#'         "Pos",
#'         sigoraExamples$comparison
#'     ), ],
#'     foundation=startingPathways,
#'     trim=TRUE,
#'     trimOrder=1
#' )
#'
pathnetCreate <- function(
        pathwayEnrichmentResult,
        columnId="pathwayId",
        columnP="pValueAdjusted",
        foundation,
        trim=TRUE,
        trimOrder=1
) {

    stopifnot(is(pathwayEnrichmentResult, "data.frame"))

    stopifnot(all(
        c(columnId, columnP) %in% colnames(pathwayEnrichmentResult)
    ))

    stopifnot(is(foundation, "data.frame"))
    stopifnot(all(
        c("pathwayName1", "pathwayName2", "distance", "pathway1", "pathway2")
        %in% colnames(foundation)
    ))

    data_env <- new.env(parent = emptyenv())
    data("pathwayCategories", envir = data_env, package = "pathlinkR")
    pathwayCategories <- data_env[["pathwayCategories"]]

    if (columnId != "pathwayId") {
        pathwayEnrichmentResult <- pathwayEnrichmentResult %>%
            rename("pathwayId"=all_of(columnId))
    }

    if (columnP != "pValueAdjusted") {
        pathwayEnrichmentResult <- pathwayEnrichmentResult %>%
            rename("pValueAdjusted"=all_of(columnP))
    }

    startingNodes <- foundation %>%
        distinct(pathway1, pathwayName1) %>%
        left_join(
            pathwayEnrichmentResult,
            by=c("pathway1" = "pathwayId"),
            multiple="all"
        )

    startingEdges <- foundation %>%
        mutate(similarity = 1 / distance) %>%
        select(pathway1, pathway2, similarity, distance)

    pathwaysAsNetwork <- tbl_graph(
        nodes=startingNodes,
        edges=startingEdges,
        directed=FALSE
    ) %>% mutate(rn=row_number())

    pathwaysAsNetwork2 <- if (trim) {
        x1 <- pathwaysAsNetwork %>%
            filter(!is.na(pValueAdjusted)) %>%
            pull(rn)

        validNodes <- map(x1, ~neighborhood(
            graph=as.igraph(pathwaysAsNetwork),
            order=trimOrder,
            nodes=.x
        )) %>%
            unlist() %>%
            unique()

        pathwaysAsNetwork %>%
            filter(rn %in% c(x1, validNodes)) %>%
            select(-rn)
    } else {
        select(pathwaysAsNetwork, -rn)
    }

    pathwaysAsNetwork3 <- pathwaysAsNetwork2 %>%
        activate("edges") %>%
        distinct() %>%
        activate("nodes")

    pathwaysAsNetwork4 <- pathwaysAsNetwork3 %>%
        left_join(
            y =select(
                pathwayCategories, pathwayId,
                pathwayName, groupedPathway
            ),
            by=c("pathway1" = "pathwayId"),
            multiple="all"
        ) %>%
        select(
            "pathwayId"=pathway1,
            "pathwayName"=pathwayName1,
            everything(),
            -any_of(c(
                "pathwayName.x", "pathwayName.y", "level_1",
                "level_2", "level1", "level2"
            ))
        )

    return(pathwaysAsNetwork4)
}
