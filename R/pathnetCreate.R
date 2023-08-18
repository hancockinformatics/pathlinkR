#' Create a pathway network from enrichment results and a pathway
#' interaction foundation
#'
#' @param pathwayEnrichmentResult Data frame of results from Sigora or
#'   ReactomePA (should be based on Reactome data). Must minimally contain the
#'   columns "pathwayId" and "pValueAdjusted".
#' @param foundation List of pathway pairs to use in constructing a network.
#'   Typically this will be the output from `createFoundation`.
#' @param trim Remove independent subgraphs which don't contain any enriched
#'   pathways (default is `TRUE`).
#' @param trimOrder Order to use when removing subgraphs; Higher values will
#'   keep more non-enriched pathway nodes. Defaults to `1`.
#'
#' @return A pathway network as a tidygraph object
#' @export
#'
#' @import dplyr
#' @import purrr
#' @import stringr
#' @import tidygraph
#' @importFrom igraph neighborhood as.igraph
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
#' startingPathways <- pathnetFoundation(
#'     mat=pathwayDistancesJaccard,
#'     maxDistance=0.8
#' )
#'
#' pathnetCreate(
#'     pathwayEnrichmentResult=sigoraExamples[
#'         grepl("Pos", sigoraExamples$comparison),
#'     ],
#'     foundation=startingPathways,
#'     trim=TRUE,
#'     trimOrder=1
#' )
#'
pathnetCreate <- function(
        pathwayEnrichmentResult,
        foundation,
        trim=TRUE,
        trimOrder=1
) {

    stopifnot(is(pathwayEnrichmentResult, "data.frame"))
    stopifnot(all(
        c("pathwayId", "pValueAdjusted") %in% colnames(pathwayEnrichmentResult)
    ))

    stopifnot(is(foundation, "data.frame"))
    stopifnot(all(
        c("pathwayName1", "pathwayName2", "distance", "pathway1", "pathway2")
        %in% colnames(foundation)
    ))

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
                topPathwaysMore, pathwayId,
                pathwayName, groupedPathway
            ),
            by=c("pathway1" = "pathwayId"),
            multiple="all"
        ) %>%
        select(
            pathway1,
            pathwayName1,
            everything(),
            -any_of(c(
                "pathwayName.x", "pathwayName.y", "level_1",
                "level_2", "level1", "level2"
            ))
        )

    return(pathwaysAsNetwork4)
}
