#' createPathnet
#'
#' @param sigoraResult Data frame of results from Sigora. Must contain columns
#' "pathwayId" and "pValueAdjusted".
#' @param foundation List of pathway pairs to use in constructing a network,
#'   output from `createFoundation`.
#' @param trim Remove subgraphs which don't contain any enriched pathways
#'   (default is `TRUE`).
#' @param trimOrder Order to use when removing subgraphs; Higher values will
#'   keep more non-enriched pathway nodes. Defaults to `1`.
#'
#' @return A tidygraph network object
#' @export
#'
#' @import dplyr
#' @import stringr
#' @importFrom igraph neighborhood as.igraph
#' @importFrom tidygraph tbl_graph activate
#' @importFrom purrr map map_chr
#'
#' @description Creates a tidygraph network object from the pathway information,
#'   ready to be visualized with `pathnetGGraph` or `pathnetVisNetwork`.
#'
#' @details With the "trim" option enabled, nodes (pathways) and subgraphs which
#'   are not sufficiently connected to enriched pathways will be removed. How
#'   aggressively this is done can be controlled via the `trimOrder` argument.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' startingPathways <- createFoundation(
#'     mat = pathwayDistancesJaccard,
#'     max_distance = 0.8
#' )
#'
#' exPathnet <- createPathnet(
#'     sigoraResult = sigoraExamples,
#'     foundation = startingPathways,
#'     trim = TRUE,
#'     trimOrder = 1
#' )
#'
createPathnet <- function(
        sigoraResult,
        foundation,
        trim = TRUE,
        trimOrder = 1
) {

    stopifnot(
        all(c("pathwayId", "pValueAdjusted") %in% colnames(sigoraResult))
    )

    startingNodes <- foundation %>%
        select(pathway1, pathwayName1) %>%
        distinct() %>%
        left_join(
            sigoraResult,
            by = c("pathway1" = "pathwayId"),
            multiple = "all"
        )

    startingEdges <- foundation %>%
        mutate(similarity = 1 / distance) %>%
        select(pathway1, pathway2, similarity, distance)

    pathwaysAsNetwork <- tbl_graph(
        nodes = startingNodes,
        edges = startingEdges,
        directed = FALSE
    ) %>%
        mutate(rn = row_number())

    pathwaysAsNetwork2 <- if (trim) {
        x1 <- pathwaysAsNetwork %>%
            filter(!is.na(pValueAdjusted)) %>%
            pull(rn)

        validNodes <- map(x1, ~neighborhood(
            graph = as.igraph(pathwaysAsNetwork),
            order = trimOrder,
            nodes = .x
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
            y  = select(
                topPathwaysMore, pathwayId,
                pathwayName, groupedPathway
            ),
            by = c("pathway1" = "pathwayId"),
            multiple = "all"
        ) %>%
        select(
            pathway1,
            pathwayName1,
            everything(),
            -any_of(c("pathwayName", "level_1", "level_2", "level1", "level2"))
        )

    return(pathwaysAsNetwork4)
}
