#' Clean GraphML or JSON input
#'
#' @param network `tidygraph` object from a GraphML or JSON file
#'
#' @return A Protein-Protein Interaction (PPI) network; a "tidygraph" object,
#'  with the minimal set of columns (other from the input are also included):
#'   \item{name}{Identifier for the node}
#'   \item{degree}{Degree of the node, i.e. the number of interactions}
#'   \item{betweenness}{Betweenness measure for the node}
#'   \item{seed}{TRUE when the node was part of the input list of genes}
#'   \item{hubScore}{Special hubScore for each node. The suffix denotes the
#'   measure being used; e.g. "hubScoreBtw" is for betweenness}
#'   \item{hgncSymbol}{HGNC gene name for the node}
#'
#' Additionally the following columns are provided for edges:
#'   \item{from}{Starting node for the interaction/edge as a row number}
#'   \item{to}{Ending node for the interaction/edge as a row number}
#'
#' @export
#'
#' @import dplyr
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph centrality_betweenness centrality_degree
#'
#' @description Takes network file (GraphML or JSON) and process it into a
#'   tidygraph object, adding network statistics along the way.
#'
#' @details This function was designed so that networks created by other
#'   packages or websites (e.g. <https://networkanalyst.ca>) could be imported
#'   and visualized with `ppiPlotNetwork`.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR/>
#'
#' @examples
#' tj1 <- jsonlite::read_json(
#'     system.file("extdata/networkAnalystExample.json", package="pathlinkR"),
#'     simplifyVector=TRUE
#' )
#'
#' tj2 <- igraph::graph_from_data_frame(
#'     d=dplyr::select(tj1$edges, source, target),
#'     directed=FALSE,
#'     vertices=dplyr::select(
#'         tj1$nodes,
#'         id,
#'         label,
#'         x,
#'         y,
#'         "types"=molType,
#'         expr
#'     )
#' )
#'
#' tj3 <- ppiCleanNetwork(tidygraph::as_tbl_graph(tj2))
#'
ppiCleanNetwork <- function(network) {

    stopifnot(is(as_tbl_graph(network), "tbl_graph"))

    network %>%
        mutate(
            degree=centrality_degree(),
            betweenness=centrality_betweenness(),
            seed=if_else(types == "Seed", TRUE, FALSE),
            hubScoreBtw=centrality_betweenness()
        ) %>%
        select(
            name,
            degree,
            betweenness,
            seed,
            hubScoreBtw,
            "hgncSymbol"=label,
            everything()
        ) %>%
        ppiRemoveSubnetworks() %>%
        as_tbl_graph()
}
