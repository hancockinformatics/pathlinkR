#' Clean GraphML or JSON input
#'
#' @param network `tidygraph` object from a GraphML or JSON file
#'
#' @return `tidygraph` object which can be forwarded to other `pathlinkR`
#'   functions such as `ppiPlotNetwork`
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
