#' Clean GraphML or JSON input
#'
#' @param network `tidygraph` object from a GraphML or JSON file
#'
#' @return `tidygraph` object which can be forwarded to other `pathnet`
#' functions such as `ppiPlotNetwork`
#'
#' @export
#'
#' @import dplyr
#' @import tidygraph
#'
#' @seealso <https://github.com/hancockinformatics/pathnet/>
#'
#' @examples
#' tj1 <- jsonlite::read_json(
#'     system.file("extdata/networkAnalystExample.json", package = "pathnet"),
#'     simplifyVector = TRUE
#' )
#'
#' tj2 <- igraph::graph_from_data_frame(
#'     d = dplyr::select(tj1$edges, source, target),
#'     directed = FALSE,
#'     vertices = dplyr::select(
#'         tj1$nodes,
#'         id,
#'         label,
#'         x,
#'         y,
#'         "types" = molType,
#'         expr
#'     )
#' )
#'
#' tj3 <- ppiCleanNetwork(tidygraph::as_tbl_graph(tj2))
#'
ppiCleanNetwork <- function(network) {

    stopifnot(is(tidygraph::as_tbl_graph(network), "tbl_graph"))

    network %>%
        mutate(
            degree = centrality_degree(),
            betweenness = centrality_betweenness(),
            seed = if_else(types == "Seed", TRUE, FALSE),
            hubScoreBtw = centrality_betweenness()
        ) %>%
        select(
            name,
            degree,
            betweenness,
            seed,
            hubScoreBtw,
            "geneName" = label,
            everything()
        ) %>%
        ppiRemoveSubnetworks() %>%
        as_tbl_graph()
}
