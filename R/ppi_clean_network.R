#' Clean GraphML or JSON input
#'
#' @param network `tidygraph` object from a GraphML or JSON file
#'
#' @return `tidygraph` object which can be forwarded to other `networker`
#' functions such as `plot_network`
#'
#' @export
#'
#' @import dplyr
#' @import tidygraph
#'
#' @seealso <https://github.com/hancockinformatics/pathnet/>
#'
ppi_clean_network <- function(network) {
    network %>%
        mutate(
            degree = centrality_degree(),
            betweenness = centrality_betweenness(),
            seed = if_else(types == "Seed", TRUE, FALSE),
            hub_score_btw = centrality_betweenness()
        ) %>%
        dplyr::select(
            name,
            degree,
            betweenness,
            seed,
            hub_score_btw,
            "gene_name" = label,
            everything()
        ) %>%
        ppi_remove_subnetworks() %>%
        as_tbl_graph()
}
