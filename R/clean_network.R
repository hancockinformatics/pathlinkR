#' Clean graphml input
#'
#' @param network `tidygraph` object from graphml file
#'
#' @return `tidygraph` object which can be forwarded to other `networker`
#' functions such as `plot_network`
#'
#' @export
#'
#' @import dplyr
#' @import tidygraph
#'
#' @seealso <https://github.com/travis-m-blimkie/networker>
#'
clean_network <- function(network) {
  network %>%
    janitor::clean_names() %>%
    mutate(
      degree = centrality_degree(),
      betweenness = centrality_betweenness(),
      seed = if_else(types == "Seed", TRUE, FALSE),
      hub_score_btw = centrality_betweenness()
    ) %>%
    select(
      name,
      degree,
      betweenness,
      seed,
      hub_score_btw,
      "gene_name" = label,
      everything()
    ) %>%
    remove_subnetworks() %>%
    as_tbl_graph()
}
