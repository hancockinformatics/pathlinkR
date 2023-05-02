#' pathnet_visNetwork
#'
#' @param network Tidygraph network object as output by `create_pathnet`
#' @param net_layout Desired layout for the network visualization. Defaults to
#'   "layout_nicely", should support most igraph layouts. See
#'   `?visIgraphLayout()` for more details.
#' @param edge_colour Colour of network edges; defaults to "#848484".
#' @param edge_width Edge width, defaults to 6.
#' @param node_size Size of the nodes, defaults to 35.
#' @param set_seed Random seed to use for reproducible node placement, defaults
#'   to 123
#'
#' @return Interactive visNetwork plot
#' @export
#'
#' @import dplyr
#' @import visNetwork
#' @importFrom igraph as.igraph
#' @importFrom purrr map_chr map2
#'
#'
#' @references <https://datastorm-open.github.io/visNetwork/>
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
pathnet_visNetwork <- function(
    network,
    net_layout = "layout_nicely",
    edge_colour = "#848484",
    edge_width = 6,
    node_size = 35,
    set_seed = 123
  ) {

  my_igraph <- network %>%
    mutate(
      title = pathway_name_1,
      background = map_chr(grouped_pathway, ~top_pathway_colours[[.x]]),
      background = if_else(!is.na(bonferroni), background, "#ffffff"),
      border = map_chr(grouped_pathway, ~top_pathway_colours[[.x]]),
      color = map2(background, border, ~list("background" = .x, "border" = .y))
    ) %>%
    select(-c(background, border)) %>%
    as.igraph()

  legend_df <- top_pathway_colours %>%
    enframe("label", "icon.color") %>%
    mutate(shape = "dot", size = 15)

  visIgraph(my_igraph) %>%
    visIgraphLayout(layout = net_layout, randomSeed = set_seed) %>%
    visEdges(color = edge_colour, width = edge_width) %>%
    visNodes(size = node_size, borderWidth = 3) %>%
    visOptions(highlightNearest = TRUE, selectedBy = list(
      "variable" = "grouped_pathway",
      "style" = "width: 175px; height: 26px",
      "main" = "Select pathway group:"
    )) %>%
    visLegend(
      useGroups = FALSE,
      position = "right",
      main = "Grouped pathway",
      addNodes = legend_df
    ) %>%
    visExport(float = "right")
}
