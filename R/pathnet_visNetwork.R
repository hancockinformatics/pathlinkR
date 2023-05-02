#' pathnet_visNetwork
#'
#' @param network Tidygraph network object as output by `create_pathnet`
#' @param net_layout Desired layout for the network visualization. Defaults to
#'   "layout_nicely", should support most igraph layouts. See
#'   `?visIgraphLayout()` for more details.
#' @param edge_colour Colour of network edges; defaults to "#848484".
#' @param edge_size_range Edge width is mapped to the similarity measure (one
#'   over distance). This length-two numeric vector controls the minimum and
#'   maximum width of edges. Defaults to `c(5, 20)`.
#' @param node_size_range Node size is mapped to the negative log of the
#'   Bonferroni-adjusted p value, and this length-two numeric vector controls
#'   the minimum and maximum. Defaults to `c(20, 50)`.
#' @param node_border_width Size of the node border, defaults to 2.5
#' @param highlighting When clicking on a node, should directly neighbouring
#'   nodes be highlighted (other nodes are dimmed)? Defaults to TRUE.
#' @param set_seed Random seed to use for reproducible node placement, defaults
#'   to 123.
#'
#' @return Interactive visNetwork plot
#' @export
#'
#' @import dplyr
#' @import visNetwork
#' @importFrom igraph as.igraph
#' @importFrom purrr map_chr map2
#' @importFrom tidygraph activate
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
    edge_size_range = c(5, 20),
    node_size_range = c(20, 50),
    node_border_width = 2.5,
    highlighting = TRUE,
    set_seed = 123
) {

  my_igraph <- network %>%
    mutate(
      title = pathway_name_1,
      value = if_else(is.na(bonferroni), 1, -log10(bonferroni)),
      background = map_chr(grouped_pathway, ~top_pathway_colours[[.x]]),
      background = if_else(!is.na(bonferroni), background, "#ffffff"),
      border = map_chr(grouped_pathway, ~top_pathway_colours[[.x]]),
      color = map2(background, border, ~list("background" = .x, "border" = .y))
    ) %>%
    select(-c(background, border)) %>%
    activate("edges") %>%
    rename("value" = similarity) %>%
    activate("nodes") %>%
    as.igraph()

  legend_df <- top_pathway_colours %>%
    enframe("label", "icon.color") %>%
    mutate(shape = "dot", size = 15)

  visIgraph(my_igraph) %>%
    visIgraphLayout(layout = net_layout, randomSeed = set_seed) %>%
    visEdges(
      color = edge_colour,
      scaling = list("min" = edge_size_range[1], "max" = edge_size_range[2])
    ) %>%
    visNodes(
      borderWidth = node_border_width,
      scaling = list("min" = node_size_range[1], "max" = node_size_range[2])
    ) %>%
    visOptions(
      highlightNearest = highlighting,
      selectedBy = list(
        "variable" = "grouped_pathway",
        "style" = "width: 175px; height: 26px",
        "main" = "Select pathway group"
      )
    ) %>%
    visLegend(
      stepY = 75,
      useGroups = FALSE,
      position = "right",
      main = "Grouped pathway",
      addNodes = legend_df,
    ) %>%
    visExport(float = "right")
}
