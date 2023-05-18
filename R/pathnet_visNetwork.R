#' Visualize enriched Reactome pathways as an interactive network
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
#' @param label_nodes Boolean determining if nodes should be labeled. Note it
#'   will only ever label enriched nodes/pathways.
#' @param node_label_size Size of the node labels in pixels; defaults to 60.
#' @param node_label_colour Colour of the node labels; defaults to "black".
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
#' @description Creates a pathway network using the visNetwork library, which
#'   allows for various forms of interactivity such as including text when
#'   hovering over nodes, node selection and dragging (including multiple
#'   selections).
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
    label_nodes = TRUE,
    node_label_size = 60,
    node_label_colour = "black",
    highlighting = TRUE,
    set_seed = 123
) {

  visnet_nodes <- network %>%
    as_tibble() %>%
    mutate(
      id = row_number(),
      value = if_else(is.na(bonferroni), 1, -log10(bonferroni)),
      background = map_chr(grouped_pathway, ~grouped_pathway_colours[[.x]]),
      background = if_else(!is.na(bonferroni), background, "#ffffff"),
      border = map_chr(grouped_pathway, ~grouped_pathway_colours[[.x]]),
      color = map2(background, border, ~list("background" = .x, "border" = .y))
    ) %>%
    dplyr::select(
      id,
      "title" = pathway_name_1,
      everything(),
      -c(background, border, direction, description, pvalue)
    )

  if (label_nodes) {
    visnet_nodes <- mutate(
      visnet_nodes,
      label = map_chr(if_else(!is.na(bonferroni), title, ""), trunc_neatly, 30)
    )
  }

  visnet_edges <- network %>%
    activate("edges") %>%
    as_tibble() %>%
    rename("value" = similarity) %>%
    distinct()

  legend_df <- grouped_pathway_colours %>%
    enframe("label", "icon.color") %>%
    mutate(shape = "dot", size = 15)

  out1 <- visNetwork(nodes = visnet_nodes, edges = visnet_edges) %>%
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

  if (label_nodes) {
    out1 %>%
      visNodes(font = paste0(node_label_size, "px arial ", node_label_colour))
  } else {
    out1
  }
}
