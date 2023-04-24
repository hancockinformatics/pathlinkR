#' plot_pathnet
#'
#' @param network Tidygraph network object, output from `create_pathnet`. See
#'  Details for specific requirements.
#' @param net_layout Desired layout for the network visualization. Defaults to
#'   "nicely", but supports any method found in `?layout_tbl_graph_igraph`
#' @param edge_alpha Alpha value for edges
#'
#' @return An object of class "gg"
#' @export
#'
#' @import ggplot2
#' @import ggraph
#' @import dplyr
#'
#' @description Plots the network object
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
plot_pathnet <- function(network, net_layout = "nicely", edge_alpha = 0.67) {

  # Check column names for both nodes and edges
  stopifnot(all(
    c(
      "pathway_name_1",
      "description",
      "grouped_pathway",
      "bonferroni"
    ) %in% colnames(as_tibble(network))
  ))

  stopifnot(all(
    "similarity" %in% colnames(as_tibble(tidygraph::activate(network, "edges")))
  ))

  network <- network %>%
    mutate(
      node_fill = if_else(!is.na(description), grouped_pathway, NA_character_),
      pathway_name_1_wrap = map_chr(
        as.character(pathway_name_1),
        ~str_wrap(.x, width = 20)
      ),
      pathway_name_1_wrap = str_replace(pathway_name_1_wrap, "^$", NA_character_)
    )

  ggraph(network, layout = net_layout) +
    geom_edge_link(aes(edge_width = log10(similarity)), alpha = edge_alpha) +
    geom_node_point(
      aes(size = -log10(bonferroni), fill = node_fill, colour = grouped_pathway),
      pch = 21,
      stroke = 1.5
    ) +
    geom_node_label(
      aes(label = pathway_name_1_wrap),
      repel = TRUE,
      alpha = 0.67,
      min.segment.length = 0,
      max.overlaps = 6
    ) +
    scale_edge_width(range = c(0.33, 3), name = "Similarity") +
    scale_size_continuous(
      labels = scales::label_math(expr = 10^-~.x),
      range = c(4, 8)
    ) +
    scale_colour_manual(values = top_pathway_colours) +
    scale_fill_manual(
      values = top_pathway_colours,
      na.value = "white",
      guide = NULL
    ) +
    scale_shape_manual(values = c("up" = 21), na.value = 16) +
    labs(
      size = "Bonferroni\np-value",
      colour = "Pathway type"
    ) +
    theme_void(base_size = 18) +
    theme(
      legend.text.align = 0,
      plot.margin = unit(rep(5, 4), "mm")
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = 5, pch = 19)),
      size = guide_legend(override.aes = list(
        colour = "black",
        fill = "white",
        stroke = 0.5
      ))
    )
}
