#' pathnet_ggraph
#'
#' @param network Tidygraph network object, output from `create_pathnet`. See
#'  Details for specific requirements.
#' @param net_layout Desired layout for the network visualization. Defaults to
#'   "nicely", but supports any method found in `?layout_tbl_graph_igraph`
#' @param edge_colour Colour of network edges; defaults to "grey30".
#' @param edge_alpha Alpha value for edges; defaults to `1`.
#' @param node_size_range Size range for nodes, mapped to significance
#'   (Bonferroni p-value). Defaults to `c(4, 8)`.
#' @param edge_width_range Range of edge widths, mapped to `log10(similarity)`.
#'   Defaults to `c(0.33, 3)`.
#' @param label_prop Proportion of "interactor" (i.e. non-enriched) pathways
#'   that the function will attempt to label. E.g. setting this to 0.5 (the
#'   default) means half of the non-enriched pathways will *potentially* be
#'   labeled - it won't be exact because the node labeling is done with
#'   `ggrepel`.
#' @param node_label_size Size of node labels; defaults to 5.
#' @param node_label_alpha Transparency of node labels. Defaults to `0.67`.
#' @param node_label_overlaps Max overlaps for node labels, from `ggrepel`.
#'   Defaults to `6`.
#' @param seg_colour Colour of line segments connecting labels to nodes.
#'   Defaults to "black".
#' @param theme_base_size Base font size for all plot elements. Defaults to
#'   `16`.
#' @param seed Random seed, which will influence node positions and labels
#'
#' @return An object of class "gg"
#' @export
#'
#' @import ggplot2
#' @import ggraph
#' @import dplyr
#'
#' @description Plots the network object generated from `create_pathnet`
#'
#' @details A note regarding node labels: The function tries to prioritize
#'   labeling enriched pathways (filled nodes), with the `label_prop` argument
#'   determining roughly how many of the remaining interactor pathways might get
#'   labels. You'll likely need to tweak this value, and try different seeds, to
#'   get the desired effect.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
pathnet_ggraph <- function(
    network,
    net_layout = "nicely",
    node_size_range = c(4, 8),
    edge_colour = "grey30",
    edge_alpha = 1,
    edge_width_range = c(0.33, 3),
    label_prop = 0.25,
    node_label_size = 5,
    node_label_alpha = 0.67,
    node_label_overlaps = 6,
    seg_colour = "black",
    theme_base_size = 16,
    seed = 123
  ) {

  set.seed(seed)

  # Check column names for both nodes and edges
  stopifnot(all(
    c(
      "pathway_name_1",
      "bonferroni",
      "grouped_pathway"
    ) %in% colnames(as_tibble(network))
  ))

  stopifnot(all(
    "similarity" %in% colnames(as_tibble(tidygraph::activate(network, "edges")))
  ))

  interactors_all <- network %>%
    filter(is.na(bonferroni)) %>%
    pull(pathway_name_1)

  interactors_to_label <- sample(
    interactors_all,
    size = (length(interactors_all) * label_prop),
    replace = FALSE
  )

  network_to_plot <- network %>% mutate(
      node_fill = if_else(!is.na(bonferroni), grouped_pathway, NA_character_),
      node_label = case_when(
        !is.na(bonferroni) ~ pathway_name_1,
        pathway_name_1 %in% interactors_to_label ~ pathway_name_1,
        TRUE ~ NA_character_
      ),
      node_label = map_chr(
        node_label,
        ~trunc_neatly(.x, l = 40) %>% str_wrap(width = 20)
      ),
      bonferroni = if_else(!is.na(bonferroni), bonferroni, 1)
    )

  ggraph(network_to_plot, layout = net_layout) +
    # Edges
    geom_edge_link(
      aes(edge_width = log10(similarity)),
      colour = edge_colour,
      alpha = edge_alpha
    ) +
    scale_edge_width(range = edge_width_range, name = "Similarity") +

    # Nodes
    geom_node_point(
      aes(size = -log10(bonferroni), fill = node_fill, colour = grouped_pathway),
      pch = 21,
      stroke = 1.5
    ) +
    scale_size_continuous(
      labels = scales::label_math(expr = 10^-~.x),
      range = node_size_range
    ) +
    scale_fill_manual(
      values = top_pathway_colours,
      na.value = "white",
      guide = NULL
    ) +
    scale_colour_manual(values = top_pathway_colours) +

    # Node labels
    geom_node_label(
      aes(label = node_label),
      repel = TRUE,
      size = node_label_size,
      alpha = node_label_alpha,
      min.segment.length = 0,
      segment.colour = seg_colour,
      max.overlaps = node_label_overlaps
    ) +

    # Misc
    labs(
      size = "Bonferroni\np-value",
      colour = "Pathway type"
    ) +
    theme_void(base_size = theme_base_size) +
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
