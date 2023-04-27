#' plot_pathnet_plotly
#'
#' @param network Network object produced by `create_pathnet`
#'
#' @return Plotly object for interactive visualization
#' @export
#'
#' @importFrom plotly plot_ly add_segments add_trace layout
#' @importFrom igraph as.igraph layout.auto get.edgelist get.vertex.attribute
#' @import dplyr
#' @import purrr
#'
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
plot_pathnet_plotly <- function(network) {

  network_igraph <- as.igraph(network)
  network_layout <- layout.auto(network_igraph)

  plotly_node_info <- bind_cols(
    as.data.frame(get.vertex.attribute(network_igraph)),
    as.data.frame(network_layout)
  ) %>%
    relocate(pathway_name_1, pathway_1, "xn" = V1, "yn" = V2, everything())

  plotly_edge_info <- as.data.frame(get.edgelist(network_igraph))

  network_edges <- list(); for (i in 1:nrow(plotly_edge_info)) {
    i_edge_shape <- list(
      x0 = plotly_node_info[plotly_edge_info$V1[i], "xn"],
      y0 = plotly_node_info[plotly_edge_info$V1[i], "yn"],
      x1 = plotly_node_info[plotly_edge_info$V2[i], "xn"],
      y1 = plotly_node_info[plotly_edge_info$V2[i], "yn"]
    )
    network_edges[[i]] <- i_edge_shape
  }
  network_edges_df <- map(network_edges, as.data.frame) %>% bind_rows()

  plotly_axis <- list(
    title = "",
    showgrid = FALSE,
    showticklabels = FALSE,
    zeroline = FALSE
  )

  p <- plot_ly() %>%
    add_segments(
      inherit = FALSE,
      type = "scatter",
      mode = "lines",
      data = network_edges_df,
      x = ~x0,
      y = ~y0,
      xend = ~x1,
      yend = ~y1,
      line = list(color = "black", size = 2),
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "scatter",
      mode = "markers",
      data = plotly_node_info,
      x = ~xn,
      y = ~yn,
      color = ~grouped_pathway,
      marker = list(
        size = 20,
        line = list(color = "grey30", width = 1),
        colorscale = top_pathway_colours
      ),
      text = ~pathway_name_1,
      hoverinfo = "text",
      inherit = FALSE
    ) %>%
    layout(xaxis = plotly_axis, yaxis = plotly_axis)

  return(p)
}
