#' Plot an undirected PPI network using ggraph
#'
#' @param network `tidygraph` object, output from `build_network`
#' @param fill_column Tidy-select column for mapping node colour. Designed to
#'   handle continuous numeric mappings (either positive/negative only, or
#'   both), and categorical mappings, plus a special case for displaying fold
#'   changes from, for example, RNA-Seq data. See `fill_type` for more details
#'   on how to set this up.
#' @param fill_type String denoting type of fill mapping to perform for nodes.
#'   Options are: "fold_change", "two_sided", "one_sided", or "categorical".
#' @param cat_fill_colours Colour palette to be used when `fill_type` is set
#'   to "categorical." Defaults to "Set1" from RColorBrewer. Will otherwise be
#'   passed as the "values" argument in `scale_fill_manual()`.
#' @param layout Layout of nodes in the network. Supports all layouts from
#'   `ggraph`/`igraph`, or a data frame of x and y coordinates for each
#'   node (order matters!).
#' @param legend Should a legend be included? Defaults to FALSE.
#' @param fontfamily Font to use for labels and legend (if present). Defaults to
#'   "Helvetica".
#' @param edge_colour_ Edge colour, defaults to "grey40"
#' @param edge_alpha_ Transparency of edges, defaults to 0.5
#' @param edge_width_ Thickness of edges connecting nodes. Defaults to 0.5
#' @param node_size Numeric vector of length two, specifying size range of node
#'   sizes (maps to node degree). Default is `c(3, 9)`.
#' @param node_colour Colour (stroke or outline) of all nodes in the network.
#'   Defaults to "grey30".
#' @param int_colour Fill colour for non-seed nodes, i.e. interactors. Defaults
#'   to "grey70".
#' @param fc_up_col Colour to use for up regulated nodes when `fill_type` is
#'   set to "fold_change". Defaults to "firebrick3".
#' @param fc_down_col Colour to use for down regulated nodes when `fill_type`
#'   is set to "fold_change". Defaults to "#188119".
#' @param label Boolean, whether labels should be added to nodes. Defaults to
#'   FALSE.
#' @param label_column Tidy-select column of the network/data to be used in
#'   labeling nodes. Recommend setting to `gene_name`, which contains HGNC
#'   symbols mapped from the input Ensembl IDs via biomaRt.
#' @param label_filter Degree filter used to determine which nodes should be
#'   labeled. Defaults to 0. This value can be increased to reduce the number of
#'   node labels, to prevent the network from being too crowded.
#' @param label_size Size of node labels, defaults to 5.
#' @param label_colour Colour of node labels, defaults to "black"
#' @param hub_colour Colour of node labels for hubs. The top 2% of nodes
#'   (based on calculated hub score) are highlighted with this colour, if
#'   `label = TRUE`.
#' @param label_face Font face for node labels, defaults to "bold"
#' @param label_padding Padding around the label, defaults to 0.25 lines.
#' @param min_seg_length Minimum length of lines to be drawn from labels to
#'   points. The default specified here is 0.25, half of the normal default
#'   value.
#' @param subnet Logical determining if networks produced by
#'   `extract_subnetwork` should be treated as such, or just as a normal network
#'   from `build_network`.
#' @param ... Further parameters can be passed on to `ggplot2::theme()`, e.g.
#'   `legend.position`
#'
#' @return An object of class "gg"
#'
#' @export
#'
#' @import ggplot2
#' @import ggraph
#' @import dplyr
#'
#' @details Any layout supported by ggraph can be specified here - see
#'   `?layout_tbl_graph_igraph` for a list of options. Or you can supply a data
#'   frame containing coordinates for each node. The first and second columns
#'   will be used for x and y, respectively. Note that having columns named "x"
#'   and "y" in the input network will generate a warning message when supplying
#'   custom coordinates.
#'
#'   Since this function returns a standard ggplot object, you can tweak the
#'   final appearance using the normal array of ggplot2 function, e.g. `labs()`
#'   and `theme()` to further customize the final appearance.
#'
#'   The `fill_type` argument will determine how the node colour is mapped to
#'   the desired column. "fold_change" represents a special case, where the fill
#'   column is numeric and whose values should be mapped to up (> 0) or down (<
#'   0). "two_sided" and "one_sided" are designed for numeric data that contains
#'   either positive and negative values, or only positive/negative values,
#'   respectively. "categorical" handles any other non-numeric colour mapping,
#'   and uses "Set1" from RColorBrewer.
#'
#'   Node statistics (degree, betweenness, and hub score) are calculated using
#'   the respective functions from the `tidygraph` package.
#'
#'   If plotting a network created from the `extract_subnetwork` function, the
#'   genes belonging to the extracted pathway (i.e. the contents of the
#'   "gene_id" column in the enrichment results) will be highlighted, instead of
#'   hub nodes. The colour used here can be controlled via the `hub_colour`
#'   argument, and this behaviour can be turned off altogether by setting the
#'   `subnet` argument to FALSE to return an "standard" network. Additionally,
#'   one can set the "starters" attribute of the network being plotted to a
#'   custom character vector of Ensembl gene IDs, to allow highlighting a custom
#'   selection of nodes.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet/>
#'
#' @examples
#' library(dplyr)
#'
#' ex_de_genes <- deseq_example_list[[2]] %>%
#'     tibble::rownames_to_column("gene") %>%
#'     as_tibble() %>%
#'     filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
#'
#' ex_network <- ppi_build_network(
#'     df = ex_de_genes,
#'     col = "gene",
#'     order = "zero"
#' )
#'
#' ppi_plot_network(
#'     ex_network,
#'     fill_column = log2FoldChange,
#'     fill_type = "fold_change",
#'     layout = "lgl",
#'     label = TRUE,
#'     label_column = gene_name,
#'     label_filter = 5,
#'     legend = TRUE
#' )
#'
ppi_plot_network <- function(
        network,
        fill_column,
        fill_type,
        cat_fill_colours = "Set1",
        layout           = "kk",
        legend           = FALSE,
        fontfamily       = "Helvetica",
        edge_colour_     = "grey40",
        edge_alpha_      = 0.5,
        edge_width_      = 0.5,
        node_size        = c(3, 9),
        node_colour      = "grey30",
        int_colour       = "grey70",
        fc_up_col        = "firebrick3",
        fc_down_col      = "#188119",
        label            = FALSE,
        label_column,
        label_filter     = 0,
        label_size       = 5,
        label_colour     = "black",
        hub_colour       = "blue2",
        label_face       = "bold",
        label_padding    = 0.25,
        min_seg_length   = 0.25,
        subnet           = TRUE,
        ...
) {

    ## Set up fill scaling based on argument `fill_type`
    if (fill_type == "fold_change") {
        stopifnot(is.numeric(pull(network, {{fill_column}})))

        network <- network %>% mutate(
            new_fill_col = case_when(
                {{fill_column}} < 0 ~ "Down",
                {{fill_column}} > 0 ~ "Up",
                TRUE ~ NA_character_
            )
        )

        network_fill_geom <- scale_fill_manual(
            values   = c("Up" = fc_up_col, "Down" = fc_down_col),
            na.value = int_colour
        )

        network_fill_guide <- guides(
            fill = guide_legend(
                title = "Direction",
                override.aes = list(size = 5)
            )
        )

    } else if (fill_type == "two_sided") {

        network <- mutate(network, new_fill_col = {{fill_column}})
        network_fill_geom <- scale_fill_gradient2(
            low  = "#313695",
            mid  = "white",
            high = "#a50026",
            midpoint = 0,
            na.value = int_colour,
            guide    = ifelse(legend, "colourbar", "none")
        )
        network_fill_guide <- NULL

    } else if (fill_type == "one_sided") {

        network <- mutate(network, new_fill_col = {{fill_column}})
        network_fill_geom <- scale_fill_viridis_c(
            option = "plasma", begin = 0.2
        )
        network_fill_guide <- NULL

    } else if (fill_type == "categorical") {

        network <- mutate(network, new_fill_col = {{fill_column}})

        if (all(cat_fill_colours == "Set1")) {
            network_fill_geom <- scale_fill_brewer(
                palette  = "Set1",
                na.value = int_colour,
                guide    = ifelse(legend, "legend", "none")
            )
        } else {
            network_fill_geom <- scale_fill_manual(values = cat_fill_colours)
        }
        network_fill_guide <-
            guides(fill = guide_legend(override.aes = list(size = 5)))

    } else {
        stop(
            "Argument 'fill_type' must be one of 'fold_change', 'two_sided', ",
            "'one_sided', or 'categorical'"
        )
    }

    if (is.data.frame(layout)) {
        message("Using user-supplied node coordinates...")
        ## By converting the layout object to a matrix, we no longer need to
        ## worry about column names. The first and second column will be "x" and
        ## "y", respectively.
        layout_object <- as.matrix(layout)
    } else {
        layout_object <- layout
    }

    ## Set a plain white background
    set_graph_style(foreground = "white")

    ## Get the fill_column as a string, so we can clean it up if the legend is
    ## included
    legend_name <- match.call()$fill_column

    ## Theme tweaks for all plot types
    theme_tweaks <- theme(
        text         = element_text(family = fontfamily),
        plot.margin  = unit(rep(0, 4), "cm"),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14),
        ...
    )

    if (!label) {
        ggraph(network, layout = layout_object) +
            geom_edge_link(
                show.legend = FALSE,
                edge_alpha = edge_alpha_,
                edge_colour = edge_colour_,
                edge_width = edge_width_
            ) +
            geom_node_point(
                aes(size = degree, fill = new_fill_col),
                pch = 21,
                colour = node_colour
            ) +
            network_fill_geom +
            scale_size_continuous(range = node_size, guide = "none") +
            labs(fill = NULL) +
            theme_tweaks +
            network_fill_guide

    } else {

        hub_nodes <- as_tibble(network) %>%
            rename("hub_score" = starts_with("hub_score")) %>%
            arrange(desc(hub_score)) %>%
            slice_head(n = 3 + ceiling(nrow(as_tibble(network)) * 0.01)) %>%
            pull(name)

        network_w_labels <- network %>% mutate(
            node_label = case_when(
                degree > label_filter ~ {{label_column}},
                TRUE ~ NA_character_
            ),
            is_hub = case_when(
                name %in% hub_nodes ~ "y",
                TRUE ~ "n"
            )
        )

        ggraph(network_w_labels, layout = layout_object) +
            geom_edge_link(
                show.legend = FALSE,
                edge_alpha = edge_alpha_,
                edge_colour = edge_colour_,
                edge_width = edge_width_
            ) +
            geom_node_point(
                aes(size = degree, fill = new_fill_col),
                pch = 21,
                colour = node_colour
            ) +
            network_fill_geom +
            scale_size_continuous(range = node_size, guide = "none") +
            labs(fill = NULL) +
            theme_tweaks +
            network_fill_guide +
            geom_node_text(
                aes(label = node_label, colour = is_hub),
                size          = label_size,
                repel         = TRUE,
                family        = fontfamily,
                fontface      = label_face,
                check_overlap = TRUE,
                show.legend   = FALSE,
                box.padding   = label_padding,
                min.segment.length = min_seg_length
            ) +
            scale_colour_manual(
                values = c("y" = hub_colour, "n" = label_colour)
            )
    }
}
