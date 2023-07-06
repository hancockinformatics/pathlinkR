#' Visualize enriched Reactome pathways as a network
#'
#' @param network Tidygraph network object, output from `createPathnet`. See
#'  Details for specific requirements.
#' @param networkLayout Desired layout for the network visualization. Defaults to
#'   "nicely", but supports any method found in `?layout_tbl_graph_igraph`
#' @param edgeColour Colour of network edges; defaults to "grey30".
#' @param edgeAlpha Alpha value for edges; defaults to `1`.
#' @param nodeSizeRange Size range for nodes, mapped to significance
#'   (Bonferroni p-value). Defaults to `c(4, 8)`.
#' @param edgeWidthRange Range of edge widths, mapped to
#' `log10(similarity)`. Defaults to `c(0.33, 3)`.
#' @param labelProp Proportion of "interactor" (i.e. non-enriched) pathways
#'   that the function will attempt to label. E.g. setting this to 0.5 (the
#'   default) means half of the non-enriched pathways will *potentially* be
#'   labeled - it won't be exact because the node labeling is done with
#'   `ggrepel`.
#' @param nodeLabelSize Size of node labels; defaults to 5.
#' @param nodeLabelAlpha Transparency of node labels. Defaults to `0.67`.
#' @param nodeLabelOverlaps Max overlaps for node labels, from `ggrepel`.
#'   Defaults to `6`.
#' @param segColour Colour of line segments connecting labels to nodes.
#'   Defaults to "black".
#' @param themeBaseSize Base font size for all plot elements. Defaults to
#'   `16`.
#'
#' @return An object of class "gg"
#' @export
#'
#' @import ggplot2
#' @import ggraph
#' @import dplyr
#'
#' @description Plots the network object generated from `createPathnet`
#'
#' @details A note regarding node labels: The function tries to prioritize
#'   labeling enriched pathways (filled nodes), with the `labelProp` argument
#'   determining roughly how many of the remaining interactor pathways might get
#'   labels. You'll likely need to tweak this value, and try different seeds, to
#'   get the desired effect.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' startingPathways <- createFoundation(
#'     mat = pathwayDistancesJaccard,
#'     maxDistance = 0.8
#' )
#'
#' exPathnet <- createPathnet(
#'     sigoraResult = sigoraExamples,
#'     foundation = startingPathways,
#'     trim = TRUE,
#'     trimOrder = 1
#' )
#'
#' pathnetGGraph(
#'     exPathnet,
#'     labelProp = 0.1,
#'     nodeLabelSize = 4,
#'     nodeLabelOverlaps = 8,
#'     segColour = "red"
#' )
#'
pathnetGGraph <- function(
        network,
        networkLayout = "nicely",
        nodeSizeRange = c(4, 8),
        edgeColour = "grey30",
        edgeAlpha = 1,
        edgeWidthRange = c(0.33, 3),
        labelProp = 0.25,
        nodeLabelSize = 5,
        nodeLabelAlpha = 0.67,
        nodeLabelOverlaps = 6,
        segColour = "black",
        themeBaseSize = 16
) {
    # Check column names for both nodes and edges
    stopifnot(all(
        c(
            "pathwayName1",
            "pValueAdjusted",
            "groupedPathway"
        ) %in% colnames(as_tibble(network))
    ))

    stopifnot(all(
        "similarity" %in% colnames(
            as_tibble(tidygraph::activate(network, "edges"))
        )
    ))

    interactorsAll <- network %>%
        filter(is.na(pValueAdjusted)) %>%
        pull(pathwayName1)

    interactorsToLabel <- sample(
        interactorsAll,
        size = (length(interactorsAll) * labelProp),
        replace = FALSE
    )

    networkToPlot <- network %>% mutate(
        nodeFill = if_else(
            !is.na(pValueAdjusted),
            groupedPathway,
            NA_character_
        ),
        nodeLabel = case_when(
            !is.na(pValueAdjusted) ~ pathwayName1,
            pathwayName1 %in% interactorsToLabel ~ pathwayName1,
            TRUE ~ NA_character_
        ),
        nodeLabel = map_chr(
            nodeLabel,
            ~.truncNeatly(.x, l = 40) %>% str_wrap(width = 20)
        ),
        pValueAdjusted = if_else(
            !is.na(pValueAdjusted),
            pValueAdjusted,
            1
        )
    )

    ggraph(networkToPlot, layout = networkLayout) +
        # Edges
        geom_edge_link(
            aes(edge_width = log10(similarity)),
            colour = edgeColour,
            alpha = edgeAlpha
        ) +
        scale_edge_width(range = edgeWidthRange, name = "Similarity") +

        # Nodes
        geom_node_point(
            aes(
                size = -log10(pValueAdjusted),
                fill = nodeFill,
                colour = groupedPathway
            ),
            pch = 21,
            stroke = 1.5
        ) +
        scale_size_continuous(
            labels = scales::label_math(expr = 10^-~.x),
            range = nodeSizeRange
        ) +
        scale_fill_manual(
            values = groupedPathwayColours,
            na.value = "white",
            guide = NULL
        ) +
        scale_colour_manual(values = groupedPathwayColours) +

        # Node labels
        geom_node_label(
            aes(label = nodeLabel),
            repel = TRUE,
            size = nodeLabelSize,
            alpha = nodeLabelAlpha,
            min.segment.length = 0,
            segment.colour = segColour,
            max.overlaps = nodeLabelOverlaps
        ) +

        # Misc
        labs(
            size = "Bonferroni\np-value",
            colour = "Pathway type"
        ) +
        theme_void(base_size = themeBaseSize) +
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
