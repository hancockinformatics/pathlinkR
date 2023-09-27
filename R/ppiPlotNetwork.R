#' Plot an undirected PPI network using ggraph
#'
#' @param network A `tidygraph` object, output from `ppiBuildNetwork`
#' @param networkLayout Layout of nodes in the network. Supports all layouts
#'   from `ggraph`/`igraph`, or a data frame of x and y coordinates for
#'   each node (order matters!).
#' @param title Optional title for the plot (NA)
#' @param nodeSize Length-two numeric vector, specifying size range of node
#'   sizes (maps to node degree). Default is `c(2, 6)`.
#' @param fillColumn Tidy-select column for mapping node colour. Designed to
#'   handle continuous numeric mappings (either positive/negative only, or
#'   both), and categorical mappings, plus a special case for displaying fold
#'   changes from, for example, RNA-Seq data. See `fillType` for more details
#'   on how to set this up.
#' @param fillType String denoting type of fill mapping to perform for nodes.
#'   Options are: "foldChange", "twoSided", "oneSided", or "categorical".
#' @param catFillColours Colour palette to be used when `fillType` is set
#'   to "categorical." Defaults to "Set1" from RColorBrewer. Will otherwise be
#'   passed as the "values" argument in `scale_fill_manual()`.
#' @param foldChangeColours A two-length character vector containing colours
#'   for up and down regulated genes. Defaults to `c("firebrick3", "#188119")`.
#' @param intColour Fill colour for non-seed nodes, i.e. interactors. Defaults
#'   to "grey70".
#' @param nodeBorder Colour (stroke or outline) of all nodes in the network.
#'   Defaults to "grey30".
#' @param hubColour Colour of node labels for hubs. The top 2% of nodes
#'   (based on calculated hub score) are highlighted with this colour, if
#'   `label=TRUE`.
#' @param subnetwork Logical determining if networks from
#'   `ppiExtractSubnetwork()` should be treated as such. Defaults to TRUE.
#' @param legend Should a legend be included? Defaults to FALSE.
#' @param legendTitle Optional title for the legend, defaults to `NULL`.
#' @param edgeColour Edge colour, defaults to "grey40"
#' @param edgeAlpha Transparency of edges, defaults to 0.5
#' @param edgeWidth Thickness of edges connecting nodes. Defaults to 0.5
#' @param label Boolean, whether labels should be added to nodes. Defaults to
#'   FALSE.
#' @param labelColumn Tidy-select column of the network/data to be used in
#'   labeling nodes. Recommend setting to `hgncSymbol`, which contains HGNC
#'   symbols mapped from the input Ensembl IDs via biomaRt.
#' @param labelFilter Degree filter used to determine which nodes should be
#'   labeled. Defaults to 0. This value can be increased to reduce the number of
#'   node labels, to prevent the network from being too crowded.
#' @param labelSize Size of node labels, defaults to 5.
#' @param labelColour Colour of node labels, defaults to "black"
#' @param labelFace Font face for node labels, defaults to "bold"
#' @param labelPadding Padding around the label, defaults to 0.25 lines.
#' @param minSegLength Minimum length of lines to be drawn from labels to
#'   points. The default specified here is 0.25, half of the normal default
#'   value.
#'
#' @return A Protein-Protein Interaction (PPI) network plot; an object of class
#'   "ggplot"
#' @export
#'
#' @import dplyr
#' @import ggplot2
#'
#' @importFrom ggraph geom_edge_link geom_node_point geom_node_text ggraph
#'   set_graph_style
#'
#' @description Visualize a protein-protein interaction (PPI) network using
#'   `ggraph` functions, output from `ppiBuildNetwork`.
#'
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
#'   The `fillType` argument will determine how the node colour is mapped to
#'   the desired column. "foldChange" represents a special case, where the fill
#'   column is numeric and whose values should be mapped to up (> 0) or down (<
#'   0). "twoSided" and "oneSided" are designed for numeric data that contains
#'   either positive and negative values, or only positive/negative values,
#'   respectively. "categorical" handles any other non-numeric colour mapping,
#'   and uses "Set1" from RColorBrewer.
#'
#'   Node statistics (degree, betweenness, and hub score) are calculated using
#'   the respective functions from the `tidygraph` package.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR/>
#'
#' @examples
#' exNetwork <- ppiBuildNetwork(
#'     deseqResults=exampleDESeqResults[[1]],
#'     filterInput=TRUE,
#'     order="zero"
#' )
#'
#' ppiPlotNetwork(
#'     network=exNetwork,
#'     title="COVID positive over time",
#'     fillColumn=log2FoldChange,
#'     fillType="foldChange",
#'     legend=TRUE,
#'     label=TRUE,
#'     labelColumn=hgncSymbol,
#'     labelFilter=5
#' )
#'
ppiPlotNetwork <- function(
        network,
        networkLayout="nicely",
        title=NA,
        nodeSize=c(2, 6),
        fillColumn,
        fillType,
        catFillColours="Set1",
        foldChangeColours=c("firebrick3", "#188119"),
        intColour="grey70",
        nodeBorder="grey30",
        hubColour="blue2",
        subnetwork=TRUE,
        legend=FALSE,
        legendTitle=NULL,
        edgeColour="grey40",
        edgeAlpha=0.5,
        edgeWidth=0.5,
        label=FALSE,
        labelColumn,
        labelFilter=8,
        labelSize=4,
        labelColour="black",
        labelFace="bold",
        labelPadding=0.25,
        minSegLength=0.25
) {

    stopifnot(is(network, "tbl_graph"))

    ## Set up fill scaling based on argument `fillType`
    if (fillType == "foldChange") {
        stopifnot(is.numeric(pull(network, {{fillColumn}})))

        network <- network %>% mutate(
            newFillCol=case_when(
                {{fillColumn}} < 0 ~ "Down",
                {{fillColumn}} > 0 ~ "Up",
                TRUE ~ NA_character_
            )
        )

        networkFillGeom <- scale_fill_manual(
            values=c(
                "Up"=foldChangeColours[1],
                "Down"=foldChangeColours[2]
            ),
            na.value=intColour
        )

        networkFillGuide <- guides(
            fill=guide_legend(
                title="Direction",
                override.aes=list(size=5)
            )
        )

    } else if (fillType == "twoSided") {

        network <- mutate(network, newFillCol={{fillColumn}})
        networkFillGeom <- scale_fill_gradient2(
            low="#313695",
            mid="white",
            high="#a50026",
            midpoint=0,
            na.value=intColour,
            guide=ifelse(legend, "colourbar", "none")
        )
        networkFillGuide <- NULL

    } else if (fillType == "oneSided") {

        network <- mutate(network, newFillCol={{fillColumn}})
        networkFillGeom <- scale_fill_viridis_c(
            option="plasma", begin=0.2
        )
        networkFillGuide <- NULL

    } else if (fillType == "categorical") {

        network <- mutate(network, newFillCol={{fillColumn}})

        if (all(catFillColours == "Set1")) {
            networkFillGeom <- scale_fill_brewer(
                palette ="Set1",
                na.value=intColour,
                guide   =ifelse(legend, "legend", "none")
            )
        } else {
            networkFillGeom <- scale_fill_manual(values=catFillColours)
        }
        networkFillGuide <-
            guides(fill=guide_legend(override.aes=list(size=5)))

    } else {
        stop(
            "Argument 'fillType' must be one of 'foldChange', 'twoSided', ",
            "'oneSided', or 'categorical'"
        )
    }

    if (is.data.frame(networkLayout)) {
        message("Using user-supplied node coordinates...")
        ## By converting the layout object to a matrix, we no longer need to
        ## worry about column names. The first and second column will be "x" and
        ## "y", respectively.
        layoutObject <- as.matrix(networkLayout)
    } else {
        layoutObject <- networkLayout
    }

    ## Set a plain white background
    set_graph_style(foreground="white", family="sans")

    ## Theme tweaks for all plot types
    themeTweaks <- theme(
        plot.margin=unit(c(0.5, 0, 0, 0.5), "cm"),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14)
    )


    ## Subnetwork stuff
    if (subnetwork & "starters" %in% names(attributes(network))) {
        message(
            "Detected this is a subnetwork generated by ",
            "'ppiExtractSubnetwork()'. Highlighted node labels indicate genes ",
            "from the extracted pathway or gene list."
        )

        starterNodes <- attr(network, "starters")

        network <- network %>%
            mutate(
                nodeLabel = case_when(
                    degree > labelFilter ~ {{labelColumn}},
                    TRUE ~ NA_character_
                ),
                isStarter = if_else(name %in% starterNodes, "y", "n")
            )

        if (hubColour == "blue2") {
            hubColour <- "red"
        }

        ggraph(network, layout = layoutObject) +
            geom_edge_link(
                show.legend=FALSE,
                edge_alpha=edgeAlpha,
                edge_colour=edgeColour,
                edge_width=edgeWidth
            ) +
            geom_node_point(
                aes(size=degree, fill=newFillCol),
                pch=21,
                colour=nodeBorder
            ) +
            networkFillGeom +
            scale_size_continuous(range=nodeSize, guide="none") +
            scale_colour_manual(values=c("y"=hubColour, "n"=labelColour)) +
            labs(fill=legendTitle) +
            { if (!is.na(title)) labs(title=title) } +
            geom_node_text(
                aes(label=nodeLabel, colour=isStarter),
                size=labelSize,
                bg.color="white",
                repel=TRUE,
                fontface=labelFace,
                check_overlap=TRUE,
                show.legend=FALSE,
                box.padding=labelPadding,
                min.segment.length=minSegLength
            ) +
            themeTweaks +
            networkFillGuide

    } else if (label) {
        hubNodes <- tibble::as_tibble(network) %>%
            rename("hubScore"=starts_with("hubScore")) %>%
            arrange(desc(hubScore)) %>%
            slice_head(
                n=3 + ceiling(nrow(tibble::as_tibble(network)) * 0.01)
            ) %>%
            pull(name)

        networkLabeled <- network %>% mutate(
            nodeLabel=case_when(
                degree > labelFilter ~ {{labelColumn}},
                TRUE ~ NA_character_
            ),
            isHub=case_when(
                name %in% hubNodes ~ "y",
                TRUE ~ "n"
            )
        )

        ggraph(networkLabeled, layout=layoutObject) +
            geom_edge_link(
                show.legend=FALSE,
                edge_alpha=edgeAlpha,
                edge_colour=edgeColour,
                edge_width=edgeWidth
            ) +
            geom_node_point(
                aes(size=degree, fill=newFillCol),
                pch=21,
                colour=nodeBorder
            ) +
            networkFillGeom +
            scale_size_continuous(range=nodeSize, guide="none") +
            labs(fill=legendTitle) +
            { if (!is.na(title)) labs(title=title) } +
            themeTweaks +
            networkFillGuide +
            geom_node_text(
                aes(label=nodeLabel, colour=isHub),
                size=labelSize,
                bg.color="white",
                repel=TRUE,
                fontface=labelFace,
                check_overlap=TRUE,
                show.legend=FALSE,
                box.padding=labelPadding,
                min.segment.length=minSegLength
            ) +
            scale_colour_manual(values=c("y"=hubColour, "n"=labelColour))

    } else {
        ggraph(network, layout=layoutObject) +
            geom_edge_link(
                show.legend=FALSE,
                edge_alpha=edgeAlpha,
                edge_colour=edgeColour,
                edge_width=edgeWidth
            ) +
            geom_node_point(
                aes(size=degree, fill=newFillCol),
                pch=21,
                colour=nodeBorder
            ) +
            networkFillGeom +
            scale_size_continuous(range=nodeSize, guide="none") +
            labs(fill=legendTitle) +
            { if (!is.na(title)) labs(title=title) } +
            themeTweaks +
            networkFillGuide
    }
}
