#' Visualize enriched Reactome pathways as an interactive network
#'
#' @param network Tidygraph network object as output by `createPathnet`
#' @param networkLayout Desired layout for the network visualization. Defaults
#'   to "layout_nicely", should support most igraph layouts. See
#'   `?visIgraphLayout()` for more details.
#' @param edgeColour Colour of network edges; defaults to "#848484".
#' @param edgeSizeRange Edge width is mapped to the similarity measure (one over
#'   distance). This length-two numeric vector controls the minimum and maximum
#'   width of edges. Defaults to `c(5, 20)`.
#' @param nodeSizeRange Node size is mapped to the negative log of the
#'   Bonferroni-adjusted p value, and this length-two numeric vector controls
#'   the minimum and maximum. Defaults to `c(20, 50)`.
#' @param nodeBorderWidth Size of the node border, defaults to 2.5
#' @param highlighting When clicking on a node, should directly neighbouring
#'   nodes be highlighted (other nodes are dimmed)? Defaults to TRUE.
#' @param labelNodes Boolean determining if nodes should be labeled. Note it
#'   will only ever label enriched nodes/pathways.
#' @param nodeLabelSize Size of the node labels in pixels; defaults to 60.
#' @param nodeLabelColour Colour of the node labels; defaults to "black".
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
#' pathnetVisNetwork(exPathnet)
#'
pathnetVisNetwork <- function(
        network,
        networkLayout = "layout_nicely",
        edgeColour = "#848484",
        edgeSizeRange = c(5, 20),
        nodeSizeRange = c(20, 50),
        nodeBorderWidth = 2.5,
        labelNodes = TRUE,
        nodeLabelSize = 60,
        nodeLabelColour = "black",
        highlighting = TRUE
) {
    visNetNodes <- network %>%
        as_tibble() %>%
        mutate(
            id = row_number(),
            value = if_else(
                is.na(pValueAdjusted),
                1,
                -log10(pValueAdjusted)
            ),
            background = map_chr(
                groupedPathway, ~groupedPathwayColours[[.x]]
            ),
            background = if_else(
                !is.na(pValueAdjusted),
                background,
                "#ffffff"
            ),
            border = map_chr(groupedPathway, ~groupedPathwayColours[[.x]]),
            color = map2(
                background,
                border,
                ~list("background" = .x, "border" = .y)
            )
        ) %>%
        select(
            id,
            "title" = pathwayName1,
            everything(),
            -any_of(c(
                "background", "border", "direction", "pathwayName", "pValue"
            ))
        )

    if (labelNodes) {
        visNetNodes <- mutate(
            visNetNodes,
            label = map_chr(
                if_else(!is.na(pValueAdjusted), title, ""),
                .truncNeatly,
                30
            )
        )
    }

    visnetEdges <- network %>%
        activate("edges") %>%
        as_tibble() %>%
        rename("value" = similarity) %>%
        distinct()

    legendDf <- groupedPathwayColours %>%
        enframe("label", "icon.color") %>%
        mutate(shape = "dot", size = 15)

    out1 <- visNetwork(nodes = visNetNodes, edges = visnetEdges) %>%
        visIgraphLayout(layout = networkLayout) %>%
        visEdges(
            color = edgeColour,
            scaling = list(
                "min" = edgeSizeRange[1],
                "max" = edgeSizeRange[2]
            )
        ) %>%
        visNodes(
            borderWidth = nodeBorderWidth,
            scaling = list(
                "min" = nodeSizeRange[1],
                "max" = nodeSizeRange[2]
            )
        ) %>%
        visOptions(
            highlightNearest = highlighting,
            selectedBy = list(
                "variable" = "groupedPathway",
                "style" = "width: 175px; height: 26px",
                "main" = "Select pathway group"
            )
        ) %>%
        visLegend(
            stepY = 75,
            useGroups = FALSE,
            position = "right",
            main = "Grouped pathway",
            addNodes = legendDf,
        ) %>%
        visExport(float = "right")

    if (labelNodes) {
        out1 %>%
            visNodes(font = paste0(
                nodeLabelSize,
                "px arial ",
                nodeLabelColour
            ))
    } else {
        out1
    }
}
