#' Visualize enriched Reactome pathways as an interactive network
#'
#' @param network Tidygraph network object as output by `createPathnet`
#' @param networkLayout Desired layout for the network visualization. Defaults
#'   to "layout_nicely", and should support most igraph layouts. See
#'   `?visIgraphLayout` for more details.
#' @param nodeSizeRange Node size is mapped to the negative log of the
#'   Bonferroni-adjusted p value, and this length-two numeric vector controls
#'   the minimum and maximum. Defaults to `c(20, 50)`.
#' @param nodeBorderWidth Size of the node border, defaults to 2.5
#' @param labelNodes Boolean determining if nodes should be labeled. Note it
#'   will only ever label enriched nodes/pathways.
#' @param nodeLabelSize Size of the node labels in pixels; defaults to 60.
#' @param nodeLabelColour Colour of the node labels; defaults to "black".
#' @param edgeColour Colour of network edges; defaults to "#848484".
#' @param edgeWidthRange Edge width is mapped to the similarity measure (one
#'   over distance). This length-two numeric vector controls the minimum and
#'   maximum width of edges. Defaults to `c(5, 20)`.
#' @param highlighting When clicking on a node, should directly neighbouring
#'   nodes be highlighted (other nodes are dimmed)? Defaults to TRUE.
#'
#' @return Interactive visNetwork plot
#' @export
#'
#' @import dplyr
#'
#' @importFrom purrr map2 map_chr
#' @importFrom tibble enframe
#' @importFrom tidygraph activate
#' @importFrom visNetwork visEdges visExport visIgraphLayout visLegend
#'   visNetwork visNodes visOptions
#'
#' @description Plots the network object generated from `createPathnet`,
#'   creating a visual and interactive representation of similarities/
#'   interactions between pathways using their overlapping genes.
#'
#' @details  This function makes use of the visNetwork library, which allows for
#'   various forms of interactivity, such as including text when hovering over
#'   nodes, node selection and dragging (including multiple selections), and
#'   highlighting nodes belonging to a larger group (e.g. top-level Reactome
#'   category).
#'
#' @references <https://datastorm-open.github.io/visNetwork/>
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' startingPathways <- pathnetFoundation(
#'     mat=pathwayDistancesJaccard,
#'     maxDistance=0.8
#' )
#'
#' exPathnet <- pathnetCreate(
#'     pathwayEnrichmentResult=sigoraExamples[
#'         grepl("Pos", sigoraExamples$comparison),
#'     ],
#'     foundation=startingPathways,
#'     trim=TRUE,
#'     trimOrder=1
#' )
#'
#' pathnetVisNetwork(exPathnet)
#'
pathnetVisNetwork <- function(
        network,
        networkLayout="layout_nicely",
        nodeSizeRange=c(20, 50),
        nodeBorderWidth=2.5,
        labelNodes=TRUE,
        nodeLabelSize=60,
        nodeLabelColour="black",
        edgeColour="#848484",
        edgeWidthRange=c(5, 20),
        highlighting=TRUE
) {

    stopifnot("'nodeSizeRange' should be a length-two numeric vector" = {
        length(nodeSizeRange) == 2
    })

    stopifnot("'edgeWidthRange' should be a length-two numeric vector" = {
        length(edgeWidthRange) == 2
    })

    visNetNodes <- network %>%
        tibble::as_tibble() %>%
        mutate(
            id=row_number(),
            value=if_else(
                is.na(pValueAdjusted),
                1,
                -log10(pValueAdjusted)
            ),
            background=map_chr(groupedPathway, ~groupedPathwayColours[[.x]]),
            background=if_else(!is.na(pValueAdjusted), background, "#ffffff"),
            border=map_chr(groupedPathway, ~groupedPathwayColours[[.x]]),
            color=map2(background, border, ~list("background"=.x, "border"=.y))
        ) %>%
        select(
            id,
            "title"=pathwayName,
            everything(),
            -any_of(c(
                "background", "border", "direction", "pValue"
            ))
        )

    if (labelNodes) {
        visNetNodes <- mutate(
            visNetNodes,
            label=map_chr(
                if_else(!is.na(pValueAdjusted), title, ""),
                .truncNeatly,
                30
            )
        )
    }

    visnetEdges <- network %>%
        activate("edges") %>%
        tibble::as_tibble() %>%
        rename("value"=similarity) %>%
        distinct()

    legendDf <- groupedPathwayColours %>%
        enframe("label", "icon.color") %>%
        mutate(shape="dot", size=15)

    out1 <- visNetwork(nodes=visNetNodes, edges=visnetEdges) %>%
        visIgraphLayout(layout=networkLayout) %>%
        visEdges(
            color=edgeColour,
            scaling=list("min"=edgeWidthRange[1], "max"=edgeWidthRange[2])
        ) %>%
        visNodes(
            borderWidth=nodeBorderWidth,
            scaling=list("min"=nodeSizeRange[1], "max"=nodeSizeRange[2])
        ) %>%
        visOptions(
            highlightNearest=highlighting,
            selectedBy=list(
                "variable"="groupedPathway",
                "style"="width: 175px; height: 26px",
                "main"="Select pathway group"
            )
        ) %>%
        visLegend(
            stepY=75,
            useGroups=FALSE,
            position="right",
            main="Grouped pathway",
            addNodes=legendDf,
        ) %>%
        visExport(float="right")

    if (labelNodes) {
        out1 %>%
            visNodes(font=paste0(nodeLabelSize, "px arial ", nodeLabelColour))
    } else {
        out1
    }
}
