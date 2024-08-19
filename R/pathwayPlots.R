#' Plot pathway enrichment results
#'
#' @param pathwayEnrichmentResults Data frame of results from the function
#'   `enrichPathway`
#' @param columns Number of columns to split the pathways across, particularly
#'   relevant if there are many significant pathways. Can specify up to 3
#'   columns, with a default of 1.
#' @param specificTopPathways Only plot pathways from a specific vector of
#'   "topLevelPathway". Defaults to "any" which includes all pathway results, or
#'   see `unique(pathwayEnrichmentResults$topLevelPathway)` (i.e. the
#'   input) for possible values.
#' @param specificPathways Only plot specific pathways. Defaults to "any".
#' @param colourValues Length-two character vector of colours to use for the
#'   scale. Defaults to `c("blue", "red")`.
#' @param nameWidth How many characters to show for pathway name before
#'   truncating? Defaults to 35.
#' @param nameRows For pathway names (y axis), how many rows (lines) should
#'   names wrap across when they're too long? Defaults to 1.
#' @param xAngle Angle of x axis labels, set to "angled" (45 degrees),
#'   "horizontal" (0 degrees), or "vertical" (90 degrees).
#' @param maxPVal P values below `10 ^ -maxPVal` will be set to that value.
#' @param intercepts Add vertical lines to separate different groupings, by
#'   providing a vector of intercepts (e.g. `c(1.5, 2.5)`).
#'   Defaults to `NA`.
#' @param includeGeneRatio Boolean (FALSE). Should the gene ratio be included as
#'   an aesthetic mapping?If so, then it is attributed to the size of the
#'   triangles.
#' @param size Size of points if not scaling to gene ratio. Defaults to 4.
#' @param legendMultiply Size of the legend, e.g. increase if there are a lot of
#'   pathways which makes the legend small and unreadable by comparison.
#'   Defaults to 1, i.e. no increase in legend size.
#' @param showNumGenes Boolean, defaults to FALSE. Show the number of genes for
#'   each comparison as brackets under the comparison's name.
#' @param pathwayPosition Whether to have the y-axis labels (pathway names) on
#'   the left or right side. Default is "right".
#' @param newGroupNames If you want to change the names of the comparisons to
#'   different names. Input a vector in the order as they appear.
#' @param fontSize Base font size for all text elements of the plot. Defaults to
#'   12.
#'
#' @return A plot of enriched pathways; a "ggplot" object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#'
#' @importFrom stringr str_remove_all str_split str_wrap
#'
#' @description Creates a plot to visualize and compare pathway enrichment
#'   results from multiple DE comparisons. Can automatically assign each pathway
#'   into an informative top-level category.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'          <https://bioconductor.org/packages/fgsea/>
#'
#' @examples
#' data("sigoraExamples")
#' pathwayPlots(sigoraExamples, columns=2)
#'
pathwayPlots <- function(
        pathwayEnrichmentResults,
        columns=1,
        specificTopPathways="any",
        specificPathways="any",
        colourValues=c("blue", "red"),
        nameWidth=35,
        nameRows=1,
        xAngle="angled",
        maxPVal=50,
        intercepts=NA,
        includeGeneRatio=FALSE,
        size=4,
        legendMultiply=1,
        showNumGenes=FALSE,
        pathwayPosition="right",
        newGroupNames=NA,
        fontSize=12
) {

    stopifnot("A maximum of three columns can be specified."=columns <= 3)
    stopifnot(
        "'xAngle' should be a string; see `?pathwayPlots` for options"=
            is.character(xAngle)
    )

    plotData <- pathwayEnrichmentResults

    fgseaColumns <- c("log2err", "ES", "NES", "size", "leadingEdge")

    if (all(fgseaColumns %in% colnames(plotData))) {
        stopifnot(
            "Option 'includeGeneRatio' is not supported for fgsea results."=
                !includeGeneRatio
        )
    }

    ## If new group names are to be used, add them in
    if (!is.na(newGroupNames[1])) {
        plotData <- left_join(
            plotData,
            tibble::tibble(
                comparison=unique(plotData$comparison),
                newNames=newGroupNames
            ),
            multiple="all"
        ) %>%
            select(-comparison) %>%
            rename("comparison"=newNames)
    }

    ## If indicated add the number of genes to each comparison
    if (showNumGenes) {
        plotData <- plotData %>% mutate(
            comparison=paste0(comparison, "\n(", totalGenes, ")")
        )
    }

    plotDataGraph <- plotData %>%
        mutate(
            logMax=case_when(
                -log10(pValueAdjusted) > maxPVal ~ maxPVal,
                -log10(pValueAdjusted) <= maxPVal ~ -log10(pValueAdjusted)
            ),
            direction=factor(
                direction,
                rev(sort(unique(plotData$direction)))
            ),
            comparison=factor(
                comparison,
                unique(plotData$comparison)
            )
        )

    if (specificTopPathways[1] != "any") {
        plotDataGraph <- plotDataGraph %>%
            filter(topLevelPathway %in% specificTopPathways)
    }
    if (specificPathways[1] != "any") {
        plotDataGraph <- plotDataGraph %>%
            filter(pathwayName %in% specificPathways)
    }

    ## How many pathways per top pathway? Add 1 to account for the extra space
    ## the facet takes up.
    numPathways <- plotDataGraph %>%
        distinct(topLevelPathway, pathwayName) %>%
        group_by(topLevelPathway) %>%
        summarise(pathways=n() + 1) %>%
        arrange(pathways)

    ## In the case where a pathway is enriched in both up- and down- regulated
    ## genes, only show the one with a lower p-value
    plotDataDups <- plotDataGraph %>%
        group_by(pathwayName, comparison) %>%
        filter(n() > 1)

    plotDataNoDups <- plotDataGraph %>%
        group_by(pathwayName, comparison) %>%
        arrange(pValueAdjusted) %>%
        distinct(pathwayName, comparison, .keep_all=TRUE)

    ## The n largest top pathways are chosen to start off the columns. For cases
    ## when you specify specific pathways, make sure that there are more
    ## pathways than columns...
    columnSplitting <- numPathways %>% tail(columns)

    if (nrow(columnSplitting) < columns) {
        columns <- nrow(columnSplitting)

    } else {
        for (i in seq(columns + 1, nrow(numPathways))) {
            ## The "n + 1" largest top pathway to add
            add <- numPathways %>% tail(i) %>% .[1,]

            columnSmallest <- columnSplitting %>%
                arrange(pathways) %>%
                head(1) %>%
                .$topLevelPathway

            ## Now, add the name of the top pathway to one of the columns and
            ## increase the number of pathways
            columnSplitting <- columnSplitting %>% mutate(
                pathways=case_when(
                    topLevelPathway == columnSmallest ~ pathways + add$pathways,
                    TRUE ~ pathways
                ),
                topLevelPathway=case_when(
                    topLevelPathway == columnSmallest ~ paste0(
                        topLevelPathway, "," ,add$topLevelPathway
                    ),
                    TRUE ~ topLevelPathway
                )
            )
        }
    }

    columnList <- columnSplitting$topLevelPathway %>%
        as.list() %>%
        str_split(",")

    ## Can be "horizontal" (0), "angled" (45 degrees), or "vertical" (90)
    if (xAngle == "angled") {
        vjust <- 1
        if (pathwayPosition == "right") {
            angle <- -45
            hjust <- 0
        } else if (pathwayPosition == "left") {
            angle <- 45
            hjust <- 1
        }
    } else if (xAngle == "horizontal") {
        angle <- 0
        hjust <- 0.5
        vjust <- 1
    } else if (xAngle == "vertical") {
        angle <- 90
        hjust <- 0.5
        vjust <- 0.5
    }

    themePathway <- theme_bw(base_size = fontSize) +
        theme(
            strip.text.x=element_text(face="bold", colour="black"),
            legend.text=element_text(size=(fontSize - 2) * legendMultiply),
            legend.title=element_text(
                size=fontSize * legendMultiply,
                margin=margin(r=10, b=7)
            ),
            axis.text.y=element_text(colour="black"),
            axis.text.x=element_text(
                colour="black",
                angle=angle,
                hjust=hjust,
                vjust=vjust
            )
        )

    ## Plot pathways
    plotList <- lapply(columnList, function(x) {
        ggplot(
            data=filter(
                plotDataNoDups,
                topLevelPathway %in% x
            ),
            mapping=aes(
                x=comparison,
                y=pathwayName,
                fill=logMax,
                shape=direction
            )
        ) +

            ggforce::facet_col(
                facets=~topLevelPathway,
                scales="free_y",
                space="free"
            ) +

            {
                if (includeGeneRatio) {
                    geom_point(aes(size=geneRatio))
                } else {
                    geom_point(size=size)
                }
            } +

            geom_point(
                data=filter(
                    plotDataDups,
                    topLevelPathway %in% x
                ),
                mapping=aes(x=comparison, y=pathwayName),
                shape=8,
                size=ifelse(size - 3 > 0, size - 3, 1),
                colour="white",
                show.legend=FALSE
            ) +

            ## Keep comparisons even if they don't enrich for any pathways
            scale_x_discrete(drop=FALSE) +

            ## Wrap and truncate pathway names if necessary
            scale_y_discrete(
                labels=~str_wrap(
                    .truncNeatly(.x, (nameWidth * nameRows) - 5),
                    width=nameWidth
                ),
                position=pathwayPosition
            ) +

            scale_shape_manual(
                values=c("Down"=25, "Up"=24, "All"=21),
                name="Regulation",
                na.value=NA,
                drop=FALSE # Keep both up/down if only one direction enriched
            ) +

            { if (requireNamespace("scales", quietly=TRUE)) {
                scale_fill_continuous(
                    name=expression(P[adjusted]),
                    labels=scales::label_math(10^-.x),
                    low=colourValues[1],
                    high=colourValues[2],
                    na.value=NA
                )
            } else {
                scale_fill_continuous(
                    name=expression(-log10(P[adjusted])),
                    low=colourValues[1],
                    high=colourValues[2],
                    na.value=NA
                )
            }} +

            ## Add optional lines to separate different groups
            {if (!is.na(intercepts[1])) geom_vline(xintercept=intercepts)} +

            labs(x=NULL, y=NULL) +
            themePathway +
            guides(
                shape=guide_legend(override.aes=list(size=4 * legendMultiply)),
                size =guide_legend(override.aes=list(shape=24, fill="black"))
            )
    })

    if (columns > 1) {
        ggpubr::ggarrange(
            plotlist=plotList,
            ncol=columns,
            common.legend=TRUE,
            legend="right",
            align="v"
        )
    } else {
        plotList[[1]]
    }
}
