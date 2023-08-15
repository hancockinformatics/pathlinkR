#' Plot Reactome pathway enrichment results
#'
#' @param enrichResults Data frame of results from the function
#'   `enrichPathway`
#' @param columns Number of columns to split the pathways across, particularly
#'   relevant if there are many significant pathways. Can specify up to 3
#'   columns, with a default of 1.
#' @param specificTopPathways Only plot pathways from a specific vector of
#'   "topPathways". Defaults to "any" which includes all pathway results, or see
#'   `unique(enrichResults$topPathways)` (i.e. the input)
#'   for possible values.
#' @param specificPathways Only plot specific pathways. Defaults to "any".
#' @param colourValues Length-two character vector of colours to use for the
#'   scale. Defaults to `c("white", "steelblue3")`.
#' @param nameWidth How many characters to show for pathway name before
#'   truncating? Defaults to 35.
#' @param nameRows How much to rows to wrap across for the pathway name?
#'   Defaults to 1.
#' @param xAngle Angle of x axis labels, set to "angled" (45 degrees),
#'   "horizontal" (0 degrees), or "vertical" (90 degrees).
#' @param maxPVal P values below `10 ^ -maxPVal` will be set to that value.
#' @param intercepts Add vertical lines to separate different groupings, by
#'   providing a vector of intercepts (e.g. `c(1.5, 2.5)`).
#'   Defaults to `NA`.
#' @param includeGeneRatio Boolean (FALSE). Should the gene ratio be included as
#'   an aesthetic mapping?If so, then it is attributed to the size of the
#'   triangles.
#' @param size Size of points if not scaling to gene ratio. Defaults to 5.
#' @param legendMultiply Size of the legend, e.g. increase if there are a lot of
#'   pathways which makes the legend small and unreadable by comparison.
#'   Defaults to 1, i.e. no increase in legend size.
#' @param showNumGenes Boolean, defaults to FALSE. Show the number of genes for
#'   each comparison as brackets under the comparison's name.
#' @param pathwayPosition Whether to have the y-axis labels (pathway names) on
#'   the left or right side. Default is "right".
#' @param newGroupNames If you want to change the names of the comparisons to
#'   different names. Input a vector in the order as they appear.
#'
#' @return A ggplot object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import stringr
#' @importFrom ggforce facet_col
#' @importFrom ggpubr ggarrange
#'
#' @description Creates a plot to visualize and compare pathway enrichment
#'   results from multiple DE comparisons. Can automatically assign each
#'   pathway into an informative top-level category.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' plotPathways(
#'     sigoraExamples,
#'     columns=2
#' )
#'
plotPathways <- function(
        enrichResults,
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
        size=5,
        legendMultiply=1,
        showNumGenes=FALSE,
        pathwayPosition="right",
        newGroupNames=NA
) {

    ## If new group names are to be used, add them in
    if (!is.na(newGroupNames[1])) {
        mapNames <- tibble(
            comparison=unique(enrichResults$comparison),
            newNames=newGroupNames
        )
        enrichResults <- left_join(
            enrichResults,
            mapNames,
            multiple="all"
        ) %>%
            select(!comparison) %>%
            mutate(comparison=newNames) %>%
            select(!newNames)
    }

    ## Convert to -log10 p value and arbitrarily set to a max -log10 p value of
    ## 50 (i.e. adj pval=10^-50), which some enrichment results surpass,
    ## especially for Sigora.
    enrichResults <- enrichResults %>% mutate(
        logMax=case_when(
            -log10(pValueAdjusted) > maxPVal ~ maxPVal,
            -log10(pValueAdjusted) <= maxPVal ~ -log10(pValueAdjusted)
        )
    )

    ## Order the directionality of results, if up and down are used
    if (!"All" %in% enrichResults$direction) {
        enrichResults$direction <- factor(
            enrichResults$direction,
            levels=c("Up", "Down"))
    } else if (
        "All" %in% enrichResults$direction &
        any(c("Up", "Down") %in% enrichResults$direction)
    ) {
        enrichResults$direction <- factor(
            enrichResults$direction,
            levels=c("Up", "Down", "All"))
    }

    ## Order the comparisons by the order they were inputted (not alphabetical)
    enrichResults$comparison <- factor(
        enrichResults$comparison,
        levels=unique(enrichResults$comparison))

    ## Add in the number of genes for each comparison if indicated
    if (showNumGenes) {
        enrichResults <- enrichResults %>% mutate(
            comparison=paste0(comparison, "\n(", totalGenes, ")")
        )
        enrichResults$comparison <- factor(
            enrichResults$comparison,
            levels=unique(enrichResults$comparison)
        )
    }

    ## If did not specify to only plot specific pathways, otherwise filter them
    if (specificTopPathways[1] == "any") {
        specificTopPathways <- unique(enrichResults$topPathways)
    }
    if (specificPathways[1] == "any") {
        specificPathways <- unique(enrichResults$pathwayName)
    }

    enrichedResultsGraph <- enrichResults %>% filter(
        topPathways %in% specificTopPathways,
        pathwayName %in% specificPathways
    )

    ## In certain cases, a pathway may be enriched by both up- and
    ## down-regulated genes. Find duplicated pathways, and only show the one
    ## that is more enriched (lower p-value)
    duplicates <- enrichedResultsGraph %>%
        group_by(pathwayName, comparison) %>%
        filter(n() > 1) %>%
        mutate(uniqueId=paste0(pathwayName, comparison))

    enrichedResultsGraph <- enrichedResultsGraph %>%
        mutate(uniqueId=paste0(pathwayName, comparison))

    enrichedResultsClean <- enrichedResultsGraph %>%
        filter(!uniqueId %in% duplicates$uniqueId)

    enrichedResultsDupes <- enrichedResultsClean[0, ]

    ## This chooses the top enriched pathway of duplicates and also adds the
    ## other pathway into enrichedResultsDupes.
    if (nrow(duplicates) > 0) {
        message(
            "\nNote: The following pathways were enriched in both directions ",
            "for the given comparisons. These are indicated with an asterisk ",
            "over the triangle, which is only shown for the lower p value ",
            "result:"
        )

        duplicates_message <- duplicates %>%
            mutate(new=paste0("\t", comparison, ": ", pathwayName)) %>%
            pull(new) %>%
            as.character() %>%
            str_remove_all("\n") %>%
            unique()

        message(paste0(
            duplicates_message,
            collapse="\n"
        ))

        for (i in seq_len(nrow(duplicates))) {
            row <- duplicates[i, ]

            choices <- enrichedResultsGraph %>%
                filter(
                    pathwayName == row$pathwayName,
                    comparison == row$comparison
                ) %>%
                arrange(pValueAdjusted) ## Choose the lowest p value

            ## Add the lower p value to the clean dataframe
            enrichedResultsClean <-
                rbind(enrichedResultsClean, choices[1, ])

            ## keep the other enrichment to the dupes dataframe
            enrichedResultsDupes <-
                rbind(enrichedResultsDupes, choices[2, ])
        }
    }

    ## Now organize pathways into multiple columns. Maximum is 3 columns to
    ## graph, if inputted larger, will be set to 3.
    if (columns > 3) {
        message("Maximum is three columns to graph. Plotting three columns.")
        columns <- 3
    }

    ## How many pathways per top pathway? Add 1 to account for the extra space
    ## the facet takes up.
    numPathways <- enrichedResultsClean %>%
        select(topPathways, pathwayName) %>%
        unique() %>%
        group_by(topPathways) %>%
        summarise(pathways=n() + 1)

    ## Arrange ascending
    numPathways <- numPathways %>% arrange(pathways)

    ## The n largest top pathways are chosen to start off the columns
    columnSplitting <- numPathways %>% tail(columns)

    ## For cases when you specify specific pathways, make sure that there are
    ## more pathways than columns...
    if (nrow(columnSplitting) < columns) {
        columns <- nrow(columnSplitting)

    } else {
        for (i in (columns + 1):nrow(numPathways)) {
            ## The "n + 1" largest top pathway to add
            add <- numPathways %>% tail(i) %>% .[1,]

            columnSmallest <- columnSplitting %>%
                arrange(pathways) %>%
                head(1) %>%
                .$topPathways

            ## Now, add the name of the top pathway to one of the columns and
            ## increase the number of pathways
            columnSplitting <- columnSplitting %>% mutate(
                pathways=case_when(
                    topPathways == columnSmallest ~ pathways + add$pathways,
                    TRUE ~ pathways
                ),
                topPathways=case_when(
                    topPathways == columnSmallest ~ paste0(
                        topPathways, "," ,add$topPathways
                    ),
                    TRUE ~ topPathways
                )
            )
        }
    }

    ## Now create a list of top pathways for each column
    columnList <- columnSplitting$topPathways %>%
        as.list() %>%
        str_split(",")

    ## Plot pathways
    plotlist <- list()
    name_trunc <- nameWidth * nameRows - 5

    ## Can be set to angled (45 degrees), "horizontal" (0 degrees), or
    ## "vertical" (90 degrees)
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


    for (n in seq_len(length(columnList))) {
        plot <-
            ggplot(
                data=filter(
                    enrichedResultsClean,
                    topPathways %in% columnList[n][[1]]
                ),
                mapping=aes(
                    x=comparison,
                    y=pathwayName,
                    fill=logMax,
                    shape=direction
                )
            ) +

            facet_col(
                facets=~topPathways,
                scales="free_y",
                space="free"
            ) +

            {if (includeGeneRatio) geom_point(aes(size=gene_ratio))} +
            {if (!includeGeneRatio) geom_point(size=size)} +

            geom_point(
                data=filter(
                    enrichedResultsDupes,
                    topPathways %in% columnList[n][[1]]
                ),
                mapping=aes(x=comparison, y=pathwayName),
                shape=8,
                size=2,
                colour="white",
                show.legend=FALSE
            ) +

            ## Wrap and truncate pathway names if necessary
            scale_y_discrete(
                labels=~str_wrap(
                    .truncNeatly(.x, name_trunc),
                    width=nameWidth
                ),
                position=pathwayPosition
            ) +

            ## Keeps comparisons even if they don"t enrich for any pathways
            scale_x_discrete(drop=FALSE) +

            scale_shape_manual(
                values=c("Down"=25 , "Up"=24, "All"=21),
                name="Regulation",
                na.value=NA,
                drop=FALSE # Keep both up/down if only one direction enriched
            ) +

            scale_fill_continuous(
                name=expression(P[adjusted]),
                limits=c(0, 50),
                breaks=c(10, 20, 30, 40, 50),
                labels=c(
                    expression(10 ^ -10),
                    expression(10 ^ -20),
                    expression(10 ^ -30),
                    expression(10 ^ -40),
                    expression(10 ^ -50)
                ),
                low=colourValues[1],
                high=colourValues[2],
                na.value=NA
            ) +

            ## Can also add lines to separate different groups
            {if (!is.na(intercepts[1])) geom_vline(xintercept=intercepts)} +
            labs(x=NULL, y=NULL) +

            theme_bw() +
            theme(
                strip.text.x=element_text(
                    size=12, face="bold", colour="black"
                ),
                legend.text=element_text(size=12 * legendMultiply),
                legend.title=element_text(size=13 * legendMultiply),
                axis.text.y=element_text(colour="black", size=12),
                axis.text.x=element_text(
                    colour="black",
                    size=12,
                    angle=angle,
                    hjust=hjust,
                    vjust=vjust
                )
            ) +
            guides(
                shape=guide_legend(
                    override.aes=list(size=5 * legendMultiply)
                ),
                size =guide_legend(
                    override.aes=list(shape=24, fill="black")
                )
            )
        plotlist <- append(plotlist, list(plot))
    }

    if (columns > 1) {
        plot <- ggarrange(
            plotlist=plotlist,
            ncol=columns,
            common.legend=TRUE,
            legend="right",
            align="v"
        )
        return(plot)
    } else {
        return(plot)
    }
}
