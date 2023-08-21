#' Create a heatmap of fold changes to visualize DESeq2 results
#'
#' @param inputList List of data frames, each the output from
#'   `DESeq2::results()`. The list names are used as the comparison name for
#'   each dataframe (e.g. "COVID vs Healthy"). Data frames should have Ensembl
#'   gene IDs as rownames.
#' @param pathName The name of a Reactome pathway to pull genes from, also used
#'   for the plot title. Alternative to `pathID`.
#' @param pathId ID of a Reactome pathway to pull genes from. Alternative to
#'   `pathName`.
#' @param genesToPlot Vector of Ensembl gene IDs you want to plot, instead of
#'   pulling the genes from a pathway, i.e. this option and
#'   `pathName`/`pathID` are mutually exclusive.
#' @param manualTitle Provide your own title, and override the use of a pathway
#'   name the title.
#' @param titleSize Font size for the title.
#' @param geneFormat Type of genes given in `genesToPlot`. Default is Ensembl
#'   gene IDs ("ensembl"), but can also input a vector of HGNC symbols ("hgnc").
#' @param pCutoff P value cutoff, default is <0.05
#' @param fcCutoff Absolute fold change cutoff, default is >1.5
#' @param cellColours Vector specifying desired colours to use for the cells in
#'   the heatmap. Defaults to `c("blue", "white", "red")`.
#' @param cellBorder A call to `grid::gpar()` to specify borders between
#'   cells in the heatmap. The default is `gpar(col="grey")`. To
#'   remove borders set to `gpar(col=NA)`
#' @param plotSignificantOnly Boolean (TRUE). Only plot genes that are
#'   differentially expressed (i.e. they pass `pCutoff` and `fcCutoff`) in
#'   any comparison from the provided list of data frames.
#' @param showStars Boolean (TRUE) show significance stars on the heatmap
#' @param hideNonsigFC Boolean (TRUE). If a gene is significant in one
#'   comparison but not in another, this will set the colour of the non-
#'   significant gene as grey to visually emphasize the significant genes. If
#'   set to FALSE, it will be set the colour to the fold change, and if the p
#'   value passes `pCutoff`, it will also display the p value (the asterisks
#'   will be grey instead of black).
#' @param vjust Adjustment of the position of the significance stars. Default
#'   is 0.75. May need to adjust if there are many genes.
#' @param rot Rotation of the position of the significance stars. Default is 0.
#' @param invert Boolean (FALSE). The default setting plots genes as rows and
#'   comparisons as columns. Setting this to `TRUE` will place genes as columns
#'   and comparisons as rows.
#' @param log2FoldChange Boolean (FALSE). Default plots the fold changes in the
#'   legend as the true fold change. Set to TRUE if you want log2 fold change.
#' @param colSplit A vector, with the same length as `inputList`, which
#'   assigns each data frame in `inputList` to a group, and splits the heatmap
#'   on these larger groupings. The order of groups in the heatmap will be
#'   carried over, so one can alter the order of `inputList` and
#'   `colSplit` to affect the heatmap. This argument will be ignored if
#'   `clusterColumns` is set to TRUE. See Details for more information.
#' @param clusterRows Boolean (TRUE). Whether to cluster the rows (genes). May
#'   need to change if `invert=TRUE`.
#' @param clusterColumns Boolean (FALSE). Whether to cluster the columns
#'   (comparisons). Will override order of `colSplit` if set to TRUE. May need
#'   to change if `invert=TRUE`.
#' @param colAngle Angle of column text. Defaults to 90.
#' @param colCenter Whether to center column text. Default is TRUE, but it
#'   should be set to FALSE if the column name is angled (e.g. `colAngle=45`).
#' @param rowAngle Angle of row text, defaults to 0.
#' @param rowCenter Whether to center column text. The default is FALSE, but it
#'   should be set to TRUE if vertical column name (e.g. `rowAngle=90`).
#'
#' @return A heatmap of fold changes for genes of interest
#' @export
#'
#' @import dplyr
#' @import purrr
#' @importFrom ComplexHeatmap draw Heatmap
#' @importFrom grid grid.text gpar
#'
#' @description Creates a heatmap of fold changes values for results from the
#'   `DESeq2::results()` function, with various parameters to tweak the
#'   appearance.
#'
#' @details The `cellColours` argument is designed to map a range of negative
#'   and positive values to the three provided colours, with zero as the middle
#'   colour. If the plotted matrix contains only positive (or negative) values,
#'   then it will become a two-colour scale, white-to-red (or blue-to-white).
#'
#' The `colSplit` argument can be used to define larger groups represented
#' in `inputList`. For example, consider an experiment comparing two
#' different treatments to an untreated control, in both wild type and mutant
#' cells. This would give the following comparisons:
#' "wildtype_treatment1_vs_untreated", "wildtype_treatment2_vs_untreated",
#' "mutant_treatment1_vs_untreated", and "mutant_treatment2_vs_untreated". One
#' could then specify `colSplit` as
#' `c("Wild type", "Wild type", "Mutant", "Mutant")` to make the wild
#' type and mutant results more visually distinct.
#'
#' @references <https://bioconductor.org/packages/ComplexHeatmap/>
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' plotFoldChange(
#'     deseqExampleList,
#'     pathName="Generation of second messenger molecules"
#' )
#'
plotFoldChange <- function(
        inputList,
        pathName=NA,
        pathId=NA,
        genesToPlot=NA,
        manualTitle=NA,
        titleSize=14,
        geneFormat="ensembl",
        pCutoff=0.05,
        fcCutoff=1.5,
        cellColours=c("blue", "white", "red"),
        cellBorder=gpar(col="grey"),
        plotSignificantOnly=TRUE,
        showStars=TRUE,
        hideNonsigFC=TRUE,
        vjust=0.75,
        rot=0,
        invert=FALSE,
        log2FoldChange=FALSE,
        colSplit=NA,
        clusterRows=TRUE,
        clusterColumns=FALSE,
        colAngle=90,
        colCenter=TRUE,
        rowAngle=0,
        rowCenter=FALSE
) {

    stopifnot(
        "Provide a named list of data frames of results, with the name
        of each item in the list as the comparison name."  = {
            is.list(inputList)
            !is.null(names(inputList))
            all(unlist(lapply(inputList, is.data.frame)))
        }
    )

    stopifnot("One of 'pathName', 'pathId', 'genesToPlot' must b provided" = {
        any(!is.na(c(pathName, pathId, genesToPlot)))
    })

    stopifnot("'geneFormat' must be either 'ensembl' or 'hgnc'" = {
        geneFormat %in% c("ensembl", "hgnc")}
    )

    if (!is.na(pathName)) {
        pathId <- sigoraDatabase %>%
            filter(pathwayName == pathName) %>%
            pull(pathwayId) %>%
            unique()

        stopifnot(
            "Specified 'pathName' was not found, please try a different pathway
            name"=length(pathId) != 0
        )

        plotTitle <- pathName

    } else if (!is.na(pathId)) {

        stopifnot(
            "Specified 'pathId' was not found, please try a different pathway
            ID"=pathId %in% unique(sigoraDatabase$pathwayId)
        )

        plotTitle <- sigoraDatabase %>%
            filter(pathwayId == pathId) %>%
            pull(pathwayName) %>%
            unique()

    } else if (any(!is.na(genesToPlot))) {
        plotTitle <- NULL
    }

    ## If a title is provided manually, overwrite any other title from above
    if (!is.na(manualTitle)) {
        plotTitle <- manualTitle
    }

    ## Get the genes to plot in the pathway or gene list of interest, and make
    ## them Ensembl IDs
    if (is.na(genesToPlot[1])) {

        genes <- sigoraDatabase %>%
            filter(pathwayId == pathId) %>%
            .$ensemblGeneId

    } else {
        genes <- genesToPlot

        if (geneFormat == "hgnc") {
            genes <- mappingFile %>%
                filter(hgncSymbol %in% genesToPlot) %>%
                pull(ensemblGeneId)
        }
        if (geneFormat == "ensembl") {
            genes <- genesToPlot
        }
    }

    dfFC <- imap(inputList, function(listItem, itemName) {
        stopifnot(
            "Rownames of data frames in 'inputList' must be Ensembl gene IDs" =
            grepl(pattern = "^ENSG", x = rownames(listItem)[1])
        )

        listItem %>%
            filter(rownames(.) %in% genes, !is.na(log2FoldChange)) %>%
            rownames_to_column("ensemblGeneId") %>%
            select(ensemblGeneId, {{itemName}} := log2FoldChange)

    }) %>% reduce(full_join, by="ensemblGeneId")


    dfP <- imap(inputList, function(listItem, itemName) {
        listItem %>%
            filter(rownames(.) %in% genes, !is.na(padj)) %>%
            rownames_to_column("ensemblGeneId") %>%
            select(ensemblGeneId, {{itemName}} := padj)

    }) %>% reduce(full_join, by="ensemblGeneId")

    ## From all the data frames, get the genes that were significant in any of
    ## them
    sigGenes <- inputList %>% map(
        ~rownames_to_column(.x, "ensemblGeneId") %>%
            filter(
                ensemblGeneId %in% genes,
                padj < pCutoff,
                abs(log2FoldChange) > log2(fcCutoff)
            ) %>%
            pull(ensemblGeneId)
    ) %>%
        unlist() %>%
        unique()

    if (plotSignificantOnly) {
        dfFC <- dfFC %>% filter(ensemblGeneId %in% sigGenes)
        dfP <- dfP %>% filter(ensemblGeneId %in% sigGenes)
    }

    ## Prepare the Heatmap matrices, and map the Ensembl IDs to HGNC symbols
    matFC <- dfFC %>%
        left_join(mappingFile, by="ensemblGeneId", multiple="all") %>%
        select(-c(ensemblGeneId, entrezGeneId)) %>%
        column_to_rownames(var="hgncSymbol") %>%
        as.matrix()
    matFC[is.na(matFC)] <- 0 ## Make any NAs into 0


    matP <- dfP %>%
        left_join(mappingFile, by="ensemblGeneId", multiple="all") %>%
        select(-c(ensemblGeneId, entrezGeneId)) %>%
        column_to_rownames(var="hgncSymbol") %>%
        as.matrix()
    matP[is.na(matP)] <- 1 ## Make any NAs into 1s

    if (hideNonsigFC) {
        matFC[abs(matFC) < log2(fcCutoff)] <- 0 ## Did not pass fcCutoff
        matFC[matP > pCutoff] <- 0 ## Did not pass pCutoff
    }

    ## Use a helper function, which properly creates the colour scale and legend
    ## breaks/labels, based on:
    ## - If we're plotting log2 or standard fold changes
    ## - If the values are one sided i.e. all greater or less than 0
    heatmapLegendInfo <- .plotFoldChangeLegend(
        .matFC = matFC,
        .log2FoldChange = log2FoldChange,
        .cellColours = cellColours
    )


    ## If columns aren't being split
    if (is.na(colSplit[1])) {
        colSplit <- rep(NA, length(inputList))
        columnTitle <- NULL
    } else {
        ## This orders the splitting into the order the data frames are in the
        ## list instead of alphabetically
        colSplit <- factor(colSplit, levels=unique(colSplit))
        columnTitle <- "%s"
    }
    rowSplit <- rep(NA, nrow(matFC))
    row_title <- NULL

    ## If plotting so that row names are conditions and column names are genes
    if (invert) {
        matFC <- t(matFC)
        matP <- t(matP)
        var <- colSplit
        colSplit <- rowSplit
        rowSplit <- var
        columnTitle <- NULL
        row_title <- "%s"
    }


    draw(
        Heatmap(
            matrix = matFC,
            border = TRUE,
            col = heatmapLegendInfo[[2]],
            rect_gp = cellBorder,
            cell_fun = function(j, i, x, y, w, h, fill) {
                if (showStars) {
                    if (abs(matFC[i, j]) > log2(1.5)) {
                        if (matP[i, j] < 0.001) {
                            grid.text("***", x, y, vjust = vjust, rot = rot)
                        }
                        else if (matP[i, j] < 0.01) {
                            grid.text("**", x, y, vjust = vjust, rot = rot)
                        }
                        else if (matP[i, j] < 0.05) {
                            grid.text("*", x, y, vjust = vjust, rot = rot)
                        }
                    }

                    ## If plotting significance values for genes that don't pass
                    ## fcCutoff
                    if (abs(matFC[i, j]) < log2(1.5) & !hideNonsigFC) {
                        if (matP[i, j] < 0.001) {
                            grid.text(
                                "***",
                                x,
                                y,
                                vjust = vjust,
                                rot = rot,
                                gp = gpar(col = "grey50")
                            )
                        }
                        else if (matP[i, j] < 0.01) {
                            grid.text(
                                "**",
                                x,
                                y,
                                vjust = vjust,
                                rot = rot,
                                gp = gpar(col = "grey50")
                            )
                        }
                        else if (matP[i, j] < 0.05) {
                            grid.text(
                                "*",
                                x,
                                y,
                                vjust = vjust,
                                rot = rot,
                                gp = gpar(col = "grey50")
                            )
                        }
                    }
                }
            },
            column_title = columnTitle,
            row_title = row_title,
            heatmap_legend_param = heatmapLegendInfo[[1]],
            column_title_gp = gpar(fontsize = titleSize),
            row_split = rowSplit,
            column_split = colSplit,
            cluster_columns = clusterColumns,
            cluster_rows = clusterRows,
            column_names_rot = colAngle,
            column_names_centered = colCenter,
            row_names_rot = rowAngle,
            row_names_centered = rowCenter
        ),
        column_title = plotTitle
    )
}
