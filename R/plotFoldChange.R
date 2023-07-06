#' Create fold change plots to visualize DESeq2 results
#'
#' @param inputList List of data frames of DESeq2 results. The list names are
#'   used as the comparison for each dataframe (e.g. COVID vs Healthy). Data
#'   frames should have Ensembl gene IDs as rownames.
#' @param pathName Name of pathway to pull genes from, will be plot title
#' @param pathId ID of pathway to pull genes from, if pathName is not given
#' @param manualTitle Set this as your title instead of pathway name
#' @param titleSize Font size for title
#' @param genesToPlot Set if you want a specific vector of genes instead of
#'   pulling from a pathway
#' @param geneFormat Default is Ensembl gene IDs ("ensg"), can also input a
#'   vector of HGNC symbols ("hgnc")
#' @param pCutoff P value cutoff. Default is p < 0.05
#' @param fcCutoff Fold change cutoff. Default is |FC| > 1.5
#' @param plotSignificantOnly Boolean (TRUE) Only plot genes that are
#'   differentially expressed (pass pCutoff and fcCutoff) in any comparison
#' @param showStars Boolean (TRUE) show significance stars on heatmap
#' @param hideNonsigFC Boolean (TRUE) If a gene is significant in one
#'   comparison but not in another, this will set the colour of the non-
#'   significant gene as grey to visually emphasize the significant genes. If
#'   set to FALSE, it will set the colour to the fold change, and if the p
#'   value passes pCutoff, it will also display the p value (the asterisks
#'   will be grey instead of black).
#' @param vjust Adjustment of the position of the significance stars. Default
#'   is 0.75. May need to adjust if there are many genes
#' @param rot Rotation of the position of the significance stars. Default is 0
#' @param invert Boolean (FALSE) Default plots genes as rows and conditions as
#'   columns, set to TRUE if you want genes as columns and conditions as rows
#' @param log2FoldChange Boolean (FALSE) Default plots the fold changes in the
#'   legend as the true fold change, set to TRUE if you want log2 fold change
#' @param colSplit Split each condition into groups for better visualization.
#'   To do so, create a vector where each item of the vector corresponds to
#'   which group the dataframe belongs to in inputList. E.g. ('Positive',
#'   'Positive', "Negative", "Negative", "Time", "Time"). The order of these
#'   groups is set in the order they appear; if you want to change the order,
#'   change the order of data frames in inputList and write colSplit in the
#'   order you want. Order will be ignored if clusterColumns is set to TRUE
#' @param clusterRows Boolean (TRUE) Whether to cluster the rows (genes). May
#'   need to change if invert = TRUE.
#' @param clusterColumns Boolean (FALSE) Whether to cluster the columns
#'   (conditions). Will override order of colSplit if set to TRUE.
#' @param colAngle angle of column text. Set default to 90
#' @param colCenter whether to center column text. Default is TRUE, should set
#'   to FALSE if angled column name (e.g. colAngle = 45)
#' @param rowAngle angle of row text. Set default to 0
#' @param rowCenter whether to center column text. Default is FALSE, should
#'   set to TRUE if vertical column name (e.g. rowAngle = 90).
#'
#'
#' @return A heatmap of fold changes for genes of interest
#' @export
#'
#' @import ComplexHeatmap
#' @import dplyr
#'
#' @description Creates a heatmap of fold changes values for results from the
#'   `DESeq2::results()` function, with various parameters to tweak the
#'   appearance.
#'
#' @references <https://bioconductor.org/packages/ComplexHeatmap/>
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' plotFoldChange(
#'     deseqExampleList,
#'     pathName = "Generation of second messenger molecules"
#' )
#'
plotFoldChange <- function(
        inputList,
        pathName = NA,
        pathId = NA,
        manualTitle = NA,
        titleSize = 14,
        genesToPlot = NA,
        geneFormat = "ensg",
        pCutoff = 0.05,
        fcCutoff = 1.5,
        plotSignificantOnly = TRUE,
        showStars = TRUE,
        hideNonsigFC = TRUE,
        vjust = 0.75,
        rot = 0,
        invert = FALSE,
        log2FoldChange = FALSE,
        colSplit = NA,
        clusterRows = TRUE,
        clusterColumns = FALSE,
        colAngle = 90,
        colCenter = TRUE,
        rowAngle = 0,
        rowCenter = FALSE) {

    # If pathway name is provided
    if (!is.na(pathName)) {
        pathId <- sigoraDatabase %>%
            filter(pathwayName == pathName) %>%
            select(pathwayId) %>%
            unlist() %>%
            .[1] %>%
            as.character()
        plotTitle <- pathName
    }

    # If pathway ID is provided
    if (!is.na(pathId)) {
        plotTitle <- sigoraDatabase %>%
            filter(pathwayId == pathId) %>%
            select(pathwayName) %>%
            unlist() %>%
            .[1] %>%
            as.character()
    }

    # If a title is provided manually, overwrite
    if(!is.na(manualTitle)){
        plotTitle <- manualTitle
    }

    # Get the genes to plot in the pathway or gene list of interest, and make
    # them Ensembl IDs
    if (is.na(genesToPlot[1])) {

        genes <- sigoraDatabase %>%
            filter(pathwayId == pathId) %>%
            .$ensemblGeneId

    } else {
        genes <- genesToPlot

        if (geneFormat == "hgnc") {
            genes <- mappingFile %>%
                filter(hgncSymbol %in% genesToPlot) %>%
                .$ensemblGeneId
        }
        if (geneFormat == "ensg") {
            genes <- genesToPlot
        }
    }

    dfFC <- imap(inputList, function(listItem, itemName) {
        listItem %>%
            filter(rownames(.) %in% genes, !is.na(log2FoldChange)) %>%
            rownames_to_column("ensemblGeneId") %>%
            select(ensemblGeneId, {{itemName}} := log2FoldChange)
    }) %>% plyr::join_all(by = "ensemblGeneId", type = "full")

    dfP <- imap(inputList, function(listItem, itemName) {
        listItem %>%
            filter(rownames(.) %in% genes, !is.na(padj)) %>%
            rownames_to_column("ensemblGeneId") %>%
            select(ensemblGeneId, {{itemName}} := padj)
    }) %>% plyr::join_all(by = "ensemblGeneId", type = "full")

    # From all the data frames, get the genes that were significant in any of
    # them
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

    # Prepare the Heatmap matrices, and map the Ensembl IDs to HGNC symbols
    matFC <- dfFC %>%
        left_join(mappingFile, by = "ensemblGeneId", multiple = "all") %>%
        select(-c(ensemblGeneId, entrezGeneId)) %>%
        column_to_rownames(var = "hgncSymbol") %>%
        as.matrix()
    matFC[is.na(matFC)] <- 0 # make any NAs into 0


    matP <- dfP %>%
        left_join(mappingFile, by = "ensemblGeneId", multiple = "all") %>%
        select(-c(ensemblGeneId, entrezGeneId)) %>%
        column_to_rownames(var = "hgncSymbol") %>%
        as.matrix()
    matP[is.na(matP)] <- 1 # make any NAs into 1s

    if (hideNonsigFC) {
        matFC[abs(matFC) < log2(fcCutoff)] <- 0 # did not pass fc cutoff
        matFC[matP > pCutoff] <- 0 # did not pass pval cutoff
    }

    # Set the limits for colours for the plotting heatmap
    limit <- ceiling(max(abs(matFC), na.rm = TRUE))

    # If plotting real fold changes instead of log2
    if(!log2FoldChange){
        if ((limit %% 2) == 0) {
            range <- c(-limit, -limit / 2, 0, limit / 2, limit)
        } else {
            limit <- limit + 1
            range <- c(-limit, -limit / 2, 0, limit / 2, limit)
        }
        foldChangeTitle <- "Fold\nChange"
        labels <- c(-2^limit, -2^(limit/2), 1, 2^(limit/2), 2^limit)
    } else {
        range <- c(-limit, -limit / 2, 0, limit / 2, limit)
        foldChangeTitle <- "log2 FC"
        labels <- range
    }
    parameters <- list(
        at = range,
        labels = labels,
        title = foldChangeTitle
    )

    # If columns aren't being split
    if (is.na(colSplit[1])) {
        colSplit <- rep(NA, length(inputList))
        columnTitle <- NULL
    } else {
        # this orders the splitting into the order the data frames are in the
        # list instead of alphabetically
        colSplit <- factor(colSplit, levels = unique(colSplit))
        columnTitle <- "%s"
    }
    rowSplit <- rep(NA, nrow(matFC))
    row_title <- NULL

    # If plotting so that row names are conditions and column names are genes
    if (invert) {
        matFC <- t(matFC)
        matP <- t(matP)
        var <- colSplit
        colSplit <- rowSplit
        rowSplit <- var
        columnTitle <- NULL
        row_title <- "%s"
    }

    draw(Heatmap(
        matFC,
        cell_fun = function(j, i, x, y, w, h, fill) {
            if (showStars) {
                if (abs(matFC[i,j]) > log2(1.5)) {
                    if (matP[i, j] < 0.001) {
                        grid::grid.text("***", x, y, vjust = vjust, rot = rot)
                    }
                    else if (matP[i, j] < 0.01) {
                        grid::grid.text("**", x, y, vjust = vjust, rot = rot)
                    }
                    else if (matP[i, j] < 0.05) {
                        grid::grid.text("*", x, y, vjust = vjust, rot = rot)
                    }
                    # as.character(expression('\u2736') # this doesn't work?
                }

                # If plotting significance values for genes that don't pass
                # fcCutoff
                if (abs(matFC[i,j]) < log2(1.5) & !hideNonsigFC) {
                    if (matP[i, j] < 0.001) {
                        grid::grid.text(
                            "***",
                            x,
                            y,
                            vjust = vjust,
                            rot = rot,
                            gp = grid::gpar(col = "grey50")
                        )
                    }
                    else if (matP[i, j] < 0.01) {
                        grid::grid.text(
                            "**",
                            x,
                            y,
                            vjust = vjust,
                            rot = rot,
                            gp = grid::gpar(col = "grey50")
                        )
                    }
                    else if (matP[i, j] < 0.05) {
                        grid::grid.text(
                            "*",
                            x,
                            y,
                            vjust = vjust,
                            rot = rot,
                            gp = grid::gpar(col = "grey50")
                        )
                    }
                }
            }
        },
        column_title = columnTitle,
        row_title = row_title,
        heatmap_legend_param = parameters,
        column_title_gp = grid::gpar(fontsize = titleSize),
        row_split = rowSplit,
        column_split = colSplit,
        cluster_columns = clusterColumns,
        cluster_rows = clusterRows,
        column_names_rot = colAngle,
        column_names_centered = colCenter,
        row_names_rot = rowAngle,
        row_names_centered = rowCenter

    ), column_title = plotTitle)
}