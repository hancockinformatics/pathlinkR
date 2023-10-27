#' Create a volcano plot of RNA-Seq results
#'
#' @param rnaseqResult Data frame of RNASeq results, with Ensembl gene IDs as
#'   rownames. Can be a "DESeqResults" or "TopTags" object, or a simple data
#'   frame. See "Details" for more information.
#' @param columnFC Character; Column to plot along the x-axis, typically log2
#'   fold change values. Only required when `rnaseqResult` is a simple data
#'   frame. Defaults to NA.
#' @param columnP Character; Column to plot along the y-axis, typically nominal
#'   or adjusted p values. Only required when `rnaseqResult` is a simple data
#'   frame. Defaults to NA.
#' @param pCutoff Adjusted p value cutoff, defaults to < 0.05
#' @param fcCutoff Absolute fold change cutoff, defaults to > 1.5
#' @param baseColour Colour of points for all significant DE genes
#'   ("steelblue4")
#' @param nonsigColour Colour of non-significant DE genes ("lightgrey")
#' @param alpha Transparency of the points (0.5)
#' @param pointSize Size of the points (1)
#' @param title Title of the plot
#' @param nonlog2 Show non-log2 fold changes instead of log2 fold change (FALSE)
#' @param xaxis Length-two numeric vector to manually specify limits of the
#'   x-axis in log2 fold change; defaults to NA which lets ggplot2 determine the
#'   best values.
#' @param yaxis Length-two numeric vector to manually specify limits of the
#'   y-axis (in -log10). Defaults to NA which lets ggplot2 determine the best
#'   values.
#' @param highlightGenes Vector of genes to emphasize by colouring differently
#'   (e.g. genes of interest). Must be Ensembl IDs.
#' @param highlightColour Colour for the genes specified in `highlightGenes`
#' @param highlightName Optional name to call the `highlightGenes` (e.g.
#'   Unique, Shared, Immune related, etc.)
#' @param label When set to "auto" (default), label the top `n` up- and
#'   down-regulated DE genes. When set to "highlight", label top `n` up- and
#'   down-regulated genes provided in `highlightGenes`. When set to "manual"
#'   label a custom selection of genes provided in `manualGenes`.
#' @param n number of top up- and down-regulated genes to label. Applies when
#'   `label` is set to "auto" or "highlight".
#' @param manualGenes If `label="manual"`, these are the genes to be
#'   specifically label. Can be HGNC symbols or Ensembl gene IDs.
#' @param removeUnannotated Boolean (TRUE): Remove genes without annotations
#'   (no HGNC symbol).
#' @param labelSize Size of font for labels
#' @param pad Padding of labels; adjust this if the labels overlap
#'
#' @return Volcano plot of genes from an RNA-Seq experiment; a "ggplot" object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom tibble rownames_to_column
#'
#' @description Creates a volcano plot of genes from RNA-Seq results, with
#'   various options for tweaking the appearance. Ensembl gene IDs should be the
#'   rownames of the input object.
#'
#' @details The input to `eruption()` can be of class "DESeqResults" (from
#'   `DESeq2`), "TopTags" (`edgeR`), or a simple data frame. When providing
#'   either of the former, the columns to plot are automatically pulled
#'   ("log2FoldChange" and "padj" for DESeqResults, or "logFC" and "FDR" for
#'   TopTags). Otherwise, the arguments "columnFC" and "columnP" must be
#'   specified. If one wishes to override the default behaviour for
#'   "DESeqResults" or "TopTags" (e.g. plot nominal p values on the y-axis),
#'   convert those objects to data frames, then supply "columnFC" and "columnP".
#'
#'   The argument `highlightGenes` can be used to draw attention to a specific
#'   set of genes, e.g. those from a pathway of interest. Setting the argument
#'   `label="highlight"` will also mean those same genes (at least some of them)
#'   will be given labels, further emphasizing them in the volcano plot.
#'
#'   Since this function returns a ggplot object, further custom changes could
#'   be applied using the standard ggplot2 functions (`labs()`,
#'   `theme()`, etc.).
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' data("exampleDESeqResults")
#' eruption(rnaseqResult=exampleDESeqResults[[1]])
#'
eruption <- function(
        rnaseqResult,
        columnFC=NA,
        columnP=NA,
        pCutoff=0.05,
        fcCutoff=1.5,
        baseColour="steelblue4",
        nonsigColour="lightgrey",
        alpha=0.5,
        pointSize=1,
        title=NA,
        nonlog2=FALSE,
        xaxis=NA,
        yaxis=NA,
        highlightGenes=c(),
        highlightColour="red",
        highlightName="Selected",
        label="auto",
        n=10,
        manualGenes=c(),
        removeUnannotated=TRUE,
        labelSize=3.5,
        pad=1.4
) {

    if (is(rnaseqResult, "DESeqResults")) {
        rnaseqResult <- rnaseqResult %>%
            as.data.frame() %>%
            rename("LogFoldChange"=log2FoldChange, "PAdjusted"=padj)

    } else if (is(rnaseqResult, "TopTags")) {
        rnaseqResult <- rnaseqResult %>%
            as.data.frame() %>%
            rename("LogFoldChange"=logFC, "PAdjusted"=FDR)

    } else {

        stopifnot(
            "If 'rnaseqResult' is a data frame, you must provide
            'columnFC' and 'columnP'"={
                !any(is.na(columnFC), is.na(columnP))
            }
        )

        rnaseqResult <- rnaseqResult %>%
            rename(
                "LogFoldChange"=all_of(columnFC),
                "PAdjusted"=all_of(columnP)
            )
    }

    if (!all(is.na(xaxis))) {
        stopifnot("'xaxis' must be a length-two numeric vector"={
            length(xaxis) == 2
        })
    }
    if (!all(is.na(yaxis))) {
        stopifnot("'yaxis' must be a length-two numeric vector"={
            length(yaxis) == 2
        })
    }


    ## Load the data we need for gene mapping
    data_env <- new.env(parent=emptyenv())
    data("mappingFile", envir=data_env, package="pathlinkR")
    mappingFile <- data_env[["mappingFile"]]


    ## If Ensembl IDs are detected, annotate them with gene names from the
    ## mapping file. If rownames are not Ensembl IDs, they will be used as-is.
    if (grepl(x=rownames(rnaseqResult)[1], pattern="^ENSG")) {
        res <- rnaseqResult %>%
            rownames_to_column("ensemblGeneId") %>%
            filter(!is.na(PAdjusted)) %>%
            left_join(mappingFile, by="ensemblGeneId", multiple="all") %>%
            mutate(geneName=ifelse(
                !is.na(hgncSymbol),
                hgncSymbol,
                ensemblGeneId
            ))
    } else {
        res <- rnaseqResult %>%
            rownames_to_column("ensemblGeneId") %>%
            mutate(geneName=ensemblGeneId)
    }

    res <- res %>% mutate(
        significant=case_when(
            PAdjusted < pCutoff & abs(LogFoldChange) > log2(fcCutoff) ~ "SIG",
            TRUE ~ "NS"
        ),
        inList=case_when(ensemblGeneId %in% highlightGenes ~ "Y", TRUE ~ "N"),
        negLogP=-log10(PAdjusted)
    )

    ## If specifying x and y axis limits, remove any genes that fall outside the
    ## specified ranges
    if (!is.na(xaxis[1])) {
        res <- filter(
            res,
            LogFoldChange > xaxis[1] & LogFoldChange < xaxis[2]
        )
    }
    if (!is.na(yaxis[1])) {
        res <- filter(
            res,
            negLogP > yaxis[1] & negLogP < yaxis[2]
        )
    }

    ## Make sure `highlightGenes` are in the input
    highlightGenes <- highlightGenes[highlightGenes %in% res$ensemblGeneId]

    ## Identify up- and down-regulated genes
    upDf <- res %>%
        filter(PAdjusted < pCutoff, LogFoldChange > log2(fcCutoff))

    downDf <- res %>%
        filter(PAdjusted < pCutoff, LogFoldChange < log2(1 / fcCutoff))

    numGenes <- c(nrow(upDf), nrow(downDf), nrow(upDf) + nrow(downDf))

    if (removeUnannotated) {
        possibleLabels <- res %>%
            filter(!grepl("ENSG", geneName)) %>%
            pull(geneName)
    } else {
        possibleLabels <- res$geneName
    }

    ## Auto-labeling (default): label the top 5 (or n) up and down-regulated
    ## genes. The "top" genes are those with the highest fold change and
    ## smallest p-value.
    if (label == "auto") {
        upGenes <- upDf %>%
            arrange(desc(LogFoldChange ^ 2 * log10(PAdjusted) ^ 2)) %>%
            filter(geneName %in% possibleLabels) %>%
            head(n) %>%
            pull(geneName)

        downGenes <- downDf %>%
            arrange(desc(LogFoldChange ^ 2 * log10(PAdjusted) ^ 2)) %>%
            filter(geneName %in% possibleLabels) %>%
            head(n) %>%
            pull(geneName)

        ## Record which genes to label
        res <- res %>% mutate(
            label=case_when(geneName %in% c(upGenes, downGenes) ~ geneName)
        )
    }

    ## Selective labeling ("highlight"): label the top n up/down genes that are
    ## in `highlightGenes`
    if (label == "highlight") {
        upGenes <- upDf %>%
            filter(inList == "Y", geneName %in% possibleLabels) %>%
            arrange(desc(LogFoldChange ^ 2 * log10(PAdjusted) ^ 2)) %>%
            head(n) %>%
            pull(geneName)

        downGenes <- downDf %>%
            filter(inList == "Y", geneName %in% possibleLabels) %>%
            arrange(desc(LogFoldChange ^ 2 * log10(PAdjusted) ^ 2)) %>%
            head(n) %>%
            pull(geneName)

        res <- res %>% mutate(
            label=case_when(geneName %in% c(upGenes, downGenes) ~ geneName)
        )
    }

    ## Manual labeling ("manual"): label the genes you provided in "manualGenes"
    if (label == "manual") {
        res <- res %>% mutate(label=case_when(
            ensemblGeneId %in% manualGenes |
                geneName %in% manualGenes ~ geneName
        ))
    }

    ## Create the plot
    p <- ggplot(res, aes(x=LogFoldChange, y=negLogP)) +

        ## Plot the non-significant genes
        geom_point(
            data=filter(res, significant == "NS", inList == "N"),
            alpha=alpha,
            show.legend=FALSE,
            size=pointSize,
            colour=nonsigColour
        ) +

        ## Plot the significant genes and genes of interest, with those in
        ## "highlightGenes" plotted on top of the others for added emphasis
        geom_point(
            data=filter(res, significant == "SIG", inList == "N"),
            mapping=aes(x=LogFoldChange, y=negLogP),
            size=pointSize,
            alpha=alpha,
            colour=baseColour
        ) +

        geom_point(
            data=filter(res, inList == "Y"),
            mapping=aes(x=LogFoldChange, y=negLogP),
            size=pointSize,
            alpha=alpha,
            colour=highlightColour
        ) +

        ## Add cutoff lines
        geom_hline(
            yintercept=-log10(pCutoff),
            linetype="dashed",
            colour="gray20"
        ) +
        geom_vline(
            xintercept=c(log2(fcCutoff), -log2(fcCutoff)),
            linetype="dashed",
            colour="gray20"
        ) +

        ## Set to a clean theme, and make the labels easier to read
        theme_bw() +
        theme(
            plot.title=element_text(face="bold", size=18),
            plot.subtitle=element_text(size=14),
            plot.background=element_blank(),
            axis.text=element_text(colour="black", size=11),
            axis.title=element_text(size=13, colour="black")
        ) +

        labs(
            x=expression(log[2]~Fold~change),
            y=expression(-~log[10]~P[adjusted])
        ) +

        ## Set axes
        {if (!is.na(yaxis[1])) ylim(yaxis)} +
        {if (!is.na(xaxis[1])) xlim(xaxis)} +

        ## Add in labels to the genes
        geom_text_repel(
            data=res %>% filter(!is.na(label)),
            mapping=aes(label=label),
            segment.color="darkgrey",
            color="black",
            bg.color="white",
            bg.r=.05,
            max.overlaps=Inf,
            fontface="bold",
            size=labelSize,
            box.padding=pad
        ) +

        ## Add in informative subtitles for number of up and down-regulated
        ## genes, the number of genes in "highlightGenes," and a title if
        ## provided
        {if (!is.na(title)) labs(title=title)} +

        ## If there are no "highlightGenes"
        {
            if (length(highlightGenes) == 0)
                labs(subtitle=paste0(
                    "Down: ", numGenes[2], ", Up: ", numGenes[1]
                ))
        } +

        ## If "highlightGenes" was given
        {
            if (length(highlightGenes) != 0)
                labs(subtitle=paste0(
                    "Down: ", numGenes[2],
                    ", Up: ", numGenes[1], ", ",
                    highlightName, ": ", length(highlightGenes)
                ))
        }

    ## For plotting with fold change, not log2 fold change, manually add the
    ## scientific notation based on the x axis range
    if (nonlog2) {

        ## If "xaxis" is not specified
        if (is.na(xaxis[1])) {
            xaxis <- c(
                min(res$LogFoldChange) - 0.1,
                max(res$LogFoldChange) + 0.1
            )
        }

        ## All the work for setting the x axis breaks/labels is done inside
        ## `.eruptionBreaks` to shorten and simplify this script
        p <- p +
            labs(x="Fold Change") +
            theme(axis.text.x=element_text(vjust=0, size=11)) +
            .eruptionBreaks(xaxis)
    }

    return(p)
}


#' INTERNAL Create manual breaks/labels for volcano plots
#'
#' @param x Length-two numeric vector to manually specify limits of the
#'   x-axis in log2 fold change; defaults to NA which lets ggplot2 determine the
#'   best values.
#'
#' @return ggplot scale object
#'
#' @import dplyr
#' @import ggplot2
#'
#' @description Internal function which is used to create even breaks for
#'   volcano plots produced by `eruption`.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
.eruptionBreaks <- function(x) {

    ## Up to log2FC=6
    if (between(max(abs(x)), 0, 6)) {
        scale_x_continuous(
            breaks=c(-6, -4, -2, 0, 2, 4, 6),
            labels=c(-32, -16, -4, 1, 4, 16, 32),
            limits=x
        )

        ## Up to log2FC=6.64
    } else if (between(max(abs(x)), 6.001, log2(100))) {
        scale_x_continuous(
            breaks=c(-log2(100), -log2(10), 0, log2(10), log2(100)),
            labels=c(
                expression(-10 ^ 2),
                expression(-10 ^ 1),
                1,
                expression(10 ^ 1),
                expression(10 ^ 2)
            ),
            limits=x
        )

        ## Up to log2FC=9.97
    } else if (between(max(abs(x)), log2(100), log2(1000))) {
        scale_x_continuous(
            breaks=c(
                -log2(1000),
                -log2(100),
                -log2(10),
                0,
                log2(10),
                log2(100),
                log2(1000)
            ),
            labels=c(
                expression(-10 ^ 3),
                expression(-10 ^ 2),
                expression(-10 ^ 1),
                1,
                expression(10 ^ 1),
                expression(10 ^ 2),
                expression(10 ^ 3)
            ),
            limits=x
        )

        ## Up to log2FC=13.29
    } else if (between(max(abs(x)), log2(1000), log2(10000))) {
        scale_x_continuous(
            breaks=c(
                -log2(10000),
                -log2(1000),
                -log2(100),
                -log2(10),
                0,
                log2(10),
                log2(100),
                log2(1000),
                log2(10000)
            ),
            labels=c(
                expression(-10 ^ 4),
                expression(-10 ^ 3),
                expression(-10 ^ 2),
                expression(-10 ^ 1),
                1,
                expression(10 ^ 1),
                expression(10 ^ 2),
                expression(10 ^ 3),
                expression(10 ^ 4)
            ),
            limits=x
        )

        ## Up to log2FC=16.61
    } else if (between(max(abs(x)), log2(10000), log2(100000))) {
        scale_x_continuous(
            breaks=c(
                -log2(10 ^ 5),
                -log2(10 ^ 3),
                -log2(10),
                0,
                log2(10),
                log2(10 ^ 3),
                log2(10 ^ 5)
            ),
            labels=c(
                expression(-10 ^ 5),
                expression(-10 ^ 3),
                expression(-10 ^ 1),
                1,
                expression(10 ^ 1),
                expression(10 ^ 3),
                expression(10 ^ 5)
            ),
            limits=x
        )

    } else if (max(abs(x)) >= log2(100000)) {
        message(
            "Something may be wrong with your DESeq model to have fold ",
            "changes > 10^5..."
        )
    }
}
