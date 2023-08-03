#' Create a volcano plot from an unfiltered DESeq2 results object
#'
#' @param deseqResults Data frame of DESeq2 results, with Ensembl gene IDs as
#'   rownames, i.e. out from `DESeq2::results()`.
#' @param pCutoff Adjusted p value cutoff, defaults to <0.05
#' @param fcCutoff Absolute fold change cutoff, defaults to >1.5
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
#' @return Volcano plot as a ggplot object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import stringr
#' @importFrom ggrepel geom_text_repel
#'
#' @description Creates a volcano plot of genes output from `DESeq2::results()`,
#'   with various options for tweaking the appearance. Ensembl gene IDs should
#'   be the rownames of the input data frame.
#'
#'   The argument `highlightGenes` can be used to draw attention to a specific
#'   set of genes, e.g. those from a pathway of interest. Setting the argument
#'   `label="highlight"` will also mean those same genes (at least some of them)
#'   will be given labels, further emphasizing them in the volcano plot.
#'
#'   Since this function returns a ggplot object, further custom changes could
#'   be applied using the standard ggplot2 functions.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' eruption(deseqResults=deseqExampleList[[1]])
#'
eruption <- function(
        deseqResults,
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

    ## Input checks
    stopifnot(is(deseqResults, "data.frame"))
    stopifnot(all(c("padj", "log2FoldChange") %in% colnames(deseqResults)))

    ## If Ensembl gene IDs are detected, annotate them with gene names from the
    ## mapping file. For Ensembl gene IDs without gene names, just use the
    ## Ensembl gene ID. If rownames are not ENSG ids, they will be used as is,
    ## and you can map it to your own ids beforehand
    if (str_detect(rownames(deseqResults)[1], "^ENSG")) {
        res <- deseqResults %>%
            rownames_to_column("ensemblGeneId") %>%
            filter(!is.na(padj)) %>%
            left_join(
                mappingFile,
                by="ensemblGeneId",
                multiple="all"
            ) %>%
            mutate(
                geneName=ifelse(
                    !is.na(hgncSymbol),
                    hgncSymbol,
                    ensemblGeneId
                )
            )
    } else {
        res <- deseqResults %>%
            rownames_to_column("ensemblGeneId") %>%
            mutate(geneName=ensemblGeneId)
    }

    res <- res %>% mutate(
        significant=case_when(
            padj < pCutoff & abs(log2FoldChange) > log2(fcCutoff) ~ "SIG",
            TRUE ~ "NS"
        ),
        inList=case_when(ensemblGeneId %in% highlightGenes ~ "Y", TRUE ~ "N"),
        negLogP=-log10(padj)
    )

    ## If specifying x and y axis limits, remove any genes that fall outside the
    ## specified ranges
    if (!is.na(xaxis[1])) {
        res <- filter(
            res,
            log2FoldChange > xaxis[1] & log2FoldChange < xaxis[2]
        )
    }
    if (!is.na(yaxis[1])) {
        res <- filter(
            res,
            negLogP > yaxis[1] & negLogP < yaxis[2]
        )
    }

    ## Make sure highlightGenes, if used, are in the data frame
    highlightGenes <- highlightGenes[highlightGenes %in% res$ensemblGeneId]

    ## Identify top up- and down-regulated genes, and output them for the user
    upDf <- res %>%
        filter(padj < pCutoff, log2FoldChange > log2(fcCutoff))

    downDf <- res %>%
        filter(padj < pCutoff, log2FoldChange < log2(1 / fcCutoff))

    ## Number of up, down, and total DE genes
    numGenes <- c(nrow(upDf), nrow(downDf), nrow(upDf) + nrow(downDf))

    message(
        "Creating volcano plot with ",
        numGenes[3], " DEGs: ",
        numGenes[2], " down-regulated and ",
        numGenes[1], " up-regulated."
    )

    if (length(highlightGenes) > 0) {
        message(
            "An additional ",
            length(highlightGenes),
            " genes will be highlighted."
        )
    }

    ## Select the genes to label. You can specify if you want to only label
    ## annotated genes (default), i.e. those that have a gene name. This is for
    ## "auto" and "highlight" labeling, i.e. it does not apply to manual
    ## labeling (e.g. if you want to manually label an unannotated gene).
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
            arrange(desc(log2FoldChange ^ 2 * log10(padj) ^ 2)) %>%
            filter(geneName %in% possibleLabels) %>%
            head(n) %>%
            select(geneName) %>%
            unlist()

        downGenes <- downDf %>%
            arrange(desc(log2FoldChange ^ 2 * log10(padj) ^ 2)) %>%
            filter(geneName %in% possibleLabels) %>%
            head(n) %>%
            select(geneName) %>%
            unlist()

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
            arrange(desc(log2FoldChange ^ 2 * log10(padj) ^ 2)) %>%
            head(n) %>%
            select(geneName) %>%
            unlist()

        downGenes <- downDf %>%
            filter(inList == "Y", geneName %in% possibleLabels) %>%
            arrange(desc(log2FoldChange ^ 2 * log10(padj) ^ 2)) %>%
            head(n) %>%
            select(geneName) %>%
            unlist()

        res <- res %>% mutate(
            label=case_when(geneName %in% c(upGenes, downGenes) ~ geneName)
        )
    }

    ## Manual labeling ("manual"): label the genes you provided in `manualGenes`
    if (label == "manual") {
        res <- res %>% mutate(label=case_when(
            ensemblGeneId %in% manualGenes |
                geneName %in% manualGenes ~ geneName
        ))
    }

    ## Create the plot
    p <- ggplot(res, aes(x=log2FoldChange, y=negLogP)) +

        ## Plot the non-significant genes
        geom_point(
            data=filter(res, significant == "NS", inList == "N"),
            alpha=alpha,
            show.legend=FALSE,
            size=pointSize,
            colour=nonsigColour
        ) +

        ## Plot the significant genes and genes of interest ("highlightGenes"),
        ## with those in `highlightGenes` plotted on top of those not in
        ## `highlightGenes` for added emphasis.
        geom_point(
            data=filter(res, significant == "SIG", inList == "N"),
            mapping=aes(x=log2FoldChange, y=negLogP),
            size=pointSize,
            alpha=alpha,
            colour=baseColour
        ) +

        geom_point(
            data=filter(res, inList == "Y"),
            mapping=aes(x=log2FoldChange, y=negLogP),
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
            axis.title=element_text(
                size=13, face="bold", colour="black"
            )
        ) +

        labs(
            x=expression(bold(log["2"]~Fold~Change)),
            y=expression(bold(-log["10"]~P[adj]))
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
                min(res$log2FoldChange) - 0.1,
                max(res$log2FoldChange) + 0.1
            )
        }

        ## All the work for setting the x axis breaks/labels is done inside
        ## `.eruptionBreaks` to shorten and simplify this script
        p <- p +
            xlab("Fold Change") +
            theme(axis.text.x=element_text(vjust=0, size=11)) +
            .eruptionBreaks(xaxis)
    }

    return(p)
}
