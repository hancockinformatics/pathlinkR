#' Create a volcano plot from DESeq2 results
#'
#' @param deseq_results Data frame of DESeq2 results with Ensembl gene IDs as
#'   rownames
#' @param p_cutoff Adjusted p value cutoff, defaults to <0.05
#' @param fc_cutoff Absolute fold change cutoff, defaults to an absolute value
#'   of >1.5
#' @param base_colour Colour of points for significant DE genes ("steelblue4")
#' @param nonsig_colour Colour of non-DE genes ("lightgrey")
#' @param alpha Transparency of the points (0.5)
#' @param point_size Size of the points (1)
#' @param title Title of the plot
#' @param nonlog2 Show non-log2 fold change instead of log2 fold change (FALSE)
#' @param xaxis Length-two numeric vector to manually specify limits of the
#'   x-axis in log2 fold change; defaults to NA which lets ggplot2 determine the
#'   best values.
#' @param yaxis Length-two numeric vector to manually specify limits of the
#'   y-axis (in -log10); defaults to NA which lets ggplot2 determine the best
#'   values.
#' @param select_genes Vector of genes to emphasize by colouring differently
#'   (e.g. genes of interest); should be Ensembl IDs.
#' @param select_colour Colour for the `select_genes`
#' @param select_name Optional name to call the `select_genes` (e.g. Unique,
#'   Shared, Immune Related, etc.)
#' @param label When set to "auto" (default), label the top n up- and
#'   down-regulated DE genes. When set to "select", label top n up- and
#'   down-regulated genes provided in `select_genes`. When set to "manual" label
#'   a custom selection of genes provided in `manual_genes`
#' @param manual_genes If `label = "manual"`, these are the genes to
#'   specifically label. Can be HGNC symbols or Ensembl gene IDs.
#' @param remove_unannotated Boolean: Remove genes without annotations (no HGNC
#'   symbol). Defaults to TRUE.
#' @param n number of top up- and down-regulated genes to label, applies when
#'   `label` is set to "auto" or "select".
#' @param label_size Size of font for labels
#' @param pad Padding of labels; adjust this if the labels overlap
#'
#' @return Volcano plot as a ggplot object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#'
#' @description Creates a volcano plot of genes output from `DESeq2::results()`,
#'   with various options for tweaking the appearance. Ensembl gene IDs should
#'   be the rownames of the input data frame.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' eruption(deseq_results = deseq_example_list[[1]])
#'
eruption <- function(
        deseq_results,
        p_cutoff = 0.05,
        fc_cutoff = 1.5,
        base_colour = "steelblue4",
        nonsig_colour = "lightgrey",
        alpha = 0.5,
        point_size = 1,
        title = NA,
        nonlog2 = FALSE,
        xaxis = NA,
        yaxis = NA,
        select_genes = c(),
        select_colour = "red",
        select_name = "Selected",
        label = "auto",
        manual_genes = c(),
        remove_unannotated = TRUE,
        n = 10,
        label_size = 3.5,
        pad = 1.4
) {

    # Input checks
    # Data frame
    stopifnot(is.data.frame(deseq_results))

    # Ensembl IDs as rownames
    stopifnot(str_detect(rownames(deseq_results)[1], "^ENSG"))

    # Columns: padj, log2FoldChange
    stopifnot(all(c("padj", "log2FoldChange") %in% colnames(deseq_results)))

    # Annotate Ensembl gene IDs with gene names from the mapping file. For
    # Ensembl gene IDs without gene names, just use the Ensembl gene ID.
    res <- deseq_results %>%
        rownames_to_column("ensg_id") %>%
        filter(!is.na(padj)) %>%
        left_join(
            mapping_file,
            by = "ensg_id"
        ) %>%
        mutate(
            gene_name = case_when(
                is.na(gene_name) ~ ensg_id,
                !is.na(gene_name) ~ gene_name
            ),
            significant = case_when(
                padj < p_cutoff & abs(log2FoldChange) > log2(fc_cutoff) ~ "SIG",
                TRUE ~ "NS"
            ),
            in_list = case_when(ensg_id %in% select_genes ~ "Y", TRUE ~ "N"),
            neglogp = -log10(padj)
        )

    # If specifying x and y axis limits, remove any genes that fall outside the
    # specified ranges
    if(!is.na(xaxis[1])) {
        res <- res %>%
            filter(log2FoldChange > xaxis[1] & log2FoldChange < xaxis[2])
    }
    if(!is.na(yaxis[1])) {
        res <- res %>% filter(neglogp > yaxis[1] & neglogp < yaxis[2])
    }

    # Make sure select_genes, if used, are in the data frame
    select_genes <- select_genes[select_genes %in% res$ensg_id]

    # Identify top up- and down-regulated genes, and output them for the user
    up_df <- res %>%
        filter(padj < p_cutoff, log2FoldChange > log2(fc_cutoff))

    down_df <- res %>%
        filter(padj < p_cutoff, log2FoldChange < log2(1/fc_cutoff))

    # number of up, down, and total DE genes
    num_genes <- c(nrow(up_df), nrow(down_df), nrow(up_df) + nrow(down_df))

    message(
        "Creating volcano plot with ",
        num_genes[3], " DEGs: ",
        num_genes[2], " down-regulated and ",
        num_genes[1], " up-regulated."
    )

    if (length(select_genes) > 0) {
        message("An additional ", length(select_genes),
                " genes will be highlighted.")
    }

    # Select the genes to label. You can specify if you want to only label
    # annotated genes (default), i.e. those that have a gene name. This is for
    # "auto" and "select" labeling, i.e. it does not apply to manual labeling
    # (e.g. if you want to manually label an unannotated gene).
    if (remove_unannotated) {
        possible_labels <- res %>%
            filter(!grepl("ENSG", gene_name)) %>%
            pull(gene_name)
    } else {
        possible_labels <- res$gene_name
    }

    # Auto-labeling (default): label the top 5 (or n) up and down-regulated
    # genes. The "top" genes are those with the highest fold change and smallest
    # p-value.
    if (label == "auto") {
        upgenes <- up_df %>%
            arrange(desc(log2FoldChange ^ 2 * log10(padj) ^ 2)) %>%
            filter(gene_name %in% possible_labels) %>%
            head(n) %>%
            dplyr::select(gene_name) %>%
            unlist()

        downgenes <- down_df %>%
            arrange(desc(log2FoldChange ^ 2 * log10(padj) ^ 2)) %>%
            filter(gene_name %in% possible_labels) %>%
            head(n) %>%
            dplyr::select(gene_name) %>%
            unlist()

        # Record which genes to label
        res <- res %>% mutate(
            label = case_when(gene_name %in% c(upgenes, downgenes) ~ gene_name)
        )
    }

    # Selective labeling ("select"): label the top n up/down genes that are in
    # `select_genes`
    if (label == "select") {
        upgenes <- up_df %>%
            filter(in_list == "Y", gene_name %in% possible_labels) %>%
            arrange(desc(log2FoldChange ^ 2 * log10(padj) ^ 2)) %>%
            head(n) %>%
            dplyr::select(gene_name) %>%
            unlist()

        downgenes <- down_df %>%
            filter(in_list == "Y", gene_name %in% possible_labels) %>%
            arrange(desc(log2FoldChange ^ 2 * log10(padj) ^ 2)) %>%
            head(n) %>%
            dplyr::select(gene_name) %>%
            unlist()

        res <- res %>% mutate(
            label = case_when(gene_name %in% c(upgenes, downgenes) ~ gene_name)
        )
    }

    # Manual labeling ("manual"): label the genes you provided in `manual_genes`
    if (label == "manual") {
        res <- res %>% mutate(label = case_when(
            ensg_id %in% manual_genes | gene_name %in% manual_genes ~ gene_name
        ))
    }

    # Create the plot
    p <- ggplot(res, aes(x = log2FoldChange, y = neglogp)) +

        # Plot the non-significant genes
        geom_point(
            data = res %>% filter(significant == "NS", in_list == "N"),
            alpha = alpha,
            show.legend = FALSE,
            size = point_size,
            colour = nonsig_colour
        ) +

        # Plot the significant genes and genes of interest ("select_genes"),
        # with those in `select_genes` overlaying those not in `select_genes`
        # for emphasis.
        geom_point(
            data = res %>% filter(significant == "SIG", in_list == "N"),
            mapping = aes(x = log2FoldChange, y = neglogp),
            size = point_size,
            alpha = alpha,
            colour = base_colour
        ) +

        geom_point(
            data = res %>% filter(in_list == "Y"),
            mapping = aes(x = log2FoldChange, y = neglogp),
            size = point_size,
            alpha = alpha,
            colour = select_colour
        ) +

        # Add cutoff lines
        geom_hline(
            yintercept = -log10(p_cutoff),
            linetype = "dashed",
            colour = "gray20"
        ) +
        geom_vline(
            xintercept = c(log2(fc_cutoff), -log2(fc_cutoff)),
            linetype = "dashed",
            colour = "gray20"
        ) +

        # Some more graph adjustments: set to a clean theme, make labels easy to
        # read
        theme_bw() +
        theme(
            axis.text = element_text(colour = "black", size = 11),
            axis.title = element_text(
                size = 13, face = "bold", colour = "black"
            ),
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(size = 14),
            plot.background = element_blank()
        ) +

        labs(
            x = expression(bold(log["2"]~Fold~Change)),
            y = expression(bold(-log["10"]~P[adj]))
        ) +

        # Set axes
        {if(!is.na(yaxis[1])) ylim(yaxis)} +
        {if(!is.na(xaxis[1])) xlim(xaxis)} +

        # Add in labels
        geom_text_repel(
            data = res %>% filter(!is.na(label)),
            mapping = aes(label = label),
            segment.color = "darkgrey",
            color = "black",
            bg.color = "white",
            bg.r = .05,
            max.overlaps = Inf,
            fontface = "bold",
            size = label_size,
            box.padding = pad
        ) +

        # Add in informative subtitles for number of up and down-regulated
        # genes, as well as number of genes in "select_genes," add a title if
        # provided
        {if(!is.na(title)) labs(title = title)} +

        # If there are no "select_genes"
        {
            if (length(select_genes) == 0)
                labs(subtitle = paste0(
                    "Down: ", num_genes[2], ", Up: ", num_genes[1]
                ))
        } +

        # If "select_genes" was given
        {
            if (length(select_genes) != 0)
                labs(subtitle = paste0(
                    "Down: ", num_genes[2],
                    ", Up: ", num_genes[1], ", ",
                    select_name, ": ", length(select_genes)
                ))
        }

    # For plotting with fold change, not log2 fold change, manually add in
    # scientific notation based on the xaxis range.
    if (nonlog2) {

        # If "xaxis" is not specified
        if (is.na(xaxis[1])) {
            xaxis <- c(
                min(res$log2FoldChange) - 0.1,
                max(res$log2FoldChange) + 0.1
            )
        }

        p <- p +
            xlab("Fold Change") +
            theme(axis.text.x = element_text(vjust = 0, size = 11)) +
            {
                if (between(max(abs(xaxis)), 0, 6))
                    # Up to log2FC = 6
                    scale_x_continuous(
                        breaks = c(-6, -4, -2, 0, 2, 4, 6),
                        labels = c(-32, -16, -4, 1, 4, 16, 32),
                        limits = xaxis
                    )
            } +
            {
                if (between(max(abs(xaxis)), 6.001, log2(100)))
                    # Up to log2FC = 6.64
                    scale_x_continuous(
                        breaks = c(
                            -log2(100), -log2(10), 0, log2(10), log2(100)
                        ),
                        labels = c(
                            expression(-10 ^ 2),
                            expression(-10 ^ 1),
                            1,
                            expression(10 ^ 1),
                            expression(10 ^ 2)
                        ),
                        limits = xaxis
                    )
            } +
            {
                if (between(max(abs(xaxis)), log2(100), log2(1000)))
                    # Up to log2FC = 9.97
                    scale_x_continuous(
                        breaks = c(
                            -log2(1000),
                            -log2(100),
                            -log2(10),
                            0,
                            log2(10),
                            log2(100),
                            log2(1000)
                        ),
                        labels = c(
                            expression(-10 ^ 3),
                            expression(-10 ^ 2),
                            expression(-10 ^ 1),
                            1,
                            expression(10 ^ 1),
                            expression(10 ^ 2),
                            expression(10 ^ 3)
                        ),
                        limits = xaxis
                    )
            } +
            {
                if (between(max(abs(xaxis)), log2(1000), log2(10000)))
                    # Up to log2FC = 13.29
                    scale_x_continuous(
                        breaks = c(
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
                        labels = c(
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
                        limits = xaxis
                    )
            } +
            {
                if (between(max(abs(xaxis)), log2(10000), log2(100000)))
                    # Up to log2FC = 16.61
                    scale_x_continuous(
                        breaks = c(
                            -log2(10 ^ 5),
                            -log2(10 ^ 3),
                            -log2(10),
                            0,
                            log2(10),
                            log2(10 ^ 3),
                            log2(10 ^ 5)
                        ),
                        labels = c(
                            expression(-10 ^ 5),
                            expression(-10 ^ 3),
                            expression(-10 ^ 1),
                            1,
                            expression(10 ^ 1),
                            expression(10 ^ 3),
                            expression(10 ^ 5)
                        ),
                        limits = xaxis
                    )
            } +
            {
                if (max(abs(xaxis)) >= log2(100000))
                    message(
                        "Something may be wrong with your DESeq model to have ",
                        "fold changes >10^5..."
                    )
            }
    }

    return(p)
}
