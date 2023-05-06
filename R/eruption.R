# Volcano plot

# load example DESeq2 outputs

library(tidyverse)
library(ggrepel)

eruption <- function(deseq_results, # dataframe of DESeq2 results with ensembl gene ids
                     p_cutoff = 0.05, # adjusted pvalue cutoff, default <0.05
                     fc_cutoff = 1.5, # absolute fold change, default abs >1.5
                     base_colour = 'steelblue4', # colour of points for DE genes
                     nonsig_colour = 'lightgrey', # colour of non-DE genes
                     alpha = 0.5, # transparency of points
                     point_size = 1, # size of points

                     title = NA, # title of plot, be informative (e.g. COVID vs Healthy)
                     absolute = FALSE, # make x-axis absolute fold change (TRUE) instead of log2fold change (FALSE)
                     xaxis = NA, # manually specify limits of x-axis (in log2FoldChange), e.g. xaxis = c(-3, 5)
                     yaxis = NA, # manually specify limits of y-axis (in -log10), e.g. yaxis = c(0, 5)

                     select_genes = c(), # vector of select genes to emphasize by colouring differently (e.g. genes of interest), should be ensg ids
                     select_colour = 'red', # colour for the select genes
                     select_name = 'Selected', # what to call select genes (e.g. Unique, Shared, Immune Related, etc.)

                     label = 'auto', # auto: label the top 5 (or n) up/down annotated DE genes,
                                     # select: label top 5 (or n) up/down genes provided in select_genes
                                     # manual: label a manual selection of genes provided in manual_genes
                     manual_genes = c(), # if label = 'manual', these are the genes to specifically label, should be ensg ids
                     remove_unannotated = TRUE, # remove genes without annotation (i.e. no hgnc symbol), this is default
                     n = 5, # number of top up and downregulated genes to label, applies to label = 'auto' or 'select'
                     label_size = 3.5, # size of font for labels
                     pad = 1.4 # padding of labels (adjust this if labels overlap)

){

  # Annotate ensembl gene ids with gene names from the mapping file
  # for ensembl gene ids without gene names, just use the ensembl gene id
  res <- left_join(deseq_results %>% rownames_to_column(var = 'ensg_id'), mapping_file) %>%
    mutate(gene_name = case_when(is.na(gene_name) ~ ensg_id,
                                 !is.na(gene_name) ~ gene_name),
           # determine which are significant
           significant = case_when(padj < p_cutoff & abs(log2FoldChange) > log2(fc_cutoff) ~ 'SIG',
                                   TRUE ~ 'NS'),
           in_list = case_when(ensg_id %in% select_genes ~ 'Y',
                               TRUE ~ 'N'),
           neglogp = -log10(padj)) %>%
    filter(!is.na(padj)) # remove any genes without a pvalue

  # If specifying x and y axis, remove any genes that fall outside the specified ranges
  if(!is.na(xaxis[1])){
    res <- res %>% filter(log2FoldChange > xaxis[1] & log2FoldChange < xaxis[2])
  }
  if(!is.na(yaxis[1])){
    res <- res %>% filter(neglogp > yaxis[1] & neglogp < yaxis[2])
  }

  # make sure select_genes, if used, are in the dataframe
  select_genes <- select_genes[select_genes %in% res$ensg_id]

  # Identify top up and down regulated genes, and output them for the user
  up_df <- res %>% filter(padj < p_cutoff, log2FoldChange > log2(fc_cutoff))
  down_df <- res %>% filter(padj < p_cutoff, log2FoldChange < log2(1/fc_cutoff))
  num_genes <- c(nrow(up_df), nrow(down_df), nrow(up_df) + nrow(down_df)) # number of up, down, and total DE genes
  print(paste0('Creating volcano plot with ',
               num_genes[3], ' DEGs: ',
               num_genes[1], ' up-regulated DEGs and ',
               num_genes[2], ' down-regulated DEGs.'))
  if(length(select_genes) > 0){
    print(paste0('An additional ', length(select_genes), ' genes are highlighted.'))}

  # Select the genes to label
  ## specify if you want to only label annotated genes (default), i.e., those that have a gene name
  ## this is for auto and select labelling, does not apply to manual labelling (e.g. if you want to manually label an unannotated gene)
  if(remove_unannotated){
    possible_labels <- res %>% filter(!grepl('ENSG', gene_name)) %>% .$gene_name
  } else {
    possible_labels <- res$gene_name
  }

  ## Auto-labelling (default): label the top 5 (or n) up and downregulated genes
  ## The "top" genes are those with the highest fold change and smallest p-values
  if(label == 'auto'){
    upgenes <- up_df %>%
      arrange(desc(log2FoldChange^2 * log10(padj)^2)) %>%
      filter(gene_name %in% possible_labels) %>% # remove unannotated genes if necessary
      head(n) %>% dplyr::select(gene_name) %>% unlist() # Label the top n genes
    downgenes <- down_df %>%
      arrange(desc(log2FoldChange^2 * log10(padj)^2)) %>%
      filter(gene_name %in% possible_labels) %>% # remove unannotated genes if necessary
      head(n) %>% dplyr::select(gene_name) %>% unlist() # Label the top n genes

    # record which genes to label
    res <- res %>% mutate(label = case_when(gene_name %in% c(upgenes, downgenes) ~ gene_name))
  }

  ## Selective labelling (select): label the top n up/down genes that are in select_genes
  if(label == 'select'){
    upgenes <- up_df %>% filter(in_list == 'Y', gene_name %in% possible_labels) %>%
      arrange(desc(log2FoldChange^2 * log10(padj)^2)) %>%
      head(n) %>% dplyr::select(gene_name) %>% unlist()
    downgenes <- down_df %>% filter(in_list == 'Y', gene_name %in% possible_labels) %>%
      arrange(desc(log2FoldChange^2 * log10(padj)^2)) %>%
      head(n) %>% dplyr::select(gene_name) %>% unlist()
    res <- res %>% mutate(label = case_when(gene_name %in% c(upgenes, downgenes) ~ gene_name))
  }

  ## Manual labelling (manual): label the genes you provided in manual_genes
  if(label == 'manual'){
    res <- res %>% mutate(label = case_when(gene_name %in% manual_genes ~ gene_name))
  }

  # Create the plot
  p <- ggplot(res, aes(x = log2FoldChange, y = neglogp)) +

    # plot non-significant genes
    geom_point(data = res %>% filter(significant == 'NS'), alpha = alpha, show.legend = FALSE, size = point_size, colour = nonsig_colour) +

    # plot significant genes, with those in select_genes overlaying those not in select_genes for emphasis
    geom_point(data = res %>% filter(significant == 'SIG', in_list == 'N'), aes(x = log2FoldChange, y = neglogp), size = point_size, alpha = alpha, colour = base_colour) +
    geom_point(data = res %>% filter(significant == 'SIG', in_list == 'Y'), aes(x = log2FoldChange, y = neglogp), size = point_size, alpha = alpha, colour = select_colour) +

    # add cutoff lines
    geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", colour = 'gray20') +
    geom_vline(xintercept = c(log2(fc_cutoff), -log2(fc_cutoff)), linetype="dashed", colour = 'gray20') +

    # some more graph adjustments
    ## set to a clean theme, make labels easy to read
    theme_bw() +
    theme(axis.text = element_text(
      colour="black",
      size=11),
      axis.title = element_text(size = 13, face = 'bold', colour = "black"),
      plot.title = element_text(face = 'bold', size = 16),
      plot.subtitle = element_text(size = 14),
      plot.background = element_blank()) +

    ylab(expression(bold(-log['10']~P[adj]))) + # rename y label
    xlab(expression(bold(log['2']~Fold~Change))) + # rename x label

    ## set axes
    {if(!is.na(yaxis[1])) ylim(yaxis)} +
    {if(!is.na(xaxis[1])) xlim(xaxis)} +

    ## add in labels
    geom_text_repel(data = res %>% filter(!is.na(label)),
                    aes(label = label),
                    segment.color = 'darkgrey',
                    color = "black",
                    bg.color = "white",
                    bg.r = .05,
                    max.overlaps = Inf,
                    fontface = 'bold',
                    size = label_size,
                    box.padding = pad) +

    ## add in informative subtitles for number of up and down-regulated genes, as well as number of genes in select_genes
    {if(!is.na(title)) labs(title = title)} + # add a title if provided
    {if(length(select_genes) == 0) labs(subtitle = paste0("Up: ", num_genes[1], ", Down: ", num_genes[2]))} + # if there are no select genes
    {if(length(select_genes) != 0) labs(subtitle = paste0("Up: ", num_genes[1], ", Down: ", num_genes[2], ", ", select_name, ": ", length(select_genes)))} # if select genes given

  # For plotting with fold change, not log2 fold change
  # manually add in scientific notation based on the xaxis range
  if(absolute){

    # if xaxis is not specified
    if(is.na(xaxis[1])){
      xaxis <- c(min(res$log2FoldChange)-0.1, max(res$log2FoldChange)+0.1)}

    p <- p +
      xlab('Fold Change') +
      theme(axis.text.x = element_text(vjust = 0, size = 11)) +

      {if(between(max(abs(xaxis)), 0, 6)) # up to log2FC = 6
        scale_x_continuous(breaks=c(-6, -4, -2, 0, 2, 4, 6),
                           labels=c(-32, -16, -4, 1, 4, 16, 32),
                           limits = xaxis)} +

      {if(between(max(abs(xaxis)), 6.001, log2(100))) #up to log2FC = 6.64
        scale_x_continuous(breaks=c(-log2(100), -log2(10), 0, log2(10), log2(100)),
                           labels=c(expression(-10^2), expression(-10^1), 1, expression(10^1), expression(10^2)),
                           limits = xaxis)} +

      {if(between(max(abs(xaxis)), log2(100), log2(1000))) #up to log2FC = 9.97
        scale_x_continuous(breaks=c(-log2(1000), -log2(100), -log2(10), 0, log2(10), log2(100), log2(1000)),
                           labels=c(expression(-10^3), expression(-10^2), expression(-10^1), 1, expression(10^1), expression(10^2), expression(10^3)),
                           limits = xaxis)} +

      {if(between(max(abs(xaxis)), log2(1000), log2(10000))) #up to log2FC = 13.29
        scale_x_continuous(breaks=c(-log2(10000), -log2(1000), -log2(100), -log2(10), 0, log2(10), log2(100), log2(1000), log2(10000)),
                           labels=c(expression(-10^4), expression(-10^3), expression(-10^2), expression(-10^1), 1, expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
                           limits = xaxis)} +

      {if(between(max(abs(xaxis)), log2(10000), log2(100000))) #up to log2FC = 16.61
        scale_x_continuous(breaks=c(-log2(10^5), -log2(10^3), -log2(10), 0, log2(10), log2(10^3), log2(10^5)),
                           labels=c(expression(-10^5), expression(-10^3), expression(-10^1), 1, expression(10^1), expression(10^3), expression(10^5)),
                           limits = xaxis)} +

      {if(max(abs(xaxis)) >= log2(100000))
        print('Something is probably wrong with your DESeq model to have huge fold changes >10^5.')}
  }

  p # return plot

}
