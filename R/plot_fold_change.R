#' Plot fold change plots
#'
#' @param input_list List of data frames of DESeq2 results. The list names
#'   are used as the comparison for each dataframe (e.g. COVID vs Healthy). Data
#'   frames should have Ensembl gene IDs as rownames.
#' @param path_name Name of pathway to pull genes from, will be plot title
#' @param path_id ID of pathway to pull genes from, if path_name is not given
#' @param manual_title Set this as your title instead of pathway name
#' @param title_size Font size for title
#' @param genes_to_plot Set if you want a specific vector of genes instead of
#'   pulling from a pathway
#' @param gene_format Default is Ensembl gene IDs ("ensg"), can also input a
#'   vector of HGNC symbols ("hgnc")
#' @param p_cutoff P value cutoff. Default is p < 0.05
#' @param fc_cutoff Fold change cutoff. Default is |FC| > 1.5
#' @param plot_significant_only Boolean (TRUE) Only plot genes that are
#'   differentially expressed (pass p_cutoff and fc_cutoff) in any comparison
#' @param hide_low_fc Boolean (TRUE) If a gene is significant in one comparison
#'   but not in another, this will set the colour of the non-significant gene
#'   as grey to visually emphasize the significant genes. If set to FALSE, it
#'   will set the colour to the fold change, and if the p value passes p_cutoff,
#'   it will also display the p value (the asterisks will be grey instead of
#'   black).
#' @param vjust Adjustment of the position of the significance stars. Default is
#'   0.75. May need to adjust if there are many genes
#' @param rot Rotation of the position of the significance stars. Default is 0
#' @param invert Boolean (FALSE) Default plots genes as rows and conditions as
#'   columns, set to TRUE if you want genes as columns and conditions as rows
#' @param log2_foldchange Boolean (FALSE) Default plots the fold changes in the
#'   legend as the true fold change, set to TRUE if you want log2 fold change
#' @param col_split Split each condition into groups for better visualization.
#'   To do so, create a vector where each item of the vector corresponds to
#'   which group the dataframe belongs to in input_list. E.g. ('Positive',
#'   'Positive', "Negative", "Negative", "Time", "Time"). The order of these
#'   groups is set in the order they appear; if you want to change the order,
#'   change the order of data frames in input_list and write col_split in the
#'   order you want. Order will be ignored if cluster_columns is set to TRUE
#' @param cluster_rows Boolean (TRUE) Whether to cluster the rows (genes). May
#'   need to change if invert = TRUE.
#' @param cluster_columns Boolean (FALSE) Whether to cluster the columns
#'   (conditions). Will override order of col_split if set to TRUE.
#' @param col_angle angle of column text. Set default to 90
#' @param col_center whether to center column text. Default is TRUE, should set
#'   to FALSE if angled column name (e.g. col_angle = 45)
#' @param row_angle angle of row text. Set default to 0
#' @param row_center whether to center column text. Default is FALSE, should set
#'   to TRUE if vertical column name (e.g. row_angle = 90)
#'
#'
#' @return A heatmap of fold changes for genes of interest
#' @export
#'
#' @import ComplexHeatmap
#' @import dplyr
#'
plot_fold_change <- function(
    input_list,
    path_name = NA,
    path_id = NA,
    manual_title = NA,
    title_size = 14,
    genes_to_plot = NA,
    gene_format = 'ensg',
    p_cutoff = 0.05,
    fc_cutoff = 1.5,
    plot_significant_only = TRUE,
    hide_low_fc = TRUE,
    vjust = 0.75,
    rot = 0,
    invert = FALSE,
    log2_foldchange = FALSE,
    col_split = NA,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    col_angle = 90,
    col_center = TRUE,
    row_angle = 0,
    row_center = FALSE) {

  # First identify the pathway to plot
  ## If pathway name is provided
  if (!is.na(path_name)) {
    path_id <- sigora_database %>%
      filter(pathway_name == path_name) %>%
      select(pathway_id) %>%
      unlist() %>% .[1] %>% as.character()
    plot_title <- path_name
  }

  ## If pathway ID is provided
  if (!is.na(path_id)) {
    plot_title <- sigora_database %>%
      filter(pathway_id == path_id) %>%
      select(pathway_name) %>%
      unlist() %>% .[1] %>% as.character()
  }

  ## If a title is provided manually, overwrite
  if(!is.na(manual_title)){
    plot_title <- manual_title
  }

  # Get the genes to plot in the pathway or genelist of interest
  # get all the genes in the pathway of interest, and make them ensg ids
  if (is.na(genes_to_plot[1])){ # get from pathway database
    genes <- sigora_database %>%
      filter(pathway_id == path_id) %>%
      .$Ensembl.Gene.ID
  } else { # get from manual gene input
    genes <- genes_to_plot
    if (gene_format == 'hgnc') {
      genes <- mapping_file %>%
        filter(gene_name %in% genes_to_plot) %>%
        .$ensg_id
    }
    if (gene_format == 'ensg') {
      genes <- genes_to_plot
    }
  }

  # Get fold changes and significance for each dataframe in input_list

  ## Collection of significant genes in each dataframe
  sig_genes <- c()

  ## Loop across each dataframe in the input_list
  for(n in 1:length(input_list)){
    ## Get fold changes for each gene of interest
    fold_change <- input_list[[n]] %>%
      filter(rownames(.) %in% genes, !is.na(log2FoldChange)) %>%
      select(log2FoldChange) %>%
      rownames_to_column()
    names(fold_change) <- c('ensg_id', names(input_list)[n])

    ## Get significance values for each gene of interest
    signif <- input_list[[n]] %>%
      filter(rownames(.) %in% genes, !is.na(padj)) %>%
      select(padj) %>%
      rownames_to_column()
    names(signif) <- c('ensg_id', names(input_list)[n])

    ## Add to the fold change and p value dataframes
    if (n == 1) {
      df_fc <- fold_change
      df_p <- signif
    } else {
      df_fc <- full_join(df_fc, fold_change, by = 'ensg_id')
      df_p <- full_join(df_p, signif, by = 'ensg_id')
    }

    ## If any are significantly DE, record them
    sig_genes <- c(
      sig_genes,
      input_list[[n]] %>%
        rownames_to_column(var = 'ensg_id') %>%
        filter(
          ensg_id %in% genes,
          padj < p_cutoff,
          abs(log2FoldChange) > log2(fc_cutoff)
        ) %>%
      .$ensg_id)
  }

  # From all the dataframes, get the genes that were significant in any
  sig_genes <- unique(sig_genes)

  if (plot_significant_only) {
    df_fc <- df_fc %>% filter(ensg_id %in% sig_genes)
    df_p <- df_p %>% filter(ensg_id %in% sig_genes)
  }

  # Prepare the Heatmap matrices
  ## Map the ensg_id to hgnc symbols
  mat_fc <- df_fc %>%
    left_join(mapping_file) %>%
    select(!ensg_id) %>%
    select(!entrez_id) %>%
    column_to_rownames(var = 'gene_name') %>%
    as.matrix()
  mat_fc[is.na(mat_fc)] <- 0 # make any NAs into 0
  if (hide_low_fc) {
    mat_fc[abs(mat_fc) < log2(fc_cutoff)] <- 0
  }

  mat_p <- df_p %>%
    left_join(mapping_file) %>%
    select(!ensg_id) %>%
    select(!entrez_id) %>%
    column_to_rownames(var = 'gene_name') %>%
    as.matrix()
  mat_p[is.na(mat_p)] <- 1 # make any NAs into 1s

  # Set the limits for colours for the plotting heatmap
  limit = max(abs(mat_fc), na.rm = TRUE) %>% ceiling()
  ## If plotting real fold changes instead of log2
  if(!log2_foldchange){
    if ((limit %% 2) == 0) {
      range = c(-limit, -limit/2, 0, limit/2, limit)
    } else {
      limit <- limit + 1
      range = c(-limit, -limit/2, 0, limit/2, limit)
    }
    foldchange_title = 'Fold\nChange'
    labels = c(-2^limit, -2^(limit/2), 1, 2^(limit/2), 2^limit)
  } else {
    range = c(-limit, -limit/2, 0, limit/2, limit)
    foldchange_title = 'log2 FC'
    labels = range
  }
  parameters <- list(
    at = range,
    labels = labels,
    title = foldchange_title)

  # If columns aren't being split
  if (is.na(col_split[1])) {
    col_split <- rep(NA, length(input_list))
    column_title <- NULL
  } else {
    # this orders the splitting into the order the dataframes are in the list
    # instead of alphabetically
    col_split <- factor(col_split, levels = unique(col_split))
    column_title <- '%s'
  }
  row_split <- rep(NA, nrow(mat_fc))
  row_title <- NULL

  # If plotting so that rownames are conditions and colnames are genes
  if (invert) {
    mat_fc <- t(mat_fc)
    mat_p <- t(mat_p)
    var <- col_split
    col_split <- row_split
    row_split <- var
    column_title <- NULL
    row_title <- '%s'
  }

  draw(Heatmap(
    mat_fc,
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (abs(mat_fc[i,j]) > log2(1.5)) {
        if (mat_p[i, j] < 0.001) {
          grid::grid.text('***', x, y, vjust = vjust, rot = rot)
        }
        else if (mat_p[i, j] < 0.01) {
          grid::grid.text('**', x, y, vjust = vjust, rot = rot)
        }
        else if (mat_p[i, j] < 0.05) {
          grid::grid.text('*', x, y, vjust = vjust, rot = rot)
        }
          # as.character(expression('\u2736') # this doesn't work?
      }

      # If plotting significance values for genes that don't pass fc_cutoff
      if (abs(mat_fc[i,j]) < log2(1.5) & !hide_low_fc) {
        if (mat_p[i, j] < 0.001) {
          grid::grid.text("***", x, y, vjust = vjust, rot = rot, gp = grid::gpar(col = 'grey50'))
        }
        else if (mat_p[i, j] < 0.01) {
          grid::grid.text("**", x, y, vjust = vjust, rot = rot, gp = grid::gpar(col = 'grey50'))
        }
        else if (mat_p[i, j] < 0.05) {
          grid::grid.text("*", x, y, vjust = vjust, rot = rot, gp = grid::gpar(col = 'grey50'))
        }
      }
      },
    column_title = column_title,
    row_title = row_title,
    heatmap_legend_param = parameters,
    column_title_gp = grid::gpar(fontsize = title_size),
    row_split = row_split,
    column_split = col_split,
    #col = circlize::colorRamp2(c(-limit, 0, limit), c("blue", "gray90", "red")),
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    column_names_rot = col_angle,
    column_names_centered = col_center,
    row_names_rot = row_angle,
    row_names_centered = row_center

  ), column_title = plot_title)

}
