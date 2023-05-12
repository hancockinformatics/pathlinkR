#' INTERNAL run_sigora
#'
#' @param enrich_genes Vector of genes to enrich
#' @param direction If up or down-regulated genes were used
#' @param gps_repo GPS object to use for testing pathways
#'
#' @return Data frame of results from Sigora
#'
#' @import dplyr
#' @import sigora
#'
run_sigora <- function(
    enrich_genes,
    direction = NA,
    gps_repo
) {

  # Run SIGORA based on default settings (GPSrepo = reaH, levels = 4, adjusted pvalue cutoff < 0.001)
  invisible(capture.output(
    sigora_data <- sigora(
      GPSrepo = gps_repo,
      level = 4,
      markers = TRUE,
      queryList = enrich_genes
    )
  ))

  sigora_results <- sigora_data$summary_results %>% filter(Bonferroni < 0.001)

  # Not all DE genes are used in the Reactome pathway genes from SIGORA. SIGORA
  # does not have the functionality to calculate gene ratio, which might be
  # useful, but it is basically a reflection of p-value and pathway size. The
  # gene ratio is somewhat confusing to calculate. Here, gene ratio = k/n where:
  #   k = analyzed genes (the genes in enrich_genes) that are found in the
  #       pathway of interest
  #   n = analyzed genes that are found in the entire database (all genes used
  #       in SIGORA in the Reactome database)

  # Calculate n
  n_genes <- enrich_genes[enrich_genes %in% idmap$Ensembl.Gene.ID] %>% length()

  if (nrow(sigora_results) > 0) {

    sigora_results$direction <- direction

    # Get the DE genes that were enriched for each pathway
    for (i in 1:length(sigora_results$pathwy.id)) {

      pathid <- sigora_results$pathwy.id[i]

      # First, get the enriched gene pairs for that pathway
      geneid <- sigora_data$detailed_results %>%
        as.data.frame() %>%
        filter(pathway == pathid)

      # The unique Entrez IDs from the gene pairs
      unique_entrez <- unique(c(geneid$gene1, geneid$gene2))

      # Convert them into HGNC symbols using SIGORA's internal mapping file,
      # keeping the genes that were input. There may be some discrepancies
      # between the number of unique Entrez IDs HGNC symbols due to duplicate
      # mapping. Also, use the ones that also were in the original ENSG input.

      unique_hgnc <- unique(
        idmap %>%
          filter(
            EntrezGene.ID %in% unique_entrez & Ensembl.Gene.ID %in% enrich_genes
          ) %>%
          .$Symbol %>%
          as.character()
      )

      # remove any that did not map to a gene symbol
      unique_hgnc <- unique_hgnc[unique_hgnc != '']

      # Old code that converted the entrez IDs to ENSG IDs, may implement
      # unique_ensg <- unique(
      #   idmap %>%
      #     filter(
      #       EntrezGene.ID %in% unique_entrez & Ensembl.Gene.ID %in% enrich_genes
      #     ) %>%
      #     .$Ensembl.Gene.ID %>%
      #     as.character()
      # )

      # Or, use our own mapping file for conversion
      # unique_ensg <- unique(
      #   mapping_file %>% filter(
      #     entrez_id %in% unique_entrez &
      #       ensg_id %in% enrich_genes
      #   ) %>% .$ensg_id
      # )
      #

      pathway_info_df <- data.frame(
        pathwy.id = pathid,
        genes = paste(unique_hgnc, collapse = ';'),
        num_candidate_genes = length(unique_hgnc)
      ) %>%
        mutate(gene_ratio = num_candidate_genes/n_genes)

      if (i == 1) {
        final_results <- pathway_info_df
      } else {
        final_results <- rbind(final_results, pathway_info_df)
      }
    }
  }

  final_output <-
    left_join(sigora_results, final_results, by = 'pathwy.id') %>%
    transmute(
      pathway_id = pathwy.id,
      pathway_description = description,
      direction,
      p_value = pvalues,
      p_value_adjusted = Bonferroni,
      genes,
      num_candidate_genes,
      num_bg_genes = n_genes,
      gene_ratio
    )
}
