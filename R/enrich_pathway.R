# pathway enrichment
# to-dos
# 1. how to make function continue to run if sigora fails?
# 2. include other pathway enrichment: hallmark, reactomePA

library(sigora)
library(ReactomePA)
library(msigdbr)
library(tidyverse)

# Function to clean-up SIGORA output ---------------------------------------
run_sigora <- function(enrich_genes, # vector of genes to enrich
                       direction = NA # if up or down-regulated genes were used
                       ){

  # run SIGORA based on default settings (GPSrepo = reaH, levels = 4, adjusted pvalue cutoff < 0.001)
  invisible(capture.output(sigora_data <- sigora(GPSrepo = reaH, level = 4, markers = TRUE, queryList = enrich_genes)))
  sigora_results <- sigora_data$summary_results %>% filter(Bonferroni < 0.001)

  # not all DE genes are used in the Reactome pathway genes from SIGORA
  # SIGORA does not have functionality to calculate gene ratio, which might be useful, but it is basically a reflection of p-value and pathway size
  # The gene ratio is somewhat confusing to calculate. Gene ratio = k/n where:
  # k = analyzed genes (the genes in enrich_genes) that are found in the pathway of interest
  # n = analyzed genes that are found in the entire database (all genes used in SIGORA in the Reactome database)
  # calculate n
  n_genes <- enrich_genes[enrich_genes %in% idmap$Ensembl.Gene.ID] %>% length()

  if(nrow(sigora_results) > 0){ # sometimes no pathways are significantly enriched
    sigora_results$direction <- direction

    # get the DE genes that were enriched for each pathway
    for(i in 1:length(sigora_results$pathwy.id)){

      pathid <- sigora_results$pathwy.id[i]

      # first, get the enriched gene pairs for that pathway
      geneid <- sigora_data$detailed_results %>% as.data.frame() %>% filter(pathway == pathid)
      unique_entrez <- unique(c(geneid$gene1, geneid$gene2)) # the unique entrez ids from the gene pairs

      # convert them into ensg_id using SIGORA's internal mapping file, keeping the genes that were inputted
      # there may be some discrepancies between the number of unique entrezIDs and ensgIDs due to duplicate mapping
      # use the number of ensgIDs as that was the input
      unique_ensg <- unique(idmap %>%
                              filter(EntrezGene.ID %in% unique_entrez & Ensembl.Gene.ID %in% enrich_genes) %>%
                              .$Ensembl.Gene.ID %>%
                              as.character())
      # or, use our own mapping file for conversion
      # unique_ensg <- unique(mapping_file %>% filter(entrez_id %in% unique_entrez & ensg_id %in% enrich_genes) %>% .$ensg_id)

      pathway_info_df <- data.frame(pathwy.id = pathid,
                                    candidate_genes = paste(unique_ensg, collapse = ';'),
                                    num_candidate_genes = length(unique_ensg)) %>%
        mutate(gene_ratio = num_candidate_genes/n_genes)

      if(i == 1){
        final_results <- pathway_info_df
      } else {
        final_results <- rbind(final_results, pathway_info_df)
      }
    }
  }

  final_output <- left_join(sigora_results, final_results, by = 'pathwy.id') %>%
    transmute(
      pathway_id = pathwy.id,
      pathway_description = description,
      direction,
      p_value = pvalues,
      p_value_adjusted = Bonferroni,
      candidate_genes,
      num_candidate_genes,
      num_bg_genes = n_genes,
      gene_ratio
    )

}

# Function to enrich pathway ----------------------------------------------
enrich_pathway <- function(deseq_result_list, # list of dataframes of DESeq2 results
                                              # list name as the comparison for each dataframe (e.g. COVID vs Healthy)
                                              # dataframes have ensembl gene ids as rownames
                           p_cutoff = 0.05, # adjusted pvalue cutoff, default <0.05
                           fc_cutoff = 1.5, # absolute fold change, default abs >1.5
                           split = TRUE, # split into up and down-regulated DEGs and do enrichment separately
                           analysis = 'sigora' # default is SIGORA ('sigora'), others: ReactomePA ('reactomepa'), mSigDB Hallmark gene sets ('hallmark')
                           ){

  if(!is.list(deseq_result_list)){
    print('Provide a list of dataframes of DESeq2 results with the name of each item in the list as the comparison name.')
    return()
  } else {
    for(i in 1:length(deseq_result_list)){
      comparison <- names(deseq_result_list[i])
      # use DE genes
      deseq_results <- deseq_result_list[i][[1]] %>% filter(padj < p_cutoff, abs(log2FoldChange) > log2(fc_cutoff))

      # SIGORA analysis
      # SIGORA uses gene pairs with the Reactome database, which often eliminates duplicate, closely related pathways (e.g. TLR7/8/9 pathways).
      # Useful for if you have a lot of DEGs that might enrich for a lot of pathways, making it difficult analyze

      if(analysis == 'sigora'){
        print('Enrichment using SIGORA')
        print(paste0('Comparison being analyzed: ', comparison))
        # if up and downregulated genes are to be enriched separately

        if(split){
          up_gns <- filter(deseq_results, log2FoldChange > 0) %>% rownames()
          dn_gns <- filter(deseq_results, log2FoldChange < 0) %>% rownames()
          print(paste0('DEGs used: Up: ', length(up_gns), ', Down: ', length(dn_gns))) # Print how many DE genes used

          # run SIGORA
          up_results <- run_sigora(up_gns, direction = 'Up')
          dn_results <- run_sigora(dn_gns, direction = 'Down')
          total_results <- rbind(up_results, dn_results)

        } else { # just use all DE genes for enrichment
          de_gns <- rownames(deseq_results)
          print(paste0('DEGs used: ', length(de_gns))) # Print how many DE genes used
          total_results <- run_sigora(de_gns, direction = 'All')
        }

        total_results_annotated <- left_join(total_results, top_pathways %>% select(pathway_id, top_pathways), by = 'pathway_id')
        total_results_annotated$comparison <- comparison
        total_results_annotated$total_genes <- nrow(deseq_results)
        print(paste0('Done! Enriched pathways: ', nrow(total_results_annotated)))
      }

    # This adds the next dataframe to the previous one
      if(i == 1){
        final_enriched_results <- total_results_annotated
      } else {
        final_enriched_results <- rbind(final_enriched_results, total_results_annotated)
      }
    }
  }

  final_enriched_results$comparison <- factor(final_enriched_results$comparison,
                                              levels = c(names(deseq_result_list)))

  # Return the results, which is a tibble of all the pathway results
  return(tibble(final_enriched_results))

}

# test_list <- list('Time 1' = deseq_example_1, 'Time 2' = deseq_example_2)
enriched_results <- enrich_pathway(deseq_example_list[c(5,6)])

# if(analysis == 'reactomePA'){
#   mapping <- final_map_condense %>% column_to_rownames(var = 'ensg_id')
#   dataframe <- merge(dataframe, mapping, by = 0, all.x = TRUE)
#
#   up <- dplyr::filter(dataframe, FoldChange > 0)
#   down <- dplyr::filter(dataframe, FoldChange < 0)
#   print(paste0('up: ', length(row.names(up)), ', down: ', length(row.names(down)))) # Print how many DE genes
#
#   up <- enrichPathway(up$entrez_id, readable = TRUE) %>% as.data.frame()
#   if(length(row.names(up) > 0)){
#     up$id=filename
#     up$direction='Up'
#   }
#   down <- enrichPathway(down$entrez_id, readable = TRUE) %>% as.data.frame()
#   if(length(row.names(down) > 0)){
#     down$id=filename
#     down$direction='Down'
#   }
#   total <- rbind(up, down)
#   total$`-log10 p-adj`=(-log((total$p.adjust+1E-300), 10))
#   names(total)[1] <- 'pathwy.id'
#   names(total)[2] <- 'description'
#
#   # Annotate the "top" pathways
#   total[, c("Top Pathway Name", 'Hierarchy')] <- NA # Generates three new columns to fill
#   for (i in c(1:length(row.names(total)))) {
#     total[i , c("Top Pathway Name", 'Hierarchy')] <- path_steps_local(total[i,1]) #Fills in column per row based on the first column (pathway id)
#     print(i)
#   }
#
#   # Save as a table
#   write.table(total, paste0('react_pathways_',filename,'.csv'), sep = ',', row.names=FALSE)
#
# }
#
# if(analysis == 'hallmark'){
#   database <- hallmark
#   msigdbr_t2g = database %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
#
#   up <- dplyr::filter(dataframe, FoldChange > 0)
#   down <- dplyr::filter(dataframe, FoldChange < 0)
#   print(paste0('up: ', length(row.names(up)), ', down: ', length(row.names(down)))) # Print how many DE genes
#
#   up <- enricher(gene = rownames(dataframe %>% filter(log2FoldChange > 0)), TERM2GENE = msigdbr_t2g, universe = universe) %>% as.data.frame()
#   if(length(rownames(up)) != 0){
#     up$direction <- 'Up'
#   }
#
#
#   down <- enricher(gene = rownames(dataframe %>% filter(log2FoldChange < 0)), TERM2GENE = msigdbr_t2g, universe = universe) %>% as.data.frame()
#   if(length(rownames(down)) != 0){
#     down$direction <- 'Down'
#   }
#
#   total <- rbind(up, down)
#   total$id <- filename
#   total$`-log10 p-adj`= -log10(total$p.adjust +1E-300)
#
# }
