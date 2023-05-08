#' enrich_pathway
#'
#' @param deseq_result_list list of data frames of DESeq2 results. The list names
#'   are used as the comparison for each dataframe (e.g. COVID vs Healthy). Data
#'   frames should have Ensembl gene IDs as rownames.
#' @param p_cutoff Adjusted p value cutoff, defaults to < 0.05.
#' @param fc_cutoff Absolute fold change, defaults to |FC| > 1.5.
#' @param split Boolean (TRUE); Split into up and down-regulated DEGs and do
#'   enrichment separately
#' @param analysis Default is SIGORA ('sigora'), others: ReactomePA
#'   ('reactomepa'), mSigDB Hallmark gene sets ('hallmark')
#'
#' @return A data frame of pathway enrichment results
#' @export
#'
#' @import dplyr
#'
enrich_pathway <- function(
    deseq_result_list,
    p_cutoff = 0.05,
    fc_cutoff = 1.5,
    split = TRUE,
    analysis = "sigora"
) {

  # TODO
  # 1. How to make function continue to run if sigora fails?
  # 2. Include other pathway enrichment methods: hallmark, reactomePA

  if (!is.list(deseq_result_list)) {
    message("Provide a list of dataframes of DESeq2 results with the name of ",
            "each item in the list as the comparison name.")
    return()

  } else {
    for (i in 1:length(deseq_result_list)) {
      comparison <- names(deseq_result_list[i])

      # Use the DE genes
      deseq_results <- deseq_result_list[i][[1]] %>%
        filter(padj < p_cutoff, abs(log2FoldChange) > log2(fc_cutoff))


      ### SIGORA analysis
      # SIGORA uses gene pairs with the Reactome database, which often
      # eliminates duplicate, closely related pathways (e.g. TLR7/8/9 pathways).
      # Useful for if you have a lot of DEGs that might enrich for a lot of
      # pathways, making it difficult analyze.

      if (analysis == "sigora") {
        message("Enrichment using SIGORA")
        message(paste0("Comparison being analyzed: ", comparison))

        # If up- and down-regulated genes are to be enriched separately
        if (split) {
          up_gns <- filter(deseq_results, log2FoldChange > 0) %>% rownames()
          dn_gns <- filter(deseq_results, log2FoldChange < 0) %>% rownames()
          message(paste0(
            "DEGs used: Up: ",
            length(up_gns),
            ", Down: ", length(dn_gns)
          ))

          # run SIGORA using the internal function
          up_results <- run_sigora(up_gns, direction = "Up")
          dn_results <- run_sigora(dn_gns, direction = "Down")

          total_results <- rbind(up_results, dn_results)

        # Or just use all DE genes for enrichment
        } else {
          de_gns <- rownames(deseq_results)
          message(paste0("DEGs used: ", length(de_gns))) # Print how many DE genes used
          total_results <- run_sigora(de_gns, direction = "All")
        }

        total_results_annotated <- left_join(
          total_results,
          top_pathways %>% select(pathway_id, top_pathways),
          by = "pathway_id"
        )

        total_results_annotated$comparison <- comparison
        total_results_annotated$total_genes <- nrow(deseq_results)
        message(paste0(
          "Done! Enriched pathways: ", nrow(total_results_annotated)
        ))
      }

      # This adds the next dataframe to the previous one
      if(i == 1){
        final_enriched_results <- total_results_annotated
      } else {
        final_enriched_results <-
          rbind(final_enriched_results, total_results_annotated)
      }
    }
  }

  final_enriched_results$comparison <- factor(
    final_enriched_results$comparison,
    levels = c(names(deseq_result_list))
  )

  ### ReactomePA
  # if(analysis == "reactomePA"){
  #   mapping <- final_map_condense %>% column_to_rownames(var = "ensg_id")
  #   dataframe <- merge(dataframe, mapping, by = 0, all.x = TRUE)
  #
  #   up <- dplyr::filter(dataframe, FoldChange > 0)
  #   down <- dplyr::filter(dataframe, FoldChange < 0)
  #   # Print how many DE genes
  #   print(paste0("up: ", length(row.names(up)), ", down: ", length(row.names(down))))
  #
  #   up <- enrichPathway(up$entrez_id, readable = TRUE) %>% as.data.frame()
  #   if(length(row.names(up) > 0)){
  #     up$id=filename
  #     up$direction="Up"
  #   }
  #   down <- enrichPathway(down$entrez_id, readable = TRUE) %>% as.data.frame()
  #   if(length(row.names(down) > 0)){
  #     down$id=filename
  #     down$direction="Down"
  #   }
  #   total <- rbind(up, down)
  #   total$`-log10 p-adj`=(-log((total$p.adjust+1E-300), 10))
  #   names(total)[1] <- "pathwy.id"
  #   names(total)[2] <- "description"
  #
  #   # Annotate the "top" pathways
  #   # Generates three new columns to fill
  #   total[, c("Top Pathway Name", "Hierarchy")] <- NA
  #
  #   # Fills in column per row based on the first column (pathway id)
  #   for (i in c(1:length(row.names(total)))) {
  #     total[i , c("Top Pathway Name", "Hierarchy")] <- path_steps_local(total[i,1])
  #     print(i)
  #   }
  #
  #   # Save as a table
  #   write.table(total, paste0("react_pathways_",filename,".csv"), sep = ",", row.names=FALSE)
  #
  # }


  ### Hallmark
  # if (analysis == "hallmark") {
  #   database <- hallmark
  #   msigdbr_t2g = database %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
  #
  #   up <- dplyr::filter(dataframe, FoldChange > 0)
  #   down <- dplyr::filter(dataframe, FoldChange < 0)
  #
  #   # Print how many DE genes
  #   print(paste0("up: ", length(row.names(up)), ", down: ", length(row.names(down))))
  #
  #   up <- enricher(
  #     gene = rownames(dataframe %>% filter(log2FoldChange > 0)),
  #     TERM2GENE = msigdbr_t2g,
  #     universe = universe
  #   ) %>% as.data.frame()
  #
  #   if (length(rownames(up)) != 0) {
  #     up$direction <- "Up"
  #   }
  #
  #
  #   down <- enricher(
  #     gene = rownames(dataframe %>% filter(log2FoldChange < 0)),
  #     TERM2GENE = msigdbr_t2g, universe = universe
  #   ) %>% as.data.frame()
  #
  #   if (length(rownames(down)) != 0) {
  #     down$direction <- "Down"
  #   }
  #
  #   total <- rbind(up, down)
  #   total$id <- filename
  #   total$`-log10 p-adj`= -log10(total$p.adjust +1E-300)
  #
  # }

  # Return the results, which is a tibble of all the pathway results
  return(as_tibble(final_enriched_results))
}
