#' Test lists of DE genes for enriched pathways
#'
#' @param deseq_result_list list of data frames of DESeq2 results. The list names
#'   are used as the comparison for each dataframe (e.g. COVID vs Healthy). Data
#'   frames should have Ensembl gene IDs as rownames.
#' @param p_cutoff Adjusted p value cutoff, defaults to < 0.05.
#' @param fc_cutoff Absolute fold change, defaults to |FC| > 1.5.
#' @param split Boolean (TRUE); Split into up and down-regulated DEGs and do
#'   enrichment separately
#' @param analysis Default is SIGORA ("sigora"), others: ReactomePA
#'   ("reactomepa"), mSigDB Hallmark gene sets ("hallmark")
#' @param gps_repo Gene Pair Signature object for Sigora to use to test for
#'   enriched pathways. We recommend using the one which ships with Sigora,
#'   which can be loaded via `data(reaH, package = "sigora")`.
#'
#' @return A data frame of pathway enrichment results
#' @export
#'
#' @import ReactomePA
#' @import clusterProfiler
#' @import dplyr
#'
enrich_pathway <- function(
    deseq_result_list,
    p_cutoff = 0.05,
    fc_cutoff = 1.5,
    split = TRUE,
    analysis = "sigora",
    gps_repo
) {

  # TODO
  # 1. Define what the "universe" is for ReactomePA and Hallmark
  # 2. Include just using a vector of genes or data frame rather than a list of
  #    DESeq2 results
  # 3. Should the candidate genes be in Ensembl IDs or HGNC symbols?
  # 4. Need to change top pathway mapping file once agreed upon

  if (!is.list(deseq_result_list)) {
    message("Provide a list of dataframes of DESeq2 results with the name of ",
            "each item in the list as the comparison name.")
    return()

  } else {
    for (i in 1:length(deseq_result_list)) {
      comparison <- names(deseq_result_list[i])
      message(paste0("Comparison being analyzed: ", comparison))

      # Use the DE genes
      deseq_results <- deseq_result_list[i][[1]] %>%
        filter(padj < p_cutoff, abs(log2FoldChange) > log2(fc_cutoff))

      # If splitting into up and down-regulated DE genes
      if (split) {
        up_gns <- filter(deseq_results, log2FoldChange > 0) %>% rownames()
        dn_gns <- filter(deseq_results, log2FoldChange < 0) %>% rownames()
        message(paste0(
          "DEGs used: Up: ", length(up_gns),
          ", Down: ", length(dn_gns)
        ))
      } else {
        all_gns <- rownames(deseq_results)
        message(paste0("DEGs used: ", length(all_gns))) # Print how many DE genes used
      }

      ### SIGORA analysis
      # SIGORA uses gene pairs with the Reactome database, which often
      # eliminates duplicate, closely related pathways (e.g. TLR7/8/9 pathways).
      # Useful for if you have a lot of DEGs that might enrich for a lot of
      # pathways, making it difficult analyze.

      if (analysis == "sigora") {
        message("Enrichment using SIGORA")

        run_sigora_safely <- possibly(run_sigora)

        # Enrich with up- and down-regulated genes separately
        if (split) {
          up_results <- run_sigora_safely(up_gns, direction = "Up", gps_repo = gps_repo)
          dn_results <- run_sigora_safely(dn_gns, direction = "Down", gps_repo = gps_repo)
          total_results <- rbind(up_results, dn_results)

        # Or just use all DE genes for enrichment
        } else {
          total_results <- run_sigora_safely(all_gns, direction = "All")
        }
      }


      ### ReactomePA analysis
      # Reactome uses regular over-representation analysis
      # Useful for finding more pathways, if you have fewer DE genes or pathways
      # Less stringent than SIGORA enrichment

      if (analysis == "reactomepa") {
        # ReactomePA needs to map to entrez_ids, does not take ensembl ids
        # can use SIGORA's mapping for consistency, but maps fewer entrez symbols
        if (split) {
          up_gns_entrez <- idmap %>%
            filter(Ensembl.Gene.ID %in% up_gns) %>%
            .$EntrezGene.ID %>%
            unique()

          up_results <- enrichPathway(up_gns_entrez, readable = TRUE) %>%
            as.data.frame() %>%
            mutate(direction = "Up")

          dn_gns_entrez <- idmap %>%
            filter(Ensembl.Gene.ID %in% dn_gns) %>%
            .$EntrezGene.ID %>%
            unique()

          dn_results <- enrichPathway(dn_gns_entrez, readable = TRUE) %>%
            as.data.frame() %>%
            mutate(direction = "Down")

          total_results <- rbind(up_results, dn_results)

        } else {
          all_gns_entrez <- idmap %>%
            filter(Ensembl.Gene.ID %in% all_gns) %>%
            .$EntrezGene.ID %>%
            unique()

          total_results <- enrichPathway(all_gns_entrez, readable = TRUE) %>%
            as.data.frame() %>%
            mutate(direction = "All")
        }

        # clean up the column names
        total_results <- total_results %>%
          separate(
            GeneRatio,
            into = c("num_candidate_genes", "num_bg_genes"),
            sep = "/"
          ) %>%
          transmute(
            pathway_id = ID,
            pathway_description = Description,
            direction,
            p_value = pvalue,
            p_value_adjusted = p.adjust,
            q_value = qvalue,
            candidate_genes = geneID,
            num_candidate_genes = as.numeric(num_candidate_genes),
            num_bg_genes = as.numeric(num_bg_genes),
            gene_ratio = num_candidate_genes/num_bg_genes
          ) %>%
          # Output genes are separated by '/', replace with ';' for consistency
          mutate_at(
            .vars = "candidate_genes",
            .funs = gsub,
            pattern = "/",
            replacement = ";"
          )
      }

      ### Hallmark enrichment from mSigDB
      # The Hallmark gene sets summarize and represent 50 specific well-defined biological
      # states or processes and display coherent expression. They are a useful tool
      # for looking at broad mechanisms, and can "validate" enrichment findings from
      # SIGORA or ReactomePA if similar mechanisms appear.

      if (analysis == "hallmark") {
        if (split) {
          up_results <- enricher(
            up_gns,
            TERM2GENE = msigdbr_t2g
          ) %>%
            as.data.frame() %>%
            mutate(direction = "Up")

          dn_results <- enricher(
            dn_gns,
            TERM2GENE = msigdbr_t2g
          ) %>% as.data.frame() %>%
            mutate(direction = "Down")

          total_results <- rbind(up_results, dn_results)

        } else {
          total_results <- enricher(
            all_gns,
            TERM2GENE = msigdbr_t2g
          ) %>% as.data.frame() %>%
            mutate(direction = "All")
        }

        total_results <- total_results %>%
          separate(
            GeneRatio,
            into = c("num_candidate_genes", "num_bg_genes"),
            sep = "/"
          ) %>%
          transmute(
            pathway_id = ID,
            pathway_description = Description,
            direction,
            p_value = pvalue,
            p_value_adjusted = p.adjust,
            q_value = qvalue,
            candidate_genes = geneID,
            num_candidate_genes = as.numeric(num_candidate_genes),
            num_bg_genes = as.numeric(num_bg_genes),
            gene_ratio = num_candidate_genes/num_bg_genes
          ) %>%
          # Output genes are separated by '/', replace with ';' for consistency
          mutate_at(
            .vars = "candidate_genes",
            .funs = gsub,
            pattern = "/",
            replacement = ";"
          )

      }

      # Annotate the  enrichment results by adding top pathway and comparisons
      total_results_annotated <- left_join(
        total_results,
        top_pathways_more %>% select(pathway_id, top_pathways),
        by = "pathway_id"
      )

      total_results_annotated$comparison <- comparison
      total_results_annotated$total_genes <- nrow(deseq_results)
      message(paste0(
        "Done! Enriched pathways: ", nrow(total_results_annotated)
      ))

      # This adds the next dataframe to the previous one
      if (i == 1) {
        final_enriched_results <- total_results_annotated
      } else {
        final_enriched_results <-
          rbind(final_enriched_results, total_results_annotated)
      }
    }
  }

  # Order the comparisons in the same order as entered (for graphing purposes
  # later)
  final_enriched_results$comparison <- factor(
    final_enriched_results$comparison,
    levels = c(names(deseq_result_list))
  )

  # Return the results, which is a tibble of all the pathway results
  return(as_tibble(final_enriched_results))
}
