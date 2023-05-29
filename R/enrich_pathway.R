#' Test lists of genes for enriched pathways
#'
#' @param input_list list of data frames of DESeq2 results. The list names
#'   are used as the comparison for each dataframe (e.g. COVID vs Healthy). Data
#'   frames should have Ensembl gene IDs as rownames.
#' @param filter_input If providing list of data frames containing the
#'   unfiltered output from `DESeq2::results()`, set this to TRUE to filter for
#'   DE genes using the thresholds set by the `p_cutoff` and `fc_cutoff`
#'   arguments. When FALSE it's assumed your passing the filtered
#'   results into `input_list` and no more filtering will be done.
#' @param p_cutoff Adjusted p value cutoff, defaults to < 0.05.
#' @param fc_cutoff Absolute fold change, defaults to |FC| > 1.5.
#' @param split Boolean (TRUE); Split into up and down-regulated DEGs and do
#'   enrichment separately
#' @param analysis Default is SIGORA ("sigora"), others: ReactomePA
#'   ("reactomepa"), mSigDB Hallmark gene sets ("hallmark")
#' @param filter_results Should the output be filtered for significance? Use
#'   `1` to return the unfiltered results, or a number between 0 and 1 for a 
#'   custom p-value cutoff. If `default`, the significance cutoffs for Sigora is
#'   <0.001, and for ReactomePA or Hallmark is <0.05.
#' @param gps_repo Gene Pair Signature object for Sigora to use to test for
#'   enriched pathways. We recommend using the one which ships with Sigora,
#'   which is already loaded as "reaH".
#' @param gene_universe Optional. Set of background genes to use when testing
#'   with ReactomePA or Hallmark gene sets. For ReactomePA this must be a
#'   character vector of Entrez genes. For Hallmark, it must be Ensembl IDs.
#'
#' @return A data frame of pathway enrichment results
#' @export
#'
#' @importFrom purrr possibly
#' @importFrom clusterProfiler enricher
#' @import dplyr
#'
#' @description This function provides a simple and consistent interface to
#'   three different pathway enrichment tools: Sigora and ReactomePA (which both
#'   test for Reactome pathways), and MSigDB Hallmark gene set enrichment. The
#'   input must be a named list of data frames, which can be pre-filtered or
#'   "raw", in which case the function can filter using user-defined cutoffs.
#'   Column names are expected to comply with those output by
#'   `DESeq2::results()` function, namely: padj and log2FoldChange. Rownames are
#'   assumed to contain the input genes to be tested.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' ex_results_sigora <- enrich_pathway(
#'     deseq_example_list[c(5, 6)],
#'     gps_repo = reaH
#' )
#'
#' plot_pathways(ex_results_sigora, columns = 2)
#'
enrich_pathway <- function(
        input_list,
        filter_input = TRUE,
        p_cutoff = 0.05,
        fc_cutoff = 1.5,
        split = TRUE,
        analysis = "sigora",
        filter_results = 'default',
        gps_repo = reaH,
        gene_universe = NULL
) {

    ### Check inputs
    stopifnot(analysis %in% c("sigora", "reactomepa", "hallmark"))

    if (
        !is.list(input_list) |
        !all(unlist(lapply(input_list, is.data.frame))) |
        is.null(names(input_list))
    ) {
        stop(
            "Provide a named list of data frames of results, with the name of ",
            "each item in the list as the comparison name."
        )
    }


    # Start looping through each data frame in "input_list"
    for (i in seq_len(length(input_list))) {

        # Make sure the rownames aren't just 1:nrow(), which can happen if the
        # input data frame is a tibble
        if (all(rownames(input_list[[i]]) == seq_len(nrow(input_list[[i]])))) {
            stop(
                "The rownames of the data frame for the element with name '",
                names(input_list[i]),
                "' don't look like gene IDs!"
            )
        }

        comparison <- names(input_list[i])
        message("Comparison being analyzed: ", comparison)

        # Use the DE genes, filtering if specified
        if (filter_input) {
            stopifnot(
                c("padj", "log2FoldChange") %in% colnames(input_list[[i]])
            )
            message("\tFiltering the results using before testing for pathways")

            deseq_results <- input_list[[i]] %>%
                filter(padj < p_cutoff, abs(log2FoldChange) > log2(fc_cutoff))
        } else {
            deseq_results <- input_list[[i]]
        }

        # If splitting into up- and down-regulated DE genes
        if (split) {
            up_gns <- filter(deseq_results, log2FoldChange > 0) %>% rownames()
            dn_gns <- filter(deseq_results, log2FoldChange < 0) %>% rownames()
            message(
                "\tDEGs used: ",
                length(up_gns), " Up, ",
                length(dn_gns), " Down"
            )
        } else {
            all_gns <- rownames(deseq_results)
            message("DEGs used: ", length(all_gns))
        }

        ### Sigora analysis
        # Sigora uses gene pairs with the Reactome database, which often
        # eliminates duplicate, closely related pathways (e.g. TLR7/8/9
        # pathways). Useful for if you have a lot of DEGs that might enrich for
        # a lot of pathways, making it difficult analyze.
        if (analysis == "sigora") {
            message("\tRunning enrichment using Sigora")

            run_sigora_safely <- purrr::possibly(run_sigora)
            
            # set default pvalue cutoff
            if (filter_results == 'default') {
                filter_results <- 0.001
            }

            # Enrich with up- and down-regulated genes separately
            if (split) {
                up_results <- run_sigora_safely(
                    up_gns,
                    direction = "Up",
                    gps_repo = gps_repo,
                    pval_filter = filter_results
                )
                dn_results <- run_sigora_safely(
                    dn_gns,
                    direction = "Down",
                    gps_repo = gps_repo,
                    pval_filter = filter_results
                )
                total_results <- rbind(up_results, dn_results)

                # Or just use all DE genes for enrichment
            } else {
                total_results <- run_sigora_safely(
                    all_gns,
                    direction = "All",
                    gps_repo = gps_repo,
                    pval_filter = filter_results
                )
            }
        }


        ### ReactomePA analysis
        # Reactome uses regular over-representation analysis. Useful for finding
        # more pathways, if you have fewer DE genes or pathways. Less stringent
        # than SIGORA enrichment. ReactomePA however is slow so we can use
        # just the enricher function that is ReactomePA is built off of without
        # loading extra packages/annotaton files
        if (analysis == "reactomepa") {

            message("\tRunning enrichment using ReactomePA")
            
            # set default pvalue cutoff
            if (filter_results == 'default') {
                filter_results <- 0.05
            }

            # ReactomePA needs to map to Entrez IDs, as it doesn't take Ensembl
            # IDs. We can use SIGORA's mapping data for consistency, though it
            # maps fewer Entrez symbols than our own mapping file.
            if (split) {
                up_gns_entrez <- idmap %>%
                    filter(Ensembl.Gene.ID %in% up_gns) %>%
                    pull(EntrezGene.ID) %>%
                    unique()

                up_results <- up_gns_entrez %>%
                    enricher(
                        TERM2GENE = select(
                            reactome_database,
                            pathway_id, entrez_id
                        ),
                        TERM2NAME = select(
                            reactome_database,
                            pathway_id,
                            pathway_name
                        ),
                        universe = gene_universe,
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff = filter_results
                    ) %>%
                    as.data.frame() %>%
                    mutate(direction = "Up")

                dn_gns_entrez <- idmap %>%
                    filter(Ensembl.Gene.ID %in% dn_gns) %>%
                    pull(EntrezGene.ID) %>%
                    unique()

                dn_results <- dn_gns_entrez %>%
                    enricher(
                        TERM2GENE = select(
                            reactome_database,
                            pathway_id, entrez_id
                        ),
                        TERM2NAME = select(
                            reactome_database,
                            pathway_id,
                            pathway_name
                        ),
                        universe = gene_universe,
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff = filter_results
                    ) %>%
                    as.data.frame() %>%
                    mutate(direction = "Down")

                total_results <- rbind(up_results, dn_results)

            } else {
                all_gns_entrez <- idmap %>%
                    filter(Ensembl.Gene.ID %in% all_gns) %>%
                    pull(EntrezGene.ID) %>%
                    unique()

                total_results <- all_gns_entrez %>%
                    enricher(
                        TERM2GENE = select(
                            reactome_database,
                            pathway_id, entrez_id
                        ),
                        TERM2NAME = select(
                            reactome_database,
                            pathway_id,
                            pathway_name
                        ),
                        universe = gene_universe,
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff = filter_results
                    ) %>%
                    as.data.frame() %>%
                    mutate(direction = "All")
            }

            # Map the entrez ids to hgnc symbols
            hgnc_gene_list <- c()
            for (r in seq_len(nrow(total_results))) {
                genelist <- total_results[r, "geneID"]
                genelist <- str_split(genelist, "/") %>% unlist()
                hgnc_genes <- mapping_file %>%
                    filter(entrez_id %in% genelist) %>%
                    .$gene_name
                hgnc_genes <- paste(hgnc_genes, collapse = ";")
                hgnc_gene_list <- c(hgnc_gene_list, hgnc_genes)
            }
            total_results$genes <- hgnc_gene_list

            # Clean up the column names
            total_results <- total_results %>%
                separate(
                    GeneRatio,
                    into = c("num_candidate_genes", "num_bg_genes"),
                    sep = "/"
                ) %>%
                # Replace "/" with ";" for consistency in "candidate_genes"
                # column
                mutate(
                    across(c(num_candidate_genes, num_bg_genes), as.numeric),
                    gene_ratio = num_candidate_genes / num_bg_genes,
                ) %>%
                select(
                    "pathway_id" = ID,
                    "pathway_description" = Description,
                    direction,
                    "p_value" = pvalue,
                    "p_value_adjusted" = p.adjust,
                    genes,
                    num_candidate_genes,
                    num_bg_genes,
                    gene_ratio
                )
        }

        ### Hallmark enrichment from mSigDB
        # The Hallmark gene sets summarize and represent 50 specific
        # well-defined biological states or processes and display coherent
        # expression. They are a useful tool for looking at broad mechanisms,
        # and can "validate" enrichment findings from SIGORA or ReactomePA if
        # similar mechanisms appear.
        if (analysis == "hallmark") {

            message("\tRunning enrichment using Hallmark")
            
            if (filter_results == 'default') {
                filter_results <- 0.05
            }

            if (split) {
                up_results <- enricher(
                    up_gns,
                    TERM2GENE = msigdbr_t2g,
                    universe = gene_universe,
                    pvalueCutoff = filter_results
                ) %>%
                    as.data.frame() %>%
                    mutate(direction = "Up")

                dn_results <- enricher(
                    dn_gns,
                    TERM2GENE = msigdbr_t2g,
                    universe = gene_universe,
                    pvalueCutoff = filter_results
                ) %>%
                    as.data.frame() %>%
                    mutate(direction = "Down")

                total_results <- rbind(up_results, dn_results)

            } else {
                total_results <- enricher(
                    all_gns,
                    TERM2GENE = msigdbr_t2g,
                    universe = gene_universe,
                    pvalueCutoff = filter_results
                ) %>%
                    as.data.frame() %>%
                    mutate(direction = "All")
            }

            # annotate the genes as hgnc symbols and add them into the table
            hgnc_gene_list <- c()
            for (r in seq_len(nrow(total_results))) {
                genelist <- total_results[r, "geneID"]
                genelist <- str_split(genelist, "/") %>% unlist()
                hgnc_genes <- mapping_file %>%
                    filter(ensg_id %in% genelist) %>%
                    .$gene_name
                hgnc_genes <- paste(hgnc_genes, collapse = ";")
                hgnc_gene_list <- c(hgnc_gene_list, hgnc_genes)
            }
            total_results$genes <- hgnc_gene_list

            total_results <- total_results %>%
                separate(
                    GeneRatio,
                    into = c("num_candidate_genes", "num_bg_genes"),
                    sep = "/"
                ) %>%
                mutate(
                    across(c(num_candidate_genes, num_bg_genes), as.numeric),
                    gene_ratio = num_candidate_genes / num_bg_genes
                ) %>%
                select(
                    "pathway_id" = ID,
                    "pathway_description" = Description,
                    direction,
                    "p_value" = pvalue,
                    "p_value_adjusted" = p.adjust,
                    genes,
                    num_candidate_genes,
                    num_bg_genes,
                    gene_ratio
                )
        }

        # Annotate the  enrichment results by adding top pathways and
        # comparisons
        total_results_annotated <- left_join(
            total_results,
            top_pathways_more %>% dplyr::select(pathway_id, top_pathways),
            by = "pathway_id"
        )

        total_results_annotated$comparison <- comparison
        total_results_annotated$total_genes <- nrow(deseq_results)

        message(
            "\tDone! Enriched pathways: ", nrow(total_results_annotated), "\n"
        )

        # This adds the next dataframe to the previous one
        if (i == 1) {
            final_enriched_results <- total_results_annotated
        } else {
            final_enriched_results <-
                rbind(final_enriched_results, total_results_annotated)
        }
    }

    # Order the comparisons in the same order as entered (for graphing purposes
    # later)
    final_enriched_results$comparison <- factor(
        final_enriched_results$comparison,
        levels = c(names(input_list))
    )

    # Return the results, which is a tibble of all the pathway results
    return(as_tibble(final_enriched_results))
}
