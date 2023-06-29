#' Test lists of genes for enriched pathways
#'
#' @param input_list list of data frames of DESeq2 results. The list names
#'   are used as the comparison for each dataframe (e.g. COVID vs Healthy). Data
#'   frames must have Ensembl gene IDs as the rownames.
#' @param filter_input If providing list of data frames containing the
#'   unfiltered output from `DESeq2::results()`, set this to TRUE to filter for
#'   DE genes using the thresholds set by the `p_cutoff` and `fc_cutoff`
#'   arguments. When FALSE it's assumed your passing the filtered
#'   results into `input_list` and no more filtering will be done.
#' @param p_cutoff Adjusted p value cutoff, defaults to 0.05.
#' @param fc_cutoff Minimum absolute fold change, defaults to 1.5.
#' @param split Boolean (TRUE); Split into up and down-regulated DEGs using the
#'   column "log2FoldChange", and do enrichment separately
#' @param analysis Default is "sigora", but can also be "reactomepa" or
#'   "hallmark"
#' @param filter_results Should the output be filtered for significance? Use
#'   `1` to return the unfiltered results, or any number less than 1 for a
#'   custom p-value cutoff. If `default`, the significance cutoff for Sigora
#'   is 0.001, and for ReactomePA or Hallmark is 0.05.
#' @param gps_repo Only applies to `analysis = "sigora"`. Gene Pair Signature
#'   object for Sigora to use to test for enriched pathways. We recommend using
#'   the one which ships with Sigora, which is already loaded as "reaH".
#' @param gene_universe Only applies when `analysis` is "reactomepa" or
#'   "hallmark". The set of background genes to use when testing with ReactomePA
#'   or Hallmark gene sets. For ReactomePA this must be a character vector of
#'   Entrez genes. For Hallmark, it must be Ensembl IDs.
#'
#' @return A data frame of pathway enrichment results for all input comparisons
#' @export
#'
#' @import dplyr
#' @importFrom purrr possibly imap_dfr
#' @importFrom clusterProfiler enricher
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
#' enrich_pathway(
#'     input_list = deseq_example_list[1],
#'     filter_input = TRUE,
#'     split = TRUE,
#'     analysis = "sigora",
#'     gps_repo = reaH,
#'     filter_results = "default",
#' )
#'
enrich_pathway <- function(
        input_list,
        filter_input = TRUE,
        p_cutoff = 0.05,
        fc_cutoff = 1.5,
        split = TRUE,
        analysis = "sigora",
        filter_results = "default",
        gps_repo = reaH,
        gene_universe = NULL
) {

    ### Check inputs
    stopifnot(analysis %in% c("sigora", "reactomepa", "hallmark"))

    stopifnot(
        "Provide a named list of data frames of results, with the name
        of each item in the list as the comparison name."  = {
            is.list(input_list)
            !is.null(names(input_list))
            all(unlist(lapply(input_list, is.data.frame)))
        }
    )

    ### Iterate through each element of "input_list"
    result_list <- imap(input_list, function(x, comparison) {

        # Basic checks before continuing
        stopifnot(
            "Elements of 'input_list' should be named" = !is.null(comparison)
        )

        if (all(rownames(x) == seq_len(nrow(x)))) {
            stop(
                "The rownames of the data frame for the element with name '",
                comparison, "' don't look like gene IDs!"
            )
        }

        message("Comparison being analyzed: ", comparison)

        # Filter the input genes if specified
        deseq_results <-
            if (filter_input) {
                stopifnot(
                    c("padj", "log2FoldChange") %in% colnames(x)
                )
                message("\tFiltering the results before testing...")
                filter(
                    x,
                    padj < p_cutoff,
                    abs(log2FoldChange) > log2(fc_cutoff)
                )
            } else {
                x
            }

        # Turn the input into a list of gene IDs, split by direction or not
        if (split) {
            prepped_genes <- list(
                "Up"   = rownames(filter(deseq_results, log2FoldChange > 0)),
                "Down" = rownames(filter(deseq_results, log2FoldChange < 0))
            )
            message(
                "\tDEGs used: ",
                length(prepped_genes$Up), " Up, ",
                length(prepped_genes$Down), " Down..."
            )
        } else {
            prepped_genes <- list("All" = rownames(deseq_results))
            message("\tDEGs used: ", length(prepped_genes$All), "...")
        }

        ### Sigora ###
        if (analysis == "sigora") {
            message("\tRunning enrichment using Sigora...")

            run_sigora_safely <- possibly(run_sigora)

            result_final <- imap_dfr(
                .x  = prepped_genes,
                .id = "direction",
                function(y, direction) {
                    run_sigora_safely(
                        enrich_genes = y,
                        gps_repo = gps_repo,
                        pval_filter = ifelse(
                            filter_results == "default",
                            0.001,
                            filter_results
                        )
                    )
                }
            )

            result_final$total_genes <- nrow(deseq_results)

            message(
                "\tDone! Found ",
                nrow(result_final),
                " enriched pathways.\n"
            )
            return(result_final)
        }


        ### ReactomePA or Hallmark ###
        if (analysis %in% c("reactomepa", "hallmark")) {

            ### ReactomePA ###
            if (analysis == "reactomepa") {
                message("\tRunning enrichment using ReactomePA")

                rpa_hall_result <- imap_dfr(
                    .x  = prepped_genes,
                    .id = "direction",
                    function(y, direction) {

                        genes_entrez <- idmap %>%
                            filter(Ensembl.Gene.ID %in% y) %>%
                            pull(EntrezGene.ID) %>%
                            unique()

                        enricher(
                            genes_entrez,
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
                            pvalueCutoff = ifelse(
                                filter_results == "default",
                                0.05,
                                filter_results
                            )
                        ) %>% as_tibble()
                    }
                ) %>%
                    mutate(geneID = as.character(geneID)) %>%
                    separate_longer_delim(geneID, delim = "/") %>%
                    left_join(mapping_file, by = c("geneID" = "entrez_id")) %>%
                    select(-any_of(c("geneID", "entrez_id", "ensg_id"))) %>%
                    group_by(ID) %>%
                    mutate(genes = paste(gene_name, collapse = ";")) %>%
                    ungroup() %>%
                    select(-gene_name) %>%
                    distinct()
            }

            ### Hallmark ###
            if (analysis == "hallmark") {
                message("\tRunning enrichment using Hallmark...")

                rpa_hall_result <- imap_dfr(
                    .x  = prepped_genes,
                    .id = "direction",
                    function(y, direction) {

                        enricher(
                            y,
                            TERM2GENE = msigdbr_t2g,
                            universe = gene_universe,
                            pvalueCutoff = ifelse(
                                filter_results == "default",
                                0.05,
                                filter_results
                            )
                        ) %>% as_tibble()
                    }
                ) %>%
                    separate_longer_delim(geneID, delim = "/") %>%
                    left_join(mapping_file, by = c("geneID" = "ensg_id")) %>%
                    select(-any_of(c("geneID", "entrez_id", "ensg_id"))) %>%
                    group_by(ID) %>%
                    mutate(genes = paste(gene_name, collapse = ";")) %>%
                    ungroup() %>%
                    select(-gene_name) %>%
                    distinct()
            }

            ### ReactomePA or Hallmark ###
            result_final <- rpa_hall_result %>%
                separate_wider_delim(
                    cols = GeneRatio,
                    delim = "/",
                    names = c("num_candidate_genes", "num_bg_genes")
                ) %>%
                mutate(
                    across(c(num_candidate_genes, num_bg_genes), as.numeric),
                    gene_ratio = num_candidate_genes / num_bg_genes,
                    total_genes = nrow(deseq_results)
                ) %>%
                select(
                    direction,
                    "pathway_id" = ID,
                    "pathway_description" = Description,
                    "p_value" = pvalue,
                    "p_value_adjusted" = p.adjust,
                    genes,
                    num_candidate_genes,
                    num_bg_genes,
                    gene_ratio,
                    total_genes
                )

            message(
                "\tDone! Found ",
                nrow(result_final),
                " enriched pathways.\n"
            )
            return(result_final)
        }
    })

    results_all_comparisons <- result_list %>%
        bind_rows(.id = "comparison") %>%
        left_join(
            .,
            select(top_pathways_more, pathway_id, top_pathways),
            by = "pathway_id"
        ) %>%
        as_tibble()

    return(results_all_comparisons)
}
