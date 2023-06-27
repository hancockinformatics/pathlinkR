#' INTERNAL run_sigora
#'
#' @param enrich_genes Vector of genes to enrich
#' @param gps_repo GPS object to use for testing pathways
#' @param pval_filter Desired threshold for filtering results
#'
#' @return Data frame of results from Sigora
#'
#' @import dplyr
#' @import sigora
#' @importFrom purrr map imap
#'
#' @description Internal wrapper function to run Sigora and return the results
#'   with desired columns
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
run_sigora <- function(
        enrich_genes,
        gps_repo,
        pval_filter = NA
) {

    # Run SIGORA based on default settings (GPSrepo = reaH, level = 4)
    invisible(capture.output(
        sigora_result_1 <- sigora(
            GPSrepo = gps_repo,
            level = 4,
            markers = TRUE,
            queryList = enrich_genes
        )
    ))

    sigora_result_2 <-
        if (is.na(pval_filter)) {
            sigora_result_1$summary_results
        } else {
            filter(sigora_result_1$summary_results, Bonferroni < pval_filter)
        }

    # Start by calculating n
    n_genes <- length(enrich_genes[enrich_genes %in% idmap$Ensembl.Gene.ID])

    # Get the DE genes that were enriched for each pathway
    sigora_detailed_list_1 <- sigora_result_1$detailed_results %>%
        as_tibble() %>%
        select(pathway, contains("gene")) %>%
        split(x = ., f = .$pathway)

    sigora_detailed_list_2 <- sigora_detailed_list_1 %>% imap(
        ~tibble(
            pathwy.id = .y,
            EntrezGene.ID = unique(c(pull(.x, gene1), pull(.x, gene2)))
        ) %>%
            left_join(idmap, by = "EntrezGene.ID", multiple = "all") %>%
            filter(Ensembl.Gene.ID %in% enrich_genes, Symbol != "^$") %>%
            select(pathwy.id, Symbol) %>%
            mutate(across(everything(), as.character)) %>%
            distinct() %>%
            arrange(pathwy.id, Symbol)
    ) %>%
        bind_rows() %>%
        group_by(pathwy.id) %>%
        summarize(genes = paste0(Symbol, collapse = ";")) %>%
        mutate(num_candidate_genes = 1 + str_count(genes, ";")) %>%
        ungroup() %>%
        mutate(gene_ratio = num_candidate_genes / n_genes)

    left_join(sigora_result_2, sigora_detailed_list_2, by = "pathwy.id") %>%
        mutate(num_bg_genes = n_genes) %>%
        select(
            "pathway_id" = pathwy.id,
            "pathway_description" = description,
            # direction,
            "p_value" = pvalues,
            "p_value_adjusted" = Bonferroni,
            genes,
            num_candidate_genes,
            num_bg_genes,
            gene_ratio
        )
}
