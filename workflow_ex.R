# Load the package --------------------------------------------------------

devtools::load_all(".")


# eruption ----------------------------------------------------------------

eruption(deseq_results = deseq_example_list[[1]])


# plot_fold_change --------------------------------------------------------

# Plot the fold change for the example DESeq2 results for the pathway
# "Generation of second messenger molecules"
plot_fold_change(
    deseq_example_list,
    path_name = "Generation of second messenger molecules"
)

# Or I want to plot the Immunoregulatory interactions pathway, but the name is
# so long. Set a manual pathway name and use the pathway id
plot_fold_change(
    deseq_example_list,
    path_id = "R-HSA-198933",
    manual_title = "Immunoregulatory interactions"
)

# Separate out the different conditions so that the names don"t take up so much
# space
deseq_example_list_renamed <- purrr::set_names(
    deseq_example_list,
    c("Pos", "Neg")
)


plot_fold_change(
    deseq_example_list_renamed,
    path_name = "Generation of second messenger molecules",
    col_split = c("Pos", "Neg"),
    col_angle = 0
)

# I just want to plot genes I like, see them in log2FC, and swap the columns and
# rows
plot_fold_change(
    deseq_example_list_renamed,
    manual_title = "Genes I like",
    genes_to_plot = c("CD4", "CD8A","CD8B", "CD28", "ZAP70"),
    gene_format = "hgnc",
    col_split = c("Pos", "Neg"),
    log2_foldchange = TRUE,
    col_angle = 45,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    invert = TRUE
)


# PPI networks ------------------------------------------------------------

ex_de_genes <- deseq_example_list[[1]] %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))

ex_network <- ppi_build_network(
    df = ex_de_genes,
    col = "gene",
    order = "zero",
    hub_measure = "hubscore",
    ppi_data = innatedb_exp
)

ppi_plot_network(
    ex_network,
    fill_column = log2FoldChange,
    fill_type = "fold_change",
    layout = "lgl",
    label = TRUE,
    label_column = gene_name,
    label_filter = 5,
    legend = TRUE
)


# enrich_pathway ----------------------------------------------------------

# Enrich with SIGORA
enriched_results_sigora <- enrich_pathway(
    input_list = deseq_example_list,
    filter_input = TRUE,
    split = TRUE,
    analysis = "sigora",
    gps_repo = reaH,
    filter_results = "default",
)

# If we only have one data frame, that's already been filtered
enriched_results_sigora_2 <- enrich_pathway(
    list("one" = ex_de_genes),
    gps_repo = reaH,
    filter_input = FALSE,
    filter_results = "default"
)

# If we want all pathways (even non-significant ones, will be slower)
# enriched_results_sigora_3 <- enrich_pathway(
#     list("one" = test_de_genes),
#     gps_repo = reaH,
#     filter_input = FALSE,
#     filter_results = 1
# )

# Enrich with ReactomePA
enriched_results_rpa <- enrich_pathway(
    deseq_example_list,
    filter_input = TRUE,
    analysis = "reactomepa"
)

# Enrich with Hallmark gene sets
enriched_results_hm <- enrich_pathway(
    deseq_example_list,
    analysis = "hallmark",
    filter_input = TRUE
)


# plot_pathways -----------------------------------------------------------

# Plot the different outputs
plot_pathways(enriched_results_sigora, columns = 2)
plot_pathways(enriched_results_rpa, columns = 2)
plot_pathways(enriched_results_hm)

# If you only want "Immune System" pathways
plot_pathways(
    enriched_results_sigora,
    specific_top_pathways = "Immune System"
)

# If you want to include gene ratio
plot_pathways(
    enriched_results_sigora,
    specific_top_pathways = "Immune System",
    include_gene_ratio = TRUE
)

# If you want change up the format (change names, change label side)
plot_pathways(
    enriched_results_sigora,
    specific_top_pathways = "Immune System",
    pathway_position = "left",
    new_group_names = c("Pos", "Neg")
)


# Pathway networks v1 -----------------------------------------------------

# dist_data <- get_pathway_distances(
#   pathway_data = sigora_database,
#   dist_method = "jaccard"
# )

starting_pathways <- create_foundation(
    mat = pathway_distances_jaccard,
    # prop_to_keep = 0.004
    max_distance = 0.8
)

my_pathway_network <- create_pathnet(
    sigora_result = filter(sigora_examples, comparison == "COVID Pos Over Time"),
    foundation = starting_pathways,
    trim = TRUE,
    trim_order = 1
)

# Static plot with ggraph
pathnet_ggraph(
    my_pathway_network,
    label_prop = 0.1,
    node_label_size = 4,
    node_label_overlaps = 8,
    seg_colour = "red"
)

# Interactive plot with visNetwork
pathnet_visNetwork(my_pathway_network)


# Pathway networks v2 -----------------------------------------------------

# Here we"ll restrict the data using the pathways/genes from the Sigora results,
# instead of using all pathways and layering our results on top.
# Start by splitting the sigora genes, and adding the needed annotations from
# sigora_database
candidate_data <- sigora_examples[[2]] %>%
    select(pathway_id, genes) %>%
    separate_rows(genes, sep = ";") %>%
    left_join(
        sigora_database,
        by = c("pathway_id", "genes" = "Symbol"),
        multiple = "all"
    ) %>%
    relocate(pathway_id, Ensembl.Gene.ID, "Symbol" = genes, pathway_name) %>%
    distinct()

# Now that we have a smaller table in the same format as sigora_database, we can
# construct our own matrix of pathway distances
candidate_dist_data <- get_pathway_distances(
    pathway_data = candidate_data,
    dist_method = "jaccard"
)

candidate_starting_pathways <- create_foundation(
    mat = candidate_dist_data,
    max_distance = 0.9
)

candidates_as_network <- create_pathnet(
    sigora_result = sigora_examples[[2]],
    foundation = candidate_starting_pathways,
    trim = FALSE
)

pathnet_visNetwork(candidates_as_network)

pathnet_ggraph(candidates_as_network) +
    facet_nodes(~direction, scales = "free", nrow = 2) +
    theme_bw(base_size = 18) +
    theme(
        axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text.align = 0
    )
