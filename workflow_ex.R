# Load the package --------------------------------------------------------

devtools::load_all(".")


# Example 1 ---------------------------------------------------------------

# |- get_pathway_distances and create_foundation --------------------------

# dist_data <- get_pathway_distances(
#   pathway_data = sigora_database,
#   dist_method = "jaccard"
# )

starting_pathways <- create_foundation(
  mat = pathway_distances_jaccard,
  # prop_to_keep = 0.004
  max_distance = 0.8
)


# |- create_pathnet -------------------------------------------------------

my_pathway_network <- create_pathnet(
  sigora_result = sigora_examples[[1]],
  foundation = starting_pathways,
  trim = TRUE,
  trim_order = 1
)


# |- plot -----------------------------------------------------------------

# Static plot with ggraph
pathnet_ggraph(
  my_pathway_network,
  label_prop = 0.1,
  node_label_size = 4,
  node_label_overlaps = 8,
  seg_colour = "red"
)


# visNetwork --------------------------------------------------------------

# Interactive plot with visNetwork
pathnet_visNetwork(my_pathway_network)


# Example 2 ---------------------------------------------------------------

# Here we'll restrict the data using the pathways/genes from the Sigora results,
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
  sigora_result = sigora_example[[2]],
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



# eruption ----------------------------------------------------------------

eruption(deseq_results = deseq_example_list[[1]])


# enrich_pathway ----------------------------------------------------------

# Enrich with SIGORA
enriched_results_sigora <- enrich_pathway(
  deseq_example_list[c(5, 6)],
  gps_repo = reaH,
  filter = TRUE
)

# Enrich with ReactomePA
enriched_results_rpa <- enrich_pathway(
  deseq_example_list[c(5, 6)],
  analysis = "reactomepa",
  filter = TRUE
)

# enrich with Hallmark gene sets
enriched_results_hm <- enrich_pathway(
  deseq_example_list[c(5, 6)],
  analysis = "hallmark",
  filter = TRUE
)


# plot_pathways -----------------------------------------------------------

# plot the different outputs
plot_pathways(enriched_results_sigora, columns = 2)
plot_pathways(enriched_results_rpa, columns = 2)
plot_pathways(enriched_results_hm)

# if you only want "immune system" pathways
plot_pathways(enriched_results_sigora,
              specific_top_pathways = 'Immune System')

# if you want to include gene ratio
plot_pathways(enriched_results_sigora,
              specific_top_pathways = 'Immune System',
              include_gene_ratio = TRUE)

# plot_fold_change --------------------------------------------------------

# plot the fold change for the example DESeq2 results for the pathway
# 'Generation of second messenger molecules'
plot_fold_change(deseq_example_list, path_name = 'Generation of second messenger molecules')

# Separate out the different conditions so that the names don't take up so much space
deseq_example_list_renamed <- deseq_example_list
names(deseq_example_list_renamed) <- c('PosT1', 'PosT2', 'NegT1', 'NegT2', 'Pos', 'Neg')
plot_fold_change(deseq_example_list_renamed,
                 path_name = 'Generation of second messenger molecules',
                 col_split = c('Pos', 'Pos', 'Neg', 'Neg', 'Time', 'Time'),
                 col_angle = 0)

# I just want to plot genes I like, and I also want to see them in log2FC,
# I also want to swap the column and rows
plot_fold_change(deseq_example_list_renamed,
                 manual_title = 'Genes I like',
                 genes_to_plot = c('CD4', 'CD8A','CD8B', 'CD28', 'ZAP70'),
                 gene_format = 'hgnc',
                 col_split = c('Pos', 'Pos', 'Neg', 'Neg', 'Time', 'Time'),
                 log2_foldchange = TRUE,
                 col_angle = 45,
                 cluster_rows = FALSE,
                 cluster_columns = TRUE,
                 invert = TRUE)
