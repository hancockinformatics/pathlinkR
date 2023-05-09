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


# enrich_pathway ----------------------------------------------------------

# test_list <- list('Time 1' = deseq_example_1, 'Time 2' = deseq_example_2)
enriched_results <- enrich_pathway(deseq_example_list[c(5,6)])


# plot_pathways -----------------------------------------------------------

# enriched_results <-
#   enrich_pathway(deseq_result_list = deseq_example_list[c(5, 6)], split = FALSE)
# plot_pathways(enriched_results, columns = 2)
