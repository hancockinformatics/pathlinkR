# Load the package --------------------------------------------------------

devtools::load_all(".")


# Example 1 ---------------------------------------------------------------

# |- get_pathway_distances and create_foundation --------------------------

dist_data <- get_pathway_distances(
  pathway_data = sigora_database,
  dist_method = "jaccard"
)

starting_pathways <- create_foundation(
  mat = pathway_distances_jaccard,
  prop_to_keep = 0.004
  # max_distance = 0.67
)


# |- create_pathnet -------------------------------------------------------

my_pathway_network <- create_pathnet(
  sigora_result = sigora_example_1,
  foundation = starting_pathways,
  trim = TRUE,
  trim_order = 2
)


# |- plot -----------------------------------------------------------------

plot_pathnet(my_pathway_network, node_label_overlaps = 8)


# Example 2 ---------------------------------------------------------------

# Here we'll restrict the data using the pathways/genes from the Sigora results,
# instead of using all pathways and layering our results on top.
# Start by splitting the sigora genes, and adding the needed annotations from
# sigora_database
candidate_data <- sigora_example_2 %>%
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
  sigora_result = sigora_example_2,
  foundation = candidate_starting_pathways,
  trim = FALSE
)

plot_pathnet(candidates_as_network) +
  facet_nodes(~direction, scales = "free", nrow = 2) +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.text.align = 0
  )
