# Load packages & example data --------------------------------------------

library(ggraph)
devtools::load_all(".")


top_pathway_colours <- set_names(
  RColorBrewer::brewer.pal("Set2", n = 7),
  c("Cell Process",
    "Cell Replication",
    "Signaling",
    "Tissue Function",
    "Immune and Hemostasis",
    "Metabolism",
    "Gene Expression")
)


# Example 1 ---------------------------------------------------------------

# |- get_pathway_distances and create_foundation --------------------------

# dist_data <- get_pathway_distances(
#   pathway_data = sigora_database,
#   dist_method = "jaccard"
# )

starting_pathways <- create_foundation(
  mat = pathway_distances_jaccard,
  max_distance = 0.67
)


# |- create_pathnet -------------------------------------------------------

my_pathway_network <- create_pathnet(
  sigora_result = sigora_example_3,
  foundation = starting_pathways,
  trim = TRUE,
  trim_order = 1
)


# |- plot -----------------------------------------------------------------

ggraph(pathways_as_network, layout = "nicely") +
  geom_edge_link(aes(edge_width = log10(similarity)), alpha = 0.67) +
  geom_node_point(
    aes(size = -log10(bonferroni), fill = node_fill, colour = grouped_pathway),
    pch = 21,
    stroke = 1.5
  ) +
  geom_node_label(
    aes(label = pathway_name),
    repel = TRUE,
    alpha = 0.67,
    min.segment.length = 0,
    max.overlaps = 6
  ) +
  scale_edge_width(range = c(0.33, 3), name = "Similarity") +
  scale_size_continuous(
    labels = scales::label_math(expr = 10^-~.x),
    range = c(4, 8)
  ) +
  scale_colour_manual(values = top_pathway_colours) +
  scale_fill_manual(values = top_pathway_colours, na.value = "white", guide = NULL) +
  scale_shape_manual(values = c("up" = 21), na.value = 16) +
  labs(
    size = "Bonferroni\np-value",
    colour = "Pathway type"
  ) +
  theme_void(base_size = 18) +
  guides(
    colour = guide_legend(override.aes = list(size = 6)),
    size = guide_legend(override.aes = list(pch = 21, colour = "black", fill = "grey"))
  )


# Example 2 ---------------------------------------------------------------

# Here we'll restrict the data using the pathways/genes from the Sigora results,
# instead of using all pathways and layering our results on top
sigora_example_2

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

ggraph(candidates_as_network, layout = "mds") +
  facet_nodes(~direction, scales = "free", nrow = 2) +
  geom_edge_link(aes(edge_width = log10(similarity)), alpha = 0.5) +
  geom_node_point(
    aes(size = -log10(bonferroni), fill = grouped_pathway),
    pch = 21,
    colour = "black"
  ) +
  geom_node_label(
    aes(label = pathway_name_1),
    size = 5,
    repel = TRUE,
    alpha = 0.7,
    max.overlaps = 3,
    min.segment.length = 0
  ) +
  scale_edge_width(range = c(0.5, 2), name = "Similarity") +
  scale_size_continuous(
    labels = scales::label_math(expr = 10^-~.x),
    range = c(4, 8)
  ) +
  scale_fill_discrete(na.value = "grey") +
  labs(
    size = "Bonferroni\np-value",
    fill = "Pathway type"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.text.align = 0
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 6))
  )
