# Load packages & example data --------------------------------------------

devtools::load_all(".")
library(ggraph)
library(tidyverse)


# Example 1 ---------------------------------------------------------------

# |- get_pathway_distances and create_foundation --------------------------

dist_data <- get_pathway_distances(
  pathway_data = sigora_database,
  dist_method = "jaccard"
)

starting_pathways <- create_foundation(
  mat = dist_data,
  max_distance = 0.75
)


# |- create_pathnet -------------------------------------------------------

pathways_as_network <- create_pathnet(
  sigora_result = sigora_example_3,
  foundation = starting_pathways,
  trim = TRUE,
  trim_order = 2
) %>%
  left_join(
    top_pathways,
    by = c("pathway_1" = "pathway_id")
  ) #%>%
  #mutate(node_label = if_else(!is.na(direction), pathway_name_1, grouped_pathway))


# |- plot -----------------------------------------------------------------

ggraph(pathways_as_network, layout = "nicely") +
  geom_edge_link(aes(edge_width = log10(similarity)), alpha = 0.3) +
  geom_node_point(
    aes(size = -log10(bonferroni), fill = direction, colour = grouped_pathway),
    pch = 21,
  ) +
  geom_node_label(
    aes(label = description),
    # size = 6,
    repel = TRUE,
    alpha = 0.5,
    # max.overlaps = 3,
    min.segment.length = 0
  ) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_size_continuous(range = c(4, 8)) +
  scale_colour_brewer(palette = "Set3") +
  scale_fill_discrete(na.value = "white") +
  theme_void()


# Example 2 ---------------------------------------------------------------

# Here we'll restrict the data using the pathways/genes from the Sigora results,
# instead of using all pathways and layering our results on top
sigora_example_2

# Start by splitting the sigora genes, and adding the needed annotations from
# sigora_database
candidate_data <- sigora_example_2 %>%
  select(pathway_id, genes) %>%
  separate_rows(genes, sep = ";") %>%
  left_join(sigora_database, by = c("pathway_id", "genes" = "Symbol")) %>%
  relocate(pathway_id, Ensembl.Gene.ID, "Symbol" = genes, pathway_name) %>%
  distinct()

# Now that we have a smaller table in the same format as sigora_database, we can
# construct the foundation of pathway interactions
candidate_dist_data <- get_pathway_distances(
  pathway_data = candidate_data,
  dist_method = "jaccard"
)

candidate_starting_pathways <- create_foundation(
  mat = candidate_dist_data,
  max_distance = 0.75
)

candidates_as_network <- create_pathnet(
  sigora_result = sigora_example_2,
  foundation = candidate_starting_pathways,
  trim = FALSE
)

ggraph(candidates_as_network, layout = "nicely") +
  geom_edge_link(aes(edge_width = log10(similarity)), alpha = 0.2) +
  geom_node_point(
    aes(size = -log10(bonferroni), fill = direction),
    pch = 21,
    colour = "black"
  ) +
  geom_node_label(
    aes(label = node_label),
    size = 6,
    repel = TRUE,
    alpha = 0.75,
    max.overlaps = 3,
    min.segment.length = 0
  ) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_size_continuous(range = c(4, 8)) +
  scale_fill_discrete(na.value = "grey") +
  theme_void()
