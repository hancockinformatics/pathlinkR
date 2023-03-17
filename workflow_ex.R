# Load packages & example data --------------------------------------------

devtools::load_all(".")
library(tidyverse)


# Approaches --------------------------------------------------------------

#' The way I see it, there are two options for how to construct these networks:
#' Option A
#' - Start with a distance matrix
#' - Define two pathways as sharing an edge based on a maximum distance cutoff
#' - Attach results from Sigora (direction, p-value, etc.)
#' - Convert this table to a network object
#'   - Remove subgraphs with no enriched nodes
#' - Plot it
#'
#' Option B
#' - Start with distance matrix
#' - Turn this into a network, including basically all pathway pairs as sharing
#'   an edge
#' - Convert this into a network object (very large!)
#' - Extract a subnetwork from this using the pathways from Sigora results as
#'   the target nodes. Attach Sigora information during this process
#' - Plot it


# Option A ----------------------------------------------------------------

# |- get_pathway_distances and create_foundation --------------------------

# Using the pre-calculated Jaccard distances, filter with a distance of 0.5,
# where any pair of pathways with distance < "max_distance" are considered
# connected.
dist_data <- get_pathway_distances(
  pathway_data = sigora_database,
  dist_method = "jaccard"
)

starting_pathways <- create_foundation(
  mat = dist_data,
  prop_to_keep = 0.002
)


# |- create_pathnet -------------------------------------------------------

pathways_as_network <-
  create_pathnet(sigora_result = sigora_result_eg, foundation = starting_pathways)


# |- Plot it! -------------------------------------------------------------

ggraph(pathways_as_network, layout = "kk") +
  geom_edge_link(aes(edge_width = similarity), alpha = 0.3) +
  geom_node_point(aes(size = -log10(bonferroni), fill = direction), pch = 21, colour = "black") +
  geom_node_label(aes(label = node_label), repel = TRUE, alpha = 0.75) +
  scale_edge_width(range = c(0.25, 1.5)) +
  scale_size_continuous(range = c(2, 6)) +
  theme_void()


# |- Remove empty subgraphs -----------------------------------------------

x1 <- pathways_as_network %>%
  filter(bonferroni < 1) %>%
  pull(rn)

valid_nodes <- map(x1, ~igraph::neighborhood(
  graph = as.igraph(pathways_as_network),
  order = 1,
  nodes = .x
)) %>%
  unlist() %>%
  unique()

pathways_as_network_trimmed <- pathways_as_network %>%
  filter(rn %in% c(x1, valid_nodes))

ggraph(pathways_as_network_trimmed, layout = "kk") +
  geom_edge_link(aes(edge_width = log10(1/distance)), alpha = 0.3) +
  geom_node_point(aes(size = -log10(bonferroni), fill = direction), pch = 21, colour = "black") +
  geom_node_label(aes(label = node_label), repel = TRUE, alpha = 0.75) +
  scale_edge_width(range = c(0.25, 1.5)) +
  scale_size_continuous(range = c(2, 6)) +
  theme_void()
