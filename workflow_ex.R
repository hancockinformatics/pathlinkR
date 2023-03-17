# Load packages & example data --------------------------------------------

devtools::load_all(".")
library(ggraph)
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

pathways_as_network <- create_pathnet(
  sigora_result = sigora_result_eg,
  foundation = starting_pathways,
  trim = TRUE,
  trim_order = 2
)


# |- Plot it! -------------------------------------------------------------

ggraph(pathways_as_network, layout = "nicely") +
  geom_edge_link(aes(edge_width = log10(similarity)), alpha = 0.3) +
  geom_node_point(
    aes(size = -log10(bonferroni), fill = direction),
    pch = 21,
    colour = "black"
  ) +
  geom_node_label(
    aes(label = node_label),
    size = 6,
    repel = TRUE,
    alpha = 0.5,
    max.overlaps = 3,
    min.segment.length = 0
  ) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_size_continuous(range = c(4, 8)) +
  scale_fill_discrete(na.value = "grey") +
  theme_void()
