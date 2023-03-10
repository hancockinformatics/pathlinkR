# Load packages -----------------------------------------------------------

devtools::load_all(".")
library(ggraph)
library(tidygraph)
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

# Using the pre-calculated Jaccard distances, filter with a distance of 0.5,
# where any pair of pathways with distance < "max_distance" are considered connected.
starting_pathways <- create_foundation(
  mat = test_jaccard,
  max_distance = 25
) %>%
  glimpse()

starting_pathways_anno <- starting_pathways %>%
  left_join(
    distinct(select(sigora_database, pathway_id, pathway_name)),
    by = c("pathway_1" = "pathway_id")
  ) %>%
  left_join(
    distinct(select(sigora_database, pathway_id, pathway_name)),
    by = c("pathway_2" = "pathway_id"),
    suffix = c("_1", "_2")
  )

# |- Add sigora results ---------------------------------------------------

# Make sure we only have one of each pathway ID, otherwise theres problems when
# objects are turned into networks.
sigora_result_eg_slim <-
  read_tsv(file.path(
    "/mnt/analysis2/HIPC_2/Ontogeny/Ontogeny_analysis/Pathways",
    "sigora_paired_locationGAM_DOLDOL7.tsv"
  )) %>%
  select(pathway_id, direction, bonferroni) %>%
  arrange(bonferroni) %>%
  distinct(pathway_id, .keep_all = TRUE)


# |- Construct network ----------------------------------------------------

starting_nodes <- starting_pathways_anno %>%
  select(pathway_1, pathway_name_1) %>%
  distinct() %>%
  left_join(., sigora_result_eg_slim, by = c("pathway_1" = "pathway_id")) %>%
  arrange(bonferroni) %>%
  replace_na(list(bonferroni = 1))

starting_edges <- starting_pathways_anno %>%
  select(pathway_1, pathway_2, similarity, distance)

pathways_as_network <- tbl_graph(
  nodes = starting_nodes,
  edges = starting_edges,
  directed = FALSE
) %>%
  mutate(
    rn = row_number(),
    node_label = if_else(is.na(direction), "", pathway_name_1),
    node_label = map_chr(node_label, ~tRavis::tr_trunc_neatly(.x, l = 50)),
    node_label = str_replace(node_label, "^$", NA_character_)
  )


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
