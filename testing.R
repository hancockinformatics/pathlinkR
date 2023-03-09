# Load packages -----------------------------------------------------------

library(ggraph)
library(tidygraph)
library(ggVennDiagram)
library(tidyverse)
devtools::load_all(".")


# Example of overlaps -----------------------------------------------------

sister <- sigora_database %>% filter(pathway_name == "Mitotic Prometaphase")
kinesins <- sigora_database %>% filter(pathway_name == "Kinesins")
rho <- sigora_database %>% filter(pathway_name == "RHO GTPase Effectors")
macro <- sigora_database %>% filter(pathway_name == "Macroautophagy")

overlaps_pathways <- list(
  Mitotic = as.character(sister$Symbol),
  Kinesins = as.character(kinesins$Symbol),
  Rho = as.character(rho$Symbol)
)

ggVennDiagram(
  overlaps_pathways,
  color = "black",
  lwd = 0.4,
  label_size = 6,
  label = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "#F4FAFE", high = "purple") +
  theme(legend.position = "none") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  geom_text(size = 10)



# Old version -------------------------------------------------------------

# |- Create a data frame of connections with Jaccard index -------------------

rownames(filtered_jaccard) <- str_wrap(rownames(filtered_jaccard), width = 20)
colnames(filtered_jaccard) <- str_wrap(colnames(filtered_jaccard), width = 20)

network_df <- jaccard_mat_id %>%
  as.data.frame() %>%
  rownames_to_column(var = "pathway1") %>%
  pivot_longer(
    contains("R-HSA"),
    names_to = "pathway2",
    values_to = "jaccard"
  ) %>%
  mutate(jaccard = 1-jaccard)


# |- Annotate the pathway names ----------------------------------------------

network_df <- network_df %>%
  left_join(reactome_names %>% transmute(pathway1 = hsa_id, name1 = name))
network_df <- network_df %>%
  left_join(reactome_names %>% transmute(pathway2 = hsa_id, name2 = name))

# Pathways of interest
pathway <- "Nucleotide biosynthesis"
filtered_net_df <- network_df %>% filter(jaccard > 0.01, jaccard != 1) %>%
  select(name1, name2, jaccard) %>%
  filter(name1 == pathway|name2 == pathway)

graph <- as_tbl_graph(filtered_net_df)

ggraph(graph, layout = "kk") +
  geom_edge_link(aes(edge_width = jaccard)) +
  geom_node_point() +
  geom_node_text(aes(label = name), repel = TRUE, max.overlaps = Inf) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_edge_width(range = c(0.4, 1.5)) +
  theme_void() +
  scale_edge_color_manual(values = c("no" = "grey40","yes" = "red"))

ggsave("test_network.png", height = 20, width = 20)



# New version -------------------------------------------------------------

pathway_foundation <- create_foundation(pathway_distances_jaccard, max_distance = 0.5) %>%
  mutate(similarity = 1 - distance)

sigora_database_slim <- sigora_database %>%
  select(pathway_id, pathway_name) %>%
  distinct()

network_df_anno <- pathway_foundation %>%
  left_join(., sigora_database_slim, by = c("pathway_1" = "pathway_id")) %>%
  left_join(., sigora_database_slim, by = c("pathway_2" = "pathway_id"), suffix = c("_1", "_2")) %>%
  select(contains("name"), similarity)

pathway_network <- as_tbl_graph(network_df_anno)

ggraph(pathway_network, layout = "kk") +
  geom_edge_link(aes(edge_width = similarity), alpha = 0.3) +
  geom_node_point(pch = 21, size = 4, colour = "black", fill = "grey") +
  geom_node_text(aes(label = name), repel = TRUE, max.overlaps = 12) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_edge_width(range = c(0.3, 1)) +
  theme_void()
