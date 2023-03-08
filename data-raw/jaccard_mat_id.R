library(tidyverse)
library(ggVennDiagram)
library(ggraph)
library(tidygraph)

# Example of overlaps
# sister <- sigora_database %>% filter(pathway_name == 'Mitotic Prometaphase')
# kinesins <- sigora_database %>% filter(pathway_name == 'Kinesins')
# rho <- sigora_database %>% filter(pathway_name == 'RHO GTPase Effectors')
# macro <- sigora_database %>% filter(pathway_name == 'Macroautophagy')
# overlaps_pathways <- list(Mitotic = as.character(sister$Symbol), Kinesins = as.character(kinesins$Symbol), Rho = as.character(rho$Symbol))
#
#
# ggVennDiagram(overlaps_pathways, color = 'black', lwd = 0.4,
#               label_size = 6, label = "count", label_alpha = 0) +
#   scale_fill_gradient(low = "#F4FAFE", high = "purple") +
#   theme(legend.position = 'none') +
#   scale_x_continuous(expand = expansion(mult = .2)) +
#   geom_text(size = 10)

# first get reactome relationships

reactome_db <- read.table('https://reactome.org/download/current/ReactomePathwaysRelation.txt', sep = '\t') #https://reactome.org/download/current/ReactomePathwaysRelation.txt
reactome_names <- read.csv('https://reactome.org/download/current/ReactomePathways.txt', sep = '\t', header = FALSE)
names(reactome_names) <- c('hsa_id', 'name', 'species')

HSA_react <- reactome_db %>% filter(grepl('HSA', V1))
names(HSA_react) <- c('Parent', 'Child')

path_steps_local <- function(pathway_id) {
  if (filter(HSA_react, Child == pathway_id)[1] %>% as.character != "character(0)") { # Looks for the top pathway
    order <- c()
    level <- pathway_id
    while (level != "character(0)") {
      hsa_id <- filter(HSA_react, Child == level)[1] %>% as.character
      path_name <- filter(reactome_names, hsa_id == level)[2] %>% as.character
      order <- c(order, path_name)
      level <- hsa_id
    }
    order <- rev(order)
    order <- paste(order, collapse="; ")
    return(list(path_name, order)) #returns the top pathway name as well as the hierarchy of pathways (arranged from highest to lowest)
  }

  if(pathway_id == 'R-HSA-194840'){
    path_name <- 'Signal Transduction'
    order <- c('Signal Transduction; Signaling by Rho GTPases, Miro GTPases and RHOBTB3; Signaling by Rho GTPases; RHO GTPase cycle')
    return(list(path_name, order))

  } else {  # If already a top pathway, just return top pathway info
    path_name <- filter(reactome_names, hsa_id == pathway_id)[2] %>% as.character
    return(list(path_name, path_name))
  }
}

get_jac_mat <- function(list) {
  list %>%
    map(~data.frame(id = .x)) %>%
    bind_rows(.id = "name") %>%
    mutate(present = 1) %>%
    distinct(name, id, .keep_all = TRUE) %>%
    pivot_wider(
      id_cols     = "id",
      names_from  = "name",
      values_from = "present"
    ) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames(var = "id")
}

# create the gene list for get-jac-mat usage

all_pathways <- sigora_database$pathway_id %>% unique()

all_pathways_df <- data.frame(
  ID = NA,
  Top = NA,
  level = NA,
  hierarchy = NA
)

for(i in all_pathways){
  results <- path_steps_local(i)
  hierarchy <- str_split(results[[2]], '; ') %>% unlist()
  df <- data.frame(
    ID = i,
    Top = results[[1]][1],
    level = length(hierarchy),
    hierarchy = results[[2]][1]
  )
  all_pathways_df <- rbind(all_pathways_df, df)
}

level_four_pathways <- all_pathways_df %>% filter(level <= 4, level != 1) #873 pathways

pathway_description <- sigora_database %>% filter(pathway_id %in% level_four_pathways$ID)

jaccard_list <- list()

for(id in level_four_pathways$ID){
  genes <- pathway_description %>% filter(pathway_id == id) %>% .$Ensembl.Gene.ID %>% as.character()
  list_genes <- list(id = genes)
  names(list_genes) <- id
  jaccard_list <- append(jaccard_list, list_genes)
}

# matrix of distances between pathways, 0 = all overlap, 1 = no overlap
jaccard_mat_id <- jaccard_list %>%
  get_jac_mat() %>%
  t() %>%
  vegan::vegdist(
    method = "jaccard",
    binary = TRUE,
    diag   = TRUE
  ) %>%
  as.matrix()

usethis::use_data(jaccard_mat_id, overwrite = TRUE)

inverse_jaccard <- 1-jaccard_mat_id

filtered_jaccard <- inverse_jaccard

filtered_jaccard[filtered_jaccard<0.2] <- 0

rownames(filtered_jaccard) <- str_wrap(rownames(filtered_jaccard), width = 20)
colnames(filtered_jaccard) <- str_wrap(colnames(filtered_jaccard), width = 20)

network_df <- data.frame(path1 = 0, path2 = 0, jaccard = 0)
for(i in 1:ncol(filtered_jaccard)){
  column <- filtered_jaccard[c(i:ncol(filtered_jaccard)),i]
  df <- data.frame(path1 = rownames(filtered_jaccard)[i:ncol(filtered_jaccard)], path2 = colnames(filtered_jaccard)[i], jaccard = unlist(column))
  network_df <- rbind(network_df, df)
}

network_df <- network_df %>% filter(jaccard != 0, jaccard != 1)

graph <- as_tbl_graph(network_df)

ggraph(graph, layout = 'kk') +
  geom_edge_link(aes(edge_width = jaccard, alpha = jaccard*2)) +
  geom_node_point() +
  geom_node_text(aes(label = name), repel = TRUE, max.overlaps = Inf) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_edge_width(range = c(0.4, 1.5)) +
  theme_void() +
  scale_edge_color_manual(values = c('no' = 'grey40','yes' = 'red'))
