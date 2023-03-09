# Load packages -----------------------------------------------------------

library(tidyverse)


# Helper functions --------------------------------------------------------

path_steps_local <- function(pathway_id) {

  # Looks for the top pathway
  if (filter(HSA_react, Child == pathway_id)[1] %>% as.character != "character(0)") {

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

    # Returns the top pathway name as well as the hierarchy of pathways
    # (arranged from highest to lowest)
    return(list(path_name, order))
  }

  if (pathway_id == "R-HSA-194840") {
    path_name <- "Signal Transduction"
    order <- paste0("Signal Transduction; Signaling by Rho GTPases, ",
                    "Miro GTPases and RHOBTB3; Signaling by Rho GTPases; ",
                    "RHO GTPase cycle")

    return(list(path_name, order))

  } else {
    # If already a top pathway, just return top pathway info
    path_name <- filter(reactome_names, hsa_id == pathway_id)[2] %>% as.character
    return(list(path_name, path_name))
  }
}

get_identity_matrix <- function(list) {
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


# First get Reactome relationships ----------------------------------------

reactome_db <- read_tsv(
  "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
  col_names = c("Parent", "Child")
)

reactome_names <- read_tsv(
  "https://reactome.org/download/current/ReactomePathways.txt",
  col_names = c("hsa_id", "name", "species")
)


# Add missing pathway -----------------------------------------------------

# R-HSA-194840 is missing (possibly updated), so we'll add it manually; it's
# "RHO GTPase cycle"
reactome_names <- rbind(
  reactome_names,
  tibble(hsa_id  = "R-HSA-194840",
         name    = "RHO GTPase cycle",
         species = "Homo sapiens")
)

HSA_react <- reactome_db %>% filter(grepl("HSA", Parent))


# Create the gene list for `get_identity_matrix()` usage --------------------------

# This script depends on the data generated from "data-raw/sigora_database.R"
load("data/sigora_database.rda")

all_pathways <- sigora_database$pathway_id %>% unique()

all_pathways_df <- data.frame(
  ID = NA,
  Top = NA,
  level = NA,
  hierarchy = NA
)

for (i in all_pathways) {
  results <- path_steps_local(i)
  hierarchy <- str_split(results[[2]], "; ") %>% unlist()
  df <- data.frame(
    ID = i,
    Top = results[[1]][1],
    level = length(hierarchy),
    hierarchy = results[[2]][1]
  )
  all_pathways_df <- rbind(all_pathways_df, df)
}

level_four_pathways <- all_pathways_df %>% filter(level <= 4, level != 1)
# 873 pathways total

pathway_description <- sigora_database %>%
  filter(pathway_id %in% level_four_pathways$ID)

jaccard_list <- list()

for (id in level_four_pathways$ID) {
  genes <- pathway_description %>% filter(pathway_id == id) %>%
    .$Ensembl.Gene.ID %>%
    as.character()
  list_genes <- list(id = genes)
  names(list_genes) <- id
  jaccard_list <- append(jaccard_list, list_genes)
}


# Matrix of distances between pathways ------------------------------------

# 0 = all overlap, 1 = no overlap
pathway_distances_jaccard <- jaccard_list %>%
  get_identity_matrix() %>%
  t() %>%
  vegan::vegdist(
    method = "jaccard",
    binary = TRUE,
    diag   = TRUE
  ) %>%
  as.matrix()

usethis::use_data(pathway_distances_jaccard, overwrite = TRUE)
