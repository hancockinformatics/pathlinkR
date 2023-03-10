# Load packages -----------------------------------------------------------

library(sigora)
library(tidyverse)


# Load data from SIGORA  ---------------------------------------------------------------

data("reaH")
data("idmap")

# Helper functions --------------------------------------------------------

# gets pathway hierarchy
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
    order <- paste(order, collapse = "; ")

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


# Table of pathway sizes --------------------------------------------------

reactome_pathway_size <- c(reaH$L1$pwyszs, reaH$L2$pwyszs, reaH$L3$pwyszs,
                           reaH$L4$pwyszs, reaH$L5$pwyszs)

reactome_pathway_size <- tibble(
  pathwy.id        = names(reactome_pathway_size),
  num_genes_sigora = reactome_pathway_size
)

reactome_pathway_size <-
  reactome_pathway_size[!duplicated(reactome_pathway_size$pathwy.id), ]


# Pathway IDs -------------------------------------------------------------

paths <- reaH$origRepo[[1]] %>%
  as.data.frame() %>%
  rownames_to_column("pwys") %>%
  mutate(pwys = as.numeric(pwys)) %>%
  as_tibble()

names(paths)[2] <- "pathway_id"


# Pathway genes -----------------------------------------------------------

genes <- reaH$origRepo[[2]] %>%
  as.data.frame() %>%
  rownames_to_column("gns") %>%
  mutate(gns = as.numeric(gns)) %>%
  as_tibble()

names(genes)[2] <- "EntrezGene.ID"

#genes <- left_join(idmap, genes)


# Combine above results ---------------------------------------------------

map <- reaH$origRepo[[3]] %>% as.data.frame()
map_new <- left_join(map, genes)
map_new <- left_join(map_new, paths)

sigora_database_full <- left_join(map_new[, 3:4], idmap)

# find the level 1-4 pathways that are from sigora
level_four_pathways <- reaH$L1$ps[!reaH$L1$ps %in% reaH$L5$ps] #1001 pathways

sigora_database_level4 <- sigora_database_full %>% filter(pathway_id %in% level_four_pathways)

# Add pathway descriptions ------------------------------------------------

pathway_names <- reaH$pathwaydescriptions
names(pathway_names) <- c("pathway_id", "pathway_name")

sigora_database_level4 <- left_join(sigora_database_level4, pathway_names)

#just use ensembl ids, and remove duplicates from entrez id redundancy
#decreased from 61045 to 60775 gene-pathway relations
sigora_database <- sigora_database_level4 %>% select(-EntrezGene.ID) %>% distinct()

#
#
# all_pathways <- sigora_database$pathway_id %>% unique()
#
# all_pathways_df <- data.frame(
#   ID = NA,
#   Top = NA,
#   level = NA,
#   hierarchy = NA
# )
#
# cli_progress_bar(
#   "Getting pathway hierarchy information:",
#   total = length(all_pathways),
#   clear = FALSE
# )
# for (i in all_pathways) {
#   results <- path_steps_local(i)
#   hierarchy <- str_split(results[[2]], "; ") %>% unlist()
#   df <- data.frame(
#     ID = i,
#     Top = results[[1]][1],
#     level = length(hierarchy),
#     hierarchy = results[[2]][1]
#   )
#   all_pathways_df <- rbind(all_pathways_df, df)
#   cli_progress_update()
# }
# cli_progress_done()
#
# sigora_database <- sigora_database %>%
#   filter(pathway_id %in% level_four_pathways$ID)

# Save the object for the package -----------------------------------------

usethis::use_data(sigora_database, overwrite = TRUE)
