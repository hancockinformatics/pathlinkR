# Load packages & data ----------------------------------------------------

library(sigora)
library(tidyverse)

data("reaH", "idmap")


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


# Combine above results ---------------------------------------------------

map <- reaH$origRepo[[3]] %>% as.data.frame()
map_new <- left_join(map, genes)
map_new <- left_join(map_new, paths)

sigora_database_full <- left_join(map_new[, 3:4], idmap)

# Find the level 1-4 pathways that are from Sigora
level_four_pathways <- reaH$L1$ps[!reaH$L1$ps %in% reaH$L5$ps] # 1001 pathways

sigora_database_level4 <- sigora_database_full %>%
  filter(pathway_id %in% level_four_pathways)


# Add pathway descriptions ------------------------------------------------

pathway_names <- reaH$pathwaydescriptions
names(pathway_names) <- c("pathway_id", "pathway_name")

sigora_database_level4 <- left_join(sigora_database_level4, pathway_names)


# Remove duplicate Entrez IDs ---------------------------------------------

# Decreased from 61045 to 60775 gene-pathway relations
sigora_database <- sigora_database_level4 %>%
  select(-EntrezGene.ID) %>%
  distinct()


# Save the object for the package -----------------------------------------

usethis::use_data(sigora_database, overwrite = TRUE)
