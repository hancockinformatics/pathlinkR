# Load packages -----------------------------------------------------------

library(sigora)
library(tidyverse)


# Load data ---------------------------------------------------------------

data("reaH")
data("idmap")


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

genes <- left_join(genes, idmap)


# Combine above results ---------------------------------------------------

map <- reaH$origRepo[[3]] %>% as.data.frame()
map_new <- left_join(map, genes)
map_new <- left_join(map_new, paths)

sigora_database <- map_new[, 3:6]


# Add pathway descriptions ------------------------------------------------

pathway_names <- reaH$pathwaydescriptions

names(pathway_names) <- c("pathway_id", "pathway_name")

sigora_database <- left_join(sigora_database, pathway_names)


# Save the object for the package -----------------------------------------

usethis::use_data(sigora_database, overwrite = TRUE)
