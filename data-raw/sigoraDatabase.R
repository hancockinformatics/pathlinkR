# Load packages & data ----------------------------------------------------

library(sigora)
library(reactome.db)
library(tidyverse)

data("reaH", "idmap")


# Get some Sigora/Reactome data -------------------------------------------

# Pathway IDs
pathwayIds <- reaH$origRepo[[1]] %>%
    enframe("pwys", "pathwayId") %>%
    mutate(pathwayId = as.character(pathwayId))

# Pathway names
pathwayNames <- reaH$pathwaydescriptions %>%
    rename("pathwayId" = pwys, "pathwayName" = nms) %>%
    as_tibble()

# Pathway genes
pathwayGenes <- reaH$origRepo[[2]] %>%
    enframe("gns", "entrezGeneId") %>%
    mutate(entrezGeneId = as.character(entrezGeneId))


# Combine above results ---------------------------------------------------

mappingTable <- reaH$origRepo[[3]] %>%
    as_tibble() %>%
    left_join(pathwayGenes) %>%
    left_join(pathwayIds) %>%
    left_join(pathwayNames) %>%
    select(pathwayId, pathwayName, entrezGeneId)

sigoraDatabaseFull <- left_join(
    mappingTable,
    idmap,
    by = c("entrezGeneId" = "EntrezGene.ID"),
    multiple = "all"
) %>%
    rename("ensemblGeneId" = Ensembl.Gene.ID, "hgncSymbol" = Symbol)


# Find the level 1-4 pathways ---------------------------------------------

levelFourPathways <- reaH$L1$ps[!reaH$L1$ps %in% reaH$L5$ps]

sigoraDatabaseLevel4 <- sigoraDatabaseFull %>%
    filter(pathwayId %in% levelFourPathways)


# Remove duplicate Ensembl IDs --------------------------------------------

# Decreased from 61045 to 60775 gene-pathway relations
sigoraDatabase <- sigoraDatabaseLevel4 %>%
    select(-entrezGeneId) %>%
    distinct() %>%
    mutate(across(everything(), as.character))


# Save the data -----------------------------------------------------------

usethis::use_data(sigoraDatabase, overwrite = TRUE)
