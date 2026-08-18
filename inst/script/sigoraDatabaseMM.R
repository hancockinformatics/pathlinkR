# Load packages & data ----------------------------------------------------

library(sigora)
library(reactome.db)
library(tidyverse)

data("reaM", "idmap")


# Get some Sigora/Reactome data -------------------------------------------

# Pathway IDs
pathwayIds <- reaM$origRepo[[1]] %>%
    enframe("pwys", "pathwayId") %>%
    mutate(pathwayId=as.character(pathwayId))

# Pathway names
pathwayNames <- reaM$pathwaydescriptions %>%
    rename("pathwayId"=pwys, "pathwayName"=nms) %>%
    as_tibble()

# Pathway genes
pathwayGenes <- reaM$origRepo[[2]] %>%
    enframe("gns", "entrezGeneId") %>%
    mutate(entrezGeneId=as.character(entrezGeneId))


# Combine above results ---------------------------------------------------

mappingTable <- reaM$origRepo[[3]] %>%
    as_tibble() %>%
    left_join(pathwayGenes) %>%
    left_join(pathwayIds) %>%
    left_join(pathwayNames) %>%
    select(pathwayId, pathwayName, entrezGeneId)

sigoraDatabaseFull <- left_join(
    mappingTable,
    idmap,
    by=c("entrezGeneId" = "EntrezGene.ID"),
    multiple="all",
    relationship="many-to-many"
) %>%
    rename("ensemblGeneId"=Ensembl.Gene.ID, "mgiSymbol"=Symbol)


# Find the level 1-4 pathways ---------------------------------------------

levelFourPathways <- reaM$L1$ps[!reaM$L1$ps %in% reaM$L5$ps]

sigoraDatabaseLevel4 <- sigoraDatabaseFull %>%
    filter(pathwayId %in% levelFourPathways)


# Remove duplicate Ensembl IDs --------------------------------------------

# Decreased from 61045 to 60775 gene-pathway relations
sigoraDatabaseMM <- sigoraDatabaseLevel4 %>%
    select(-entrezGeneId) %>%
    distinct() %>%
    drop_na(ensemblGeneId, mgiSymbol) |> 
    mutate(across(everything(), as.character))


# Save the data -----------------------------------------------------------

usethis::use_data(sigoraDatabaseMM, overwrite=TRUE)
