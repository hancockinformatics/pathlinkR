# Load packages & data ----------------------------------------------------

library(sigora)
library(reactome.db)
library(tidyverse)

data("kegH", "idmap")


# Get some Sigora/Reactome data -------------------------------------------

# Pathway IDs
pathwayIds <- kegH$origRepo[[1]] %>%
    enframe("pwys", "pathwayId") %>%
    mutate(pathwayId=as.character(pathwayId))

# Pathway names
pathwayNames <- kegH$pathwaydescriptions %>%
    rename("pathwayId"=pwys, "pathwayName"=nms) %>%
    as_tibble()

# Pathway genes
pathwayGenes <- kegH$origRepo[[2]] %>%
    enframe("gns", "entrezGeneId") %>%
    mutate(entrezGeneId=as.character(entrezGeneId))


# Combine above results ---------------------------------------------------

mappingTable <- kegH$origRepo[[3]] %>%
    as_tibble() %>%
    left_join(pathwayGenes) %>%
    left_join(pathwayIds) %>%
    left_join(pathwayNames) %>%
    select(pathwayId, pathwayName, entrezGeneId)

keggDatabaseFull <- left_join(
    mappingTable,
    idmap,
    by=c("entrezGeneId" = "EntrezGene.ID"),
    multiple="all",
    relationship="many-to-many"
) %>%
    rename("ensemblGeneId"=Ensembl.Gene.ID, "hgncSymbol"=Symbol)


# Remove duplicate Ensembl IDs --------------------------------------------

# Decreased from 33658 to 32883 gene-pathway relations
keggDatabase <- keggDatabaseFull %>%
    select(-entrezGeneId) %>%
    mutate(across(everything(), as.character)) %>%
    distinct()


# Save the data -----------------------------------------------------------

usethis::use_data(keggDatabase, overwrite=TRUE)
