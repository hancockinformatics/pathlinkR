# Load packages -----------------------------------------------------------

library(reactome.db)
library(org.Hs.eg.db)
library(tidyverse)


# Create a local Reactome database file ReactomePA ------------------------

humanGenes <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# Get mapping of pathway IDs to Entrez IDs
reactomeDb <- as.list(reactomePATHID2EXTID) %>%
    enframe() %>%
    unnest(cols = c(value)) %>%
    filter(grepl("HSA", name))

# Get mapping of pathway IDs to pathway names
reactomeNames <- as.list(reactomePATHID2NAME) %>%
    enframe() %>%
    unnest(cols = c(value)) %>%
    filter(grepl("HSA", name)) %>%
    mutate(pathwayName = str_remove(value, "Homo sapiens: ")) %>%
    dplyr::select(-value)

reactomeDb <- left_join(reactomeDb, reactomeNames) %>%
    dplyr::rename("pathwayId" = name, "entrezGeneId" = value)

# Filter out genes that aren't human; some are genes from other organisms (e.g.
# microbes for Immune System)
reactomeDatabase <- reactomeDb %>%
    filter(entrezGeneId %in% humanGenes)

usethis::use_data(reactomeDatabase, overwrite = TRUE)
