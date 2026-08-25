# Load packages -----------------------------------------------------------

library(reactome.db)
library(org.Mm.eg.db)
library(tidyverse)


# Create a local Reactome database file -----------------------------------

reactomeIds <- as.list(reactomePATHID2EXTID) %>%
    enframe("pathwayId", "entrezGeneId") %>%
    filter(str_detect(pathwayId, "^R-MMU")) %>%
    unnest(entrezGeneId)


reactomeNames <- as.list(reactomePATHID2NAME) %>%
    enframe("pathwayId", "pathwayName") %>%
    filter(str_detect(pathwayId, "^R-MMU")) %>%
    unnest(pathwayName) %>%
    mutate(pathwayName=str_remove(pathwayName, "Mus musculus: "))


# Map pathway IDs to pathway names ----------------------------------------

reactomeDb <- left_join(
    reactomeIds,
    reactomeNames,
    by="pathwayId",
    multiple="all"
)


# Filter non-mouse genes -------------------------------------------------

# Some are genes from other organisms (e.g. microbes for Immune System)
reactomeDatabaseMM <- reactomeDb %>%
    filter(entrezGeneId %in% keys(org.Mm.eg.db, keytype="ENTREZID"))


# Save the data -----------------------------------------------------------

usethis::use_data(reactomeDatabaseMM, overwrite=TRUE)
