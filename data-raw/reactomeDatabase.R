# Load packages -----------------------------------------------------------

library(reactome.db)
library(org.Hs.eg.db)
library(tidyverse)


# Create a local Reactome database file -----------------------------------

reactomeIds <- as.list(reactomePATHID2EXTID) %>%
    enframe("pathwayId", "entrezGeneId") %>%
    filter(str_detect(pathwayId, "^R-HSA")) %>%
    unnest(entrezGeneId)


reactomeNames <- as.list(reactomePATHID2NAME) %>%
    enframe("pathwayId", "pathwayName") %>%
    filter(str_detect(pathwayId, "^R-HSA")) %>%
    unnest(pathwayName) %>%
    mutate(pathwayName=str_remove(pathwayName, "Homo sapiens: "))


# Map pathway IDs to pathway names ----------------------------------------

reactomeDb <- left_join(
    reactomeIds,
    reactomeNames,
    by="pathwayId",
    multiple="all"
)


# Filter non-human genes --------------------------------------------------

# Some are genes from other organisms (e.g. microbes for Immune System)
reactomeDatabase <- reactomeDb %>%
    filter(entrezGeneId %in% keys(org.Hs.eg.db, keytype="ENTREZID"))


# Save the data -----------------------------------------------------------

usethis::use_data(reactomeDatabase, overwrite=TRUE)
