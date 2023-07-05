# Load packages & data ----------------------------------------------------

library(sigora)
library(reactome.db)
library(tidyverse)

data("reaH", "idmap")


# Table of pathway sizes --------------------------------------------------

reactomePathwaySize <- c(
    reaH$L1$pwyszs, reaH$L2$pwyszs, reaH$L3$pwyszs, reaH$L4$pwyszs,
    reaH$L5$pwyszs
)

reactomePathwaySize <- tibble(
    pathwyId = names(reactomePathwaySize),
    num_genes_sigora = reactomePathwaySize
)

reactomePathwaySize <-
    reactomePathwaySize[!duplicated(reactomePathwaySize$pathwyId), ]


# Pathway IDs -------------------------------------------------------------

paths <- reaH$origRepo[[1]] %>%
    as.data.frame() %>%
    rownames_to_column("pwys") %>%
    as_tibble() %>%
    rename("pathwayId" = 2) %>%
    mutate(
        pwys = as.numeric(pwys),
        pathwayId = as.character(pathwayId)
    )


# Pathway genes -----------------------------------------------------------

genes <- reaH$origRepo[[2]] %>%
    as.data.frame() %>%
    rownames_to_column("gns") %>%
    as_tibble() %>%
    rename("entrezGeneId" = 2) %>%
    mutate(
        gns = as.numeric(gns),
        entrezGeneId = as.character(entrezGeneId)
    )


# Combine above results ---------------------------------------------------

map <- reaH$origRepo[[3]] %>% as_tibble()
mapNew <- left_join(map, genes)
mapNew <- left_join(mapNew, paths)

sigoraDatabaseFull <- left_join(
    mapNew[, 3:4],
    idmap,
    by = c("entrezGeneId" = "EntrezGene.ID"),
    multiple = "all"
) %>%
    rename(
        "ensemblGeneId" = "Ensembl.Gene.ID",
        "hgncSymbol" = "Symbol"
    )

# Find the level 1-4 pathways that are from Sigora
levelFourPathways <- reaH$L1$ps[!reaH$L1$ps %in% reaH$L5$ps] # 1001 pathways

sigoraDatabaseLevel4 <- sigoraDatabaseFull %>%
    filter(pathwayId %in% levelFourPathways)


# Add pathway descriptions ------------------------------------------------

pathwayNames <- reaH$pathwaydescriptions %>%
    as_tibble() %>%
    rename("pathwayId" = pwys, "pathwayName" = nms)

sigoraDatabaseLevel4 <- left_join(sigoraDatabaseLevel4, pathwayNames)


# Remove duplicate Ensembl IDs --------------------------------------------

# Decreased from 61045 to 60775 gene-pathway relations
sigoraDatabase <- sigoraDatabaseLevel4 %>%
    select(-entrezGeneId) %>%
    distinct()

usethis::use_data(sigoraDatabase, overwrite = TRUE)
