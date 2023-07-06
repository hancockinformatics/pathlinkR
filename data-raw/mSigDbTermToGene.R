# Load packages -----------------------------------------------------------

library(msigdbr)
library(tidyverse)


# Top term information ----------------------------------------------------

topTerms <- list(
    cellComp = c("APICAL SURFACE", "APICAL JUNCTION", "PEROXISOME"),
    development = c(
        "MYOGENESIS",
        "ANGIOGENESIS",
        "SPERMATOGENESIS",
        "ADIPOGENESIS",
        "EPITHELIAL MESENCHYMAL TRANSITION",
        "PANCREAS BETA CELLS"
    ),
    dnaDamage = c("DNA REPAIR", "UV RESPONSE DN", "UV RESPONSE UP"),
    immune = c(
        "ALLOGRAFT REJECTION",
        "COMPLEMENT",
        "COAGULATION",
        "INFLAMMATORY RESPONSE",
        "INTERFERON ALPHA RESPONSE",
        "INTERFERON GAMMA RESPONSE"
    ),
    metabolism = c(
        "BILE ACID METABOLISM",
        "CHOLESTEROL HOMEOSTASIS",
        "GLYCOLYSIS",
        "XENOBIOTIC METABOLISM",
        "FATTY ACID METABOLISM",
        "HEME METABOLISM",
        "OXIDATIVE PHOSPHORYLATION"
    ),
    stress = c(
        "HYPOXIA",
        "APOPTOSIS",
        "UNFOLDED PROTEIN RESPONSE",
        "PROTEIN SECRETION",
        "REACTIVE OXYGEN SPECIES PATHWAY"
    ),
    proliferation = c(
        "E2F TARGETS",
        "G2M CHECKPOINT",
        "MITOTIC SPINDLE",
        "P53 PATHWAY",
        "MYC TARGETS V1",
        "MYC TARGETS V2"
    ),
    signaling = c(
        "ANDROGEN RESPONSE",
        "ESTROGEN RESPONSE EARLY",
        "ESTROGEN RESPONSE LATE",
        "KRAS SIGNALING DN",
        "KRAS SIGNALING UP",
        "MTORC1 SIGNALING",
        "NOTCH SIGNALING",
        "HEDGEHOG SIGNALING",
        "PI3K AKT MTOR SIGNALING",
        "TGF BETA SIGNALING",
        "IL6 JAK STAT3 SIGNALING",
        "IL2 STAT5 SIGNALING",
        "TNFA SIGNALING VIA NFKB",
        "WNT BETA CATENIN SIGNALING"
    )
) %>%
    enframe("topPathway", "pathwayId") %>%
    unnest(pathwayId)


# Get Hallmark info -------------------------------------------------------

hallmark <- msigdbr(category = "H") %>%
    select(gs_name, "ensemblGeneId" = ensembl_gene) %>%
    mutate(
        pathwayId = str_replace_all(gs_name, c("HALLMARK_" = "", "_" = " ")),
        pathwayName = pathwayId,
        groupedPathway = NA_character_,
        topPathwaysOriginal = NA_character_
    )


mSigDbTermToGene <- hallmark %>%
    left_join(topTerms, multiple = "all") %>%
    relocate(
        pathwayId,
        topPathway,
        pathwayName,
        ensemblGeneId,
        groupedPathway,
        topPathwaysOriginal
    ) %>%
    distinct(pathwayId, ensemblGeneId)


# Save the data -----------------------------------------------------------

usethis::use_data(mSigDbTermToGene, overwrite = TRUE)
