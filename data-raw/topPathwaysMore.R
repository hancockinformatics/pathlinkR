# Load packages -----------------------------------------------------------

library(msigdbr)
library(sigora)
library(tidyverse)
data("idmap", "reaH", package="sigora")


# Reactome databases --------------------------------------------------------

# |- Download from Reactome -------------------------------------------------

# Relationship between each pathway and the one above it (hierarchy)
reactomeDb <-
    read_tsv(
        "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
        col_names=c("parent", "child")
    ) %>%
    filter(grepl("HSA", parent))

# Names of each pathway
reactomeNames <-
    read_tsv(
        "https://reactome.org/download/current/ReactomePathways.txt",
        col_names=c("pathwayId", "pathwayName", "species")
    ) %>%
    filter(species == "Homo sapiens") %>%
    mutate(
        pathwayName=str_replace(
            pathwayName,
            "Biosynthesis of electrophilic ω-3 PUFA oxo-derivatives",
            "Biosynthesis of electrophilic omega-3 PUFA oxo-derivatives"
        )
    )

# Reactome provides some annotation on what pathways belong to which top pathways
reactomeTop <- read_tsv(
    "https://reactome.org/download/current/Complex_2_Pathway_human.txt"
) %>%
    select("pathwayId"=pathway, "topPathway"=top_level_pathway) %>%
    distinct()


# |- Annotate all the pathways with top pathways, if available ------------

reactomeLevels <- full_join(
    reactomeNames,
    reactomeTop,
    multiple="all"
) %>%
    left_join(
        select(
            reactomeNames,
            "topPathway"=pathwayId,
            "topPathwayName"=pathwayName
        )
    )


# |- Fill in missing annotations for top pathways -------------------------

# Not all appear to have been annotated, refer to the hierarchy for the
# remainder
reactomeMissing <- reactomeLevels %>% filter(is.na(topPathway))

# create a loop that keeps looping until it finds the top pathway for each
reactomeHierarchy <- list()

for (originalId in reactomeMissing$pathwayId) {

    # Look in the hierarchy table first. Make sure the pathway in question is not
    # already a top pathway, i.e. it is not found in any of the pathways under
    # "child"
    if (originalId %in% reactomeDb$child) {
        order <- c()
        level <- originalId
        originalName <-
            filter(reactomeNames, pathwayId == originalId)[2] %>% as.character

        # Loop until it hits the top pathway
        while (level != "character(0)") {
            pathwayId <- filter(reactomeDb, child == level)[1] %>% as.character
            pathName <-
                filter(reactomeNames, pathwayId == level)[2] %>% as.character
            order <- c(order, pathName)
            level <- pathwayId
        }

        # Reverse it so the top pathway is first
        order <- rev(order)
        topPath <- order[1]
        order <- paste(order, collapse="; ")

        # Add it to the data frame
        reactomeHierarchy <- reactomeHierarchy %>%
            append(list(c(
                originalId, originalName, topPath, order
            )))
    }
}

reactomeHierarchyDf <- reactomeHierarchy %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    as_tibble() %>%
    remove_rownames() %>%
    rename(
        "pathwayId"=V1,
        "pathwayName"=V2,
        "topPathwayName"=V3,
        "hierarchy"=V4
    )

reactomeHierarchyDf <- left_join(
    reactomeHierarchyDf,
    select(
        reactomeNames,
        "topPathway"=pathwayId,
        "topPathwayName"=pathwayName
    )
)


# |- Deal with missing pathways and duplicates ----------------------------

# Some of these will have character(0) as top pathway because some belong to
# multiple top pathways. Some other pathways also were originally annotated to
# multiple pathways. Annotate these 39 manually to a single top pathway.
reactomeDupe <- plyr::rbind.fill(
    filter(reactomeHierarchyDf, topPathwayName == "character(0)"),
    filter(
        reactomeLevels,
        pathwayName %in% reactomeLevels$pathwayName[duplicated(reactomeLevels$pathwayId)]
    )
)

# Six pathways from Sigora are outdated, but need to manually add them in for
# sigora to map properly
sigoraPathways <- as_tibble(reaH$pathwaydescriptions) %>%
    rename("pathwayId"=1, "pathwayName"=2)

inSigora <-
    sigoraPathways[!sigoraPathways$pathwayId %in% reactomeNames$pathwayId, ]

reactomeDupe <- plyr::rbind.fill(reactomeDupe, inSigora)


# |- Load the manually annotated duplicated pathways ----------------------

load("data/manualDupeAnnotation.rda")

# Add the missing ones back in
reactomeAllAnnotated <- plyr::rbind.fill(
    filter(
        reactomeLevels,
        !is.na(topPathway),
        !pathwayId %in% reactomeDupe$pathwayId
    ),
    select(filter(
        reactomeHierarchyDf,
        topPathwayName != "character(0)"
        ),
        -hierarchy
    ),
    manualDupeAnnotation
)


# |- Make some checks -----------------------------------------------------

# Are all the child pathways from Reactome in the final annotated dataframe?
all(reactomeDb$child %in% reactomeAllAnnotated$pathwayId)

# Are there any pathways that are not in the final annotated dataframe? These
# should only be the top pathways (29).
notInDf <-
    reactomeNames[!reactomeNames$pathwayId %in% reactomeAllAnnotated$pathwayId, ]
all(notInDf$pathwayName %in% reactomeAllAnnotated$topPathwayName)

# Add these top pathways into the dataframe
notInDf$topPathway <- notInDf$pathwayId
notInDf$topPathwayName <- notInDf$pathwayName
reactomeAllAnnotated <- rbind(reactomeAllAnnotated, notInDf)

# Check that there are no more duplicate pathway ids that belong to multiple top
# pathways?
any(!duplicated(reactomeAllAnnotated$pathwayId))

# Check that all the pathways used in Sigora are in this dataframe
all(sigoraPathways$pathwayId %in% reactomeAllAnnotated$pathwayId)

# Lastly, shrink top pathway names that are too long
reactomeAllAnnotated <- reactomeAllAnnotated %>%
    mutate(
        topPathwayNameOriginal=topPathwayName,
        topPathwayName=case_when(
            topPathwayName == "Gene expression (Transcription)" ~ "Gene expression",
            topPathwayName == "Transport of small molecules" ~ "Transport small molecules",
            topPathwayName == "Extracellular matrix organization" ~ "ECM organization",
            topPathwayName == "Cellular responses to stimuli" ~ "Cell responses to stimuli",
            topPathwayName == "Organelle biogenesis and maintenance" ~ "Organelle biogenesis",
            TRUE ~ topPathwayNameOriginal
        )
    )


# |- Add groupedTopPathways for pathway networks ---------------------------

reactomeAllGrouped <- reactomeAllAnnotated %>% mutate(
    groupedPathway=case_when(
        topPathwayName %in% c(
            "Autophagy",
            "ECM organization",
            "Organelle biogenesis",
            "Programmed Cell Death",
            "Protein localization",
            "Transport small molecules",
            "Vesicle-mediated transport"
        ) ~ "Cell Process",
        topPathwayName %in% c(
            "Cell Cycle",
            "Chromatin organization",
            "DNA Repair",
            "DNA Replication"
        ) ~ "Cell Replication",
        topPathwayName %in% c("Gene expression") ~ "Gene Expression",
        topPathwayName %in% c("Hemostasis", "Immune System") ~ "Immune/Hemostasis",
        topPathwayName %in% c(
            "Metabolism",
            "Metabolism of proteins",
            "Metabolism of RNA",
            "Drug ADME"
        ) ~ "Metabolism",
        topPathwayName %in% c(
            "Cell responses to stimuli",
            "Cell-Cell communication",
            "Signal Transduction"
        ) ~ "Signaling",
        topPathwayName %in% c(
            "Circadian Clock",
            "Developmental Biology",
            "Digestion and absorption",
            "Muscle contraction",
            "Neuronal System",
            "Reproduction",
            "Sensory Perception"
        ) ~ "Tissue Function",
        topPathwayName %in% c("Disease") ~ "Disease"
    )
)


# mSigDB Hallmark gene sets -----------------------------------------------

hallmark1 <- msigdbr(category="H")

hallmark2 <- hallmark1 %>%
    select(gs_name, ensembl_gene) %>%
    mutate(
        gs_name=str_remove(gs_name, "^HALLMARK_"),
        gs_name=str_replace_all(gs_name, "_", " ")
    )

# Annotate the gene sets to their top gene set
hallmarkCategories <- list(
    cellComp=c("APICAL SURFACE", "APICAL JUNCTION", "PEROXISOME"),
    development=c(
        "MYOGENESIS",
        "ANGIOGENESIS",
        "SPERMATOGENESIS",
        "ADIPOGENESIS",
        "EPITHELIAL MESENCHYMAL TRANSITION",
        "PANCREAS BETA CELLS"
    ),
    dnaDamage=c("DNA REPAIR", "UV RESPONSE DN", "UV RESPONSE UP"),
    immune=c(
        "ALLOGRAFT REJECTION",
        "COMPLEMENT",
        "COAGULATION",
        "INFLAMMATORY RESPONSE",
        "INTERFERON ALPHA RESPONSE",
        "INTERFERON GAMMA RESPONSE"
    ),
    metabolism=c(
        "BILE ACID METABOLISM",
        "CHOLESTEROL HOMEOSTASIS",
        "GLYCOLYSIS",
        "XENOBIOTIC METABOLISM",
        "FATTY ACID METABOLISM",
        "HEME METABOLISM",
        "OXIDATIVE PHOSPHORYLATION"
    ),
    stress=c(
        "HYPOXIA",
        "APOPTOSIS",
        "UNFOLDED PROTEIN RESPONSE",
        "PROTEIN SECRETION",
        "REACTIVE OXYGEN SPECIES PATHWAY"
    ),
    proliferation=c(
        "E2F TARGETS",
        "G2M CHECKPOINT",
        "MITOTIC SPINDLE",
        "P53 PATHWAY",
        "MYC TARGETS V1",
        "MYC TARGETS V2"
    ),
    signaling=c(
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
)

hallmarkAnnotated <- hallmark2 %>% mutate(
    topPathways=case_when(
        gs_name %in% hallmarkCategories$cellComp ~ "Cellular",
        gs_name %in% hallmarkCategories$development ~ "Development",
        gs_name %in% hallmarkCategories$dnaDamage ~ "DNA Damage",
        gs_name %in% hallmarkCategories$immune ~ "Immune",
        gs_name %in% hallmarkCategories$metabolism ~ "Metabolic",
        gs_name %in% hallmarkCategories$stress ~ "Stress",
        gs_name %in% hallmarkCategories$proliferation ~ "Growth",
        gs_name %in% hallmarkCategories$signaling ~ "Signaling"
    )
)

hallmarkDb <- hallmarkAnnotated %>%
    select(
        "pathwayId"=gs_name,
        topPathways,
        "pathwayName"=gs_name,
        "ensemblGeneId"=ensembl_gene,
    ) %>%
    mutate(
        groupedPathway=NA,
        topPathwaysOriginal=NA
    )


# Save this as another topPathways file ----------------------------------

# Difference between previous one is:
# 1. Used all Reactome pathways instead of just SIGORA pathways (2,621 pathways
#    over 1,298 originally)
# 2. Disease is now a grouped pathway instead of separating into other grouped
#    pathways
# 3. Added Hallmark gene sets
topPathwaysMore <- rbind(
    select(
        reactomeAllGrouped,
        pathwayId,
        "topPathways"=topPathwayName,
        pathwayName,
        groupedPathway,
        "topPathwaysOriginal"=topPathwayNameOriginal
    ),
    hallmarkDb %>% select(-ensemblGeneId) %>% distinct()
) %>%
    as_tibble()

usethis::use_data(topPathwaysMore, overwrite=TRUE, compress="bzip2")
