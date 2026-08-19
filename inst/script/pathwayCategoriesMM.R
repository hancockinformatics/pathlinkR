# Load packages -----------------------------------------------------------

library(sigora)
library(tidyverse)

data("idmap", "reaM", package="sigora")


# Reactome databases --------------------------------------------------------

# * Download from Reactome -------------------------------------------------

# Relationship between each pathway and the one above it (hierarchy)
reactomeDb <-
    read_tsv(
        "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
        col_names=c("parent", "child")
    ) %>%
    filter(grepl("MMU", parent))

# Names of each pathway
reactomeNames <-
    read_tsv(
        "https://reactome.org/download/current/ReactomePathways.txt",
        col_names=c("pathwayId", "pathwayName", "species")
    ) %>%
    filter(species == "Mus musculus")

# Reactome provides some annotation on what pathways belong to which top
# pathways. Had to email Reactome to get the mouse version; they have a 
# script which generated it:
# https://github.com/reactome/data-export/blob/main/src/main/java/org/reactome/server/export/tasks/Complex_2_Pathway_Human.java
reactomeTop <- read_tsv(
    "inst/extdata/Complex_2_Pathway_mouse.txt.gz"
) %>%
    select("pathwayId"=pathway, "topPathway"=top_level_pathway) %>%
    distinct()


# * Annotate all the pathways with top pathways, if available ------------

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


# * Fill in missing annotations for top pathways -------------------------

# Not all appear to have been annotated, refer to the hierarchy for the
# remainder
reactomeMissing <- reactomeLevels %>% filter(is.na(topPathway))

# Create a loop that keeps looping until it finds the top pathway for each
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
    ) %>%
    left_join(select(
        reactomeNames,
        "topPathway"=pathwayId,
        "topPathwayName"=pathwayName
    ))


# * Deal with missing pathways and duplicates ----------------------------

# Some of these will have character(0) as top pathway because some belong to
# multiple top pathways. Some other pathways also were originally annotated to
# multiple pathways. Annotate these 39 manually to a single top pathway.
reactomeDupe <- plyr::rbind.fill(
    filter(reactomeHierarchyDf, topPathwayName == "character(0)"),
    filter(
        reactomeLevels,
        pathwayName %in% reactomeLevels$pathwayName[duplicated(reactomeLevels$pathwayId)]
    )
) %>%
    as_tibble()

# Six pathways from Sigora are outdated, but need to manually add them in for
# Sigora to map properly
sigoraPathways <- as_tibble(reaM$pathwaydescriptions) %>%
    rename("pathwayId"=1, "pathwayName"=2) |> 
    mutate(across(everything(), as.character))

inSigora <-
    sigoraPathways[!sigoraPathways$pathwayId %in% reactomeNames$pathwayId, ]

reactomeDupe <- as_tibble(plyr::rbind.fill(reactomeDupe, inSigora))


# * Load the manually annotated duplicated pathways ----------------------

manualDupeAnnotation <- read_tsv("inst/extdata/manualDupeAnnotationMouse.tsv")

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
) %>% as_tibble()


# * Make some checks -----------------------------------------------------

# Are all the child pathways from Reactome in the final annotated dataframe?
all(reactomeDb$child %in% reactomeAllAnnotated$pathwayId)

# Are there any pathways that are not in the final annotated dataframe? These
# should only be the top pathways (28).
notInDf <-
    reactomeNames[!reactomeNames$pathwayId %in% reactomeAllAnnotated$pathwayId, ]

all(notInDf$pathwayName %in% reactomeAllAnnotated$topPathwayName)

# Add these top pathways into the dataframe
notInDf$topPathway <- notInDf$pathwayId
notInDf$topPathwayName <- notInDf$pathwayName
reactomeAllAnnotated <- rbind(reactomeAllAnnotated, notInDf)

# Check that there are no more duplicate pathway IDs that belong to multiple top
# pathways?
any(!duplicated(reactomeAllAnnotated$pathwayId))

# Check that all the pathways used in Sigora are in this dataframe
all(sigoraPathways$pathwayId %in% reactomeAllAnnotated$pathwayId)

# Lastly, shrink top pathway names that are too long
reactomeAllAnnotated <- reactomeAllAnnotated %>%
    mutate(
        topPathwayNameOriginal=topPathwayName,
        topPathwayName=case_when(
            topPathwayName == "Circadian Clock" ~ "Circadian clock",
            topPathwayName == "Gene expression (Transcription)" ~ "Gene expression",
            topPathwayName == "Transport of small molecules" ~ "Transport small molecules",
            topPathwayName == "Extracellular matrix organization" ~ "ECM organization",
            topPathwayName == "Cellular responses to stimuli" ~ "Cell responses to stimuli",
            topPathwayName == "Organelle biogenesis and maintenance" ~ "Organelle biogenesis",
            TRUE ~ topPathwayNameOriginal
        )
    )


# * Add groupedTopPathways for pathway networks ---------------------------

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
            "Circadian clock",
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

reactomeFinal <- reactomeAllGrouped %>%
    as_tibble() %>%
    select(
        pathwayId,
        pathwayName,
        "topLevelPathway" = topPathwayName,
        groupedPathway,
        "topLevelOriginal" = topPathwayNameOriginal
    )


# KEGG data ---------------------------------------------------------------

keggJSON <- jsonlite::fromJSON(paste0(
    "https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=",
    "json&filedir="
))

keggTidy <- keggJSON$children %>%
    as_tibble() %>%
    unnest(children, names_repair="universal") %>%
    unnest(children, names_repair="universal")

keggFinal <- keggTidy %>%
    mutate(
        pathwayId = paste0("mmu", str_extract(name...3, pattern="^[0-9]{5}")),
        pathwayName = str_trim(str_remove(name...3, pattern="^[0-9]{5}")),
        topLevelPathway = name...1,
        groupedPathway = name...1,
        topLevelOriginal = NA_character_
    ) %>%
    select(pathwayId:topLevelOriginal)


# Save this topPathways file ----------------------------------------------

pathwayCategoriesMM <- bind_rows(
    reactomeFinal,
    keggFinal
)

glimpse(pathwayCategoriesMM)

usethis::use_data(pathwayCategoriesMM, overwrite=TRUE)
