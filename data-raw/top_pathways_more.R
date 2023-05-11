# Get top pathways for Reactome and Hallmark
# temporary script, remove later

library(msigdbr)
# library(ReactomePA)
library(sigora)
library(dplyr)
data("idmap", "reaH", package = "sigora")


# Reactome databases --------------------------------------------------------


# |- Download from Reactome -------------------------------------------------

# Relationship between each pathway and the one above it (hierarchy)
reactome_db <-
  read.csv(
    'https://reactome.org/download/current/ReactomePathwaysRelation.txt',
    sep = '\t',
    header = FALSE
  ) %>%
  filter(grepl('HSA', V1)) %>%
  as_tibble()
colnames(reactome_db) <- c('parent', 'child')

# Names of each pathway
reactome_names <-
  read.csv(
    'https://reactome.org/download/current/ReactomePathways.txt',
    sep = '\t',
    header = FALSE
  ) %>%
  filter(V3 == 'Homo sapiens') %>%
  as_tibble() %>%
  mutate(
    V2 = stringr::str_replace(
      V2,
      "Biosynthesis of electrophilic ω-3 PUFA oxo-derivatives",
      "Biosynthesis of electrophilic omega-3 PUFA oxo-derivatives"
    )
  )
colnames(reactome_names) <- c('pathway_id', 'pathway_name', 'species')

# Reactome provides some annotation on what pathways belong to which top pathways
reactome_top <- read.csv(
  'https://reactome.org/download/current/Complex_2_Pathway_human.txt',
  sep = '\t'
) %>%
  transmute(pathway_id = pathway, top_pathway = top_level_pathway) %>%
  unique() %>%
  as_tibble()


# |- Annotate all the pathways with top pathways, if available ------------

reactome_levels <- full_join(
  reactome_names,
  reactome_top,
  multiple = 'all'
) %>%
  left_join(
    reactome_names %>% transmute(top_pathway = pathway_id, top_pathway_name = pathway_name)
  )


# |- Fill in missing annotations for top pathways -------------------------

# Not all appear to have been annotated, refer to the hierarchy for the remainder
reactome_missing <- reactome_levels %>% filter(is.na(top_pathway))

# create a loop that keeps looping until it finds the top pathway for each
reactome_hierarchy <- list()

for (original_id in reactome_missing$pathway_id) {

  # Look in the hierarchy table first. Make sure the pathway in question is not
  # already a top pathway, i.e. it is not found in any of the pathways under
  # 'child'
  if (original_id %in% reactome_db$child) {
    order <- c()
    level <- original_id
    original_name <- filter(reactome_names, pathway_id == original_id)[2] %>% as.character

    # Loop until it hits the top pathway
    while (level != "character(0)") {
      pathway_id <- filter(reactome_db, child == level)[1] %>% as.character
      path_name <- filter(reactome_names, pathway_id == level)[2] %>% as.character
      order <- c(order, path_name)
      level <- pathway_id
    }

    # Reverse it so the top pathway is first
    order <- rev(order)
    top_path <- order[1]
    order <- paste(order, collapse="; ")

    # Add it to the data frame
    reactome_hierarchy <- reactome_hierarchy %>%
      append(list(c(original_id, original_name, top_path, order)))
  }
}

reactome_hierarchy_df <- reactome_hierarchy %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  as_tibble()

rownames(reactome_hierarchy_df) <- NULL
colnames(reactome_hierarchy_df) <-
  c('pathway_id', 'pathway_name', 'top_pathway_name', 'hierarchy')

reactome_hierarchy_df <- reactome_hierarchy_df %>% left_join(
  reactome_names %>% transmute(top_pathway = pathway_id, top_pathway_name = pathway_name)
)


# |- Deal with missing pathways and duplicates ----------------------------

# Some of these will have character(0) as top pathway because some belong to
# multiple top pathways. Some other pathways also were originally annotated to
# multiple pathways. Annotate these 39 manually to a single top pathway.
reactome_dupe <- plyr::rbind.fill(
  reactome_hierarchy_df %>% filter(top_pathway_name == 'character(0)'),
  reactome_levels %>% filter(pathway_name %in% reactome_levels$pathway_name[duplicated(reactome_levels$pathway_id)])
)

# Six pathways from SIGORA are outdated, but need to manually add them in for
# sigora to map properly
sigora_pathways <- as_tibble(reaH$pathwaydescriptions)
names(sigora_pathways) <- c("pathway_id", "pathway_name")

in_sigora <- sigora_pathways[!sigora_pathways$pathway_id %in% reactome_names$pathway_id,]

reactome_dupe <- plyr::rbind.fill(
  reactome_dupe,
  in_sigora
)

# |- Load the manually annotated duplicated pathways ----------------------

load("data/manual_dupe_annotation.rda")

# Add the missing ones back in
reactome_all_annotated <- plyr::rbind.fill(
  reactome_levels %>%
    filter(!is.na(top_pathway), !pathway_id %in% reactome_dupe$pathway_id),
  reactome_hierarchy_df %>%
    filter(top_pathway_name != 'character(0)') %>%
    select(!hierarchy),
  manual_dupe_annotation
)


# |- Make some checks -----------------------------------------------------

# Are all the child pathways from Reactome in the final annotated dataframe?
all(reactome_db$child %in% reactome_all_annotated$pathway_id)

# Are there any pathways that are not in the final annotated dataframe? These
# should only be the top pathways (29).
not_in_df <- reactome_names[!reactome_names$pathway_id %in% reactome_all_annotated$pathway_id, ]
all(not_in_df$pathway_name %in% reactome_all_annotated$top_pathway_name)

# Add these top pathways into the dataframe
not_in_df$top_pathway <- not_in_df$pathway_id
not_in_df$top_pathway_name <- not_in_df$pathway_name
reactome_all_annotated <- rbind(reactome_all_annotated, not_in_df)

# Check that there are no more duplicate pathway ids that belong to multiple top
# pathways?
any(!duplicated(reactome_all_annotated$pathway_id))

# Check that all the pathways used in Sigora are in this dataframe
all(sigora_pathways$pathway_id %in% reactome_all_annotated$pathway_id)

# Lastly, shrink top pathway names that are too long
reactome_all_annotated <- reactome_all_annotated %>% mutate(
  top_pathway_name_original = top_pathway_name,
  top_pathway_name = case_when(
    top_pathway_name == 'Gene expression (Transcription)' ~ 'Gene expression',
    top_pathway_name == 'Transport of small molecules' ~ 'Transport small molecules',
    top_pathway_name == 'Extracellular matrix organization' ~ 'ECM organization',
    top_pathway_name == 'Cellular responses to stimuli' ~ 'Cell responses to stimuli',
    top_pathway_name == 'Organelle biogenesis and maintenance' ~ 'Organelle biogenesis',
    TRUE ~ top_pathway_name_original
  )
)

# Add "grouped top pathways" for pathway network visualization
reactome_all_grouped <- reactome_all_annotated %>% mutate(
  grouped_pathway = case_when(
    top_pathway_name %in% c(
      'Autophagy',
      'ECM organization',
      'Organelle biogenesis',
      'Programmed Cell Death',
      'Protein localization',
      'Transport small molecules',
      'Vesicle-mediated transport'
    ) ~ 'Cell Process',
    top_pathway_name %in% c(
      'Cell Cycle',
      'Chromatin organization',
      'DNA Repair',
      'DNA Replication'
    ) ~ 'Cell Replication',
    top_pathway_name %in% c('Gene expression') ~ 'Gene Expression',
    top_pathway_name %in% c('Hemostasis', 'Immune System') ~ 'Immune/Hemostasis',
    top_pathway_name %in% c(
      'Metabolism',
      'Metabolism of proteins',
      'Metabolism of RNA',
      'Drug ADME'
    ) ~ 'Metabolism',
    top_pathway_name %in% c(
      'Cell responses to stimuli',
      'Cell-Cell communication',
      'Signal Transduction'
    ) ~ 'Signaling',
    top_pathway_name %in% c(
      'Circadian Clock',
      'Developmental Biology',
      'Digestion and absorption',
      'Muscle contraction',
      'Neuronal System',
      'Reproduction',
      'Sensory Perception'
    ) ~ 'Tissue Function',
    top_pathway_name %in% c('Disease') ~ 'Disease'
  )
)


# |- Deal with Disease terms ----------------------------------------------

# # The "Disease" top pathway varies, annotate them to the closest non-disease
# # pathway's top pathway
# reactome_disease <- reactome_all_grouped %>% filter(is.na(grouped_pathway))
#
# # Use the entrez mappers, the ensg mappers aren't available for all pathways
# reactome_genes_data <- read.csv(
#   'https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt',
#   sep = '\t',
#   header = FALSE
# )
#
# reactome_genes <- reactome_genes_data %>%
#   transmute(gene_id = V1, pathway_id = V2, pathway_name = V4) %>%
#   filter(grepl('HSA', pathway_id))
#
# # Are all the disease pathways in this database? Five that are not present, deal
# # with them manually afterwards
# reactome_disease$pathway_name[!reactome_disease$pathway_id %in% reactome_genes$pathway_id]
#
# # Combine all the genes from each "grouped pathway" together
# reactome_genes_grouped <- left_join(reactome_genes, reactome_all_grouped) %>%
#   mutate(pathway_id = case_when(
#     !is.na(grouped_pathway) ~ grouped_pathway,
#     TRUE ~ pathway_id
#   ))
#
# # Closest based on Jaccard index
# identity_table <- reactome_genes %>%
#   select(gene_id, pathway_id) %>%
#   distinct() %>%
#   mutate(present = 1) %>%
#   pivot_wider(
#     id_cols     = pathway_id,
#     names_from  = gene_id,
#     values_from = "present"
#   ) %>%
#   replace(is.na(.), 0) %>%
#   column_to_rownames(var = 'pathway_id') %>%
#   as.matrix()
#
# distance_matrix <- identity_table %>%
#   vegan::vegdist(
#     method = "jaccard",
#     binary = TRUE,
#     diag = TRUE
#   ) %>%
#   as.matrix()
#
# # Pathways that have a grouped pathway annotated to them
# annotated_pathways <- unique(
#   reactome_all_grouped %>% filter(!is.na(grouped_pathway)) %>% .$pathway_id
# )
#
#
# # Rows are all the pathways that aren't already annotated. Columns are all the
# # annotated pathways
# distance_matrix_filtered <- distance_matrix[
#   !row.names(distance_matrix) %in% annotated_pathways,
#   annotated_pathways[annotated_pathways %in% row.names(distance_matrix)]
# ]
#
# # Select the closest grouped pathway for each Disease pathway
# disease_grouped_pathway <- matrix(nrow = 0, ncol = 2)
# colnames(disease_grouped_pathway) <- c('pathway_id', 'grouped_pathway')
#
# for (i in 1:nrow(distance_matrix_filtered)) {
#   values <- distance_matrix_filtered[i, ]
#   pathway_id <- row.names(distance_matrix_filtered)[i]
#
#   closest_pathway <- data.frame(values) %>%
#     arrange(values) %>%
#     head(1) %>%
#     row.names()
#
#   grouped_pathway <- reactome_all_grouped %>%
#     filter(pathway_id == closest_pathway) %>%
#     .$grouped_pathway
#
#   disease_grouped_pathway <- rbind(
#     disease_grouped_pathway,
#     data.frame(pathway_id = pathway_id, grouped_pathway = grouped_pathway)
#   )
# }
#
# # Add it to the Disease pathways
# reactome_disease_grouped <- left_join(
#   reactome_disease %>% select(!grouped_pathway),
#   disease_grouped_pathway
# )
#
# # TODO Check it over and manually change some that don't make sense


# mSigDB Hallmark gene sets -----------------------------------------------

hallmark <- msigdbr(category = 'H')

hallmark <- hallmark %>%
  tidyr::separate(gs_name, sep = 'HALLMARK_', into = c('blank', 'gs_name'))

hallmark$gs_name <- gsub('_', ' ', hallmark$gs_name)

# Annotate the gene sets to their top gene set
cell_comp <- c("APICAL SURFACE", "APICAL JUNCTION", "PEROXISOME")
development <- c(
  "MYOGENESIS",
  "ANGIOGENESIS",
  "SPERMATOGENESIS",
  "ADIPOGENESIS",
  "EPITHELIAL MESENCHYMAL TRANSITION",
  "PANCREAS BETA CELLS"
)
dna_damage <- c("DNA REPAIR", "UV RESPONSE DN", "UV RESPONSE UP")
immune <- c(
  "ALLOGRAFT REJECTION",
  "COMPLEMENT",
  "COAGULATION",
  "INFLAMMATORY RESPONSE",
  "INTERFERON ALPHA RESPONSE",
  "INTERFERON GAMMA RESPONSE"
)
metabolism <- c(
  "BILE ACID METABOLISM",
  "CHOLESTEROL HOMEOSTASIS",
  "GLYCOLYSIS",
  "XENOBIOTIC METABOLISM",
  "FATTY ACID METABOLISM",
  "HEME METABOLISM",
  "OXIDATIVE PHOSPHORYLATION"
)
stress <- c(
  "HYPOXIA",
  "APOPTOSIS",
  "UNFOLDED PROTEIN RESPONSE",
  "PROTEIN SECRETION",
  "REACTIVE OXYGEN SPECIES PATHWAY"
)
proliferation <- c(
  "E2F TARGETS",
  "G2M CHECKPOINT",
  "MITOTIC SPINDLE",
  "P53 PATHWAY",
  "MYC TARGETS V1",
  "MYC TARGETS V2"
)
signaling <- c(
  'ANDROGEN RESPONSE',
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

hallmark_annotated <- hallmark %>% mutate(top_pathways = case_when(
  gs_name %in% cell_comp ~ 'Cellular',
  gs_name %in% development ~ 'Development',
  gs_name %in% dna_damage ~ 'DNA Damage',
  gs_name %in% immune ~ 'Immune',
  gs_name %in% metabolism ~ 'Metabolic',
  gs_name %in% stress ~ 'Stress',
  gs_name %in% proliferation ~ 'Growth',
  gs_name %in% signaling ~ 'Signaling'
))

hallmark_db <- hallmark_annotated %>%
  transmute(
    pathway_id = gs_name,
    top_pathways,
    pathway_name = gs_name,
    ensg_id = ensembl_gene,
    grouped_pathway = NA,
    top_pathways_original = NA
  )


# Save this as another top_pathways file ----------------------------------

# Difference between previous one is:
# 1. Used all Reactome pathways instead of just SIGORA pathways (2,621 pathways
#    over 1,298 originally)
# 2. Disease is now a grouped pathway instead of separating into other grouped
#    pathways
# 3. Added Hallmark gene sets
top_pathways_more <- rbind(
  reactome_all_grouped %>% transmute(
    pathway_id,
    top_pathways = top_pathway_name,
    pathway_name,
    grouped_pathway,
    top_pathways_original = top_pathway_name_original),
  hallmark_db %>% select(!ensg_id) %>% unique()
) %>%
  as_tibble()

usethis::use_data(top_pathways_more, overwrite = TRUE, compress = "bzip2")
