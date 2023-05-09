# Get top pathways for Reactome and Hallmark
# temporary script, remove later

library(msigdbr)
library(tidyverse)

# Reactome databases
# Download from Reactome

# Relationship between each pathway and the one above it (hierarchy)
reactome_db <- read.csv('https://reactome.org/download/current/ReactomePathwaysRelation.txt',
                        sep = '\t', header = FALSE) %>% filter(grepl('HSA', V1))
colnames(reactome_db) <- c('parent', 'child')

# Names of each pathway
reactome_names <- read.csv('https://reactome.org/download/current/ReactomePathways.txt',
                           sep = '\t', header = FALSE) %>% filter(V3 == 'Homo sapiens')
colnames(reactome_names) <- c('pathway_id', 'pathway_name', 'species')

reactome_top <- read.csv('https://reactome.org/download/current/Complex_2_Pathway_human.txt',
                         sep = '\t') %>%
  transmute(pathway_id = pathway,
            top_pathway = top_level_pathway) %>%
  unique()

reactome_levels <- full_join(reactome_names, reactome_top) %>%
  left_join(reactome_names %>% transmute(top_pathway = pathway_id,
                                         top_pathway_name = pathway_name))

# not all appear to have been annotated, refer to hierarchy for the remainder
reactome_missing <- reactome_levels %>% filter(is.na(top_pathway))
# create a loop that keeps looping until it finds the top pathway for each

reactome_hierarchy <- list()

for(original_id in reactome_missing$pathway_id){
  # look in the hierarchy table first
  # make sure the pathway in question is not already a top pathway
  # i.e. it is not found in any of the pathways under 'child'
  if (original_id %in% reactome_db$child) {
    order <- c()
    level <- original_id
    original_name <- filter(reactome_names, pathway_id == original_id)[2] %>% as.character
    # loop until it hits the top pathway
    while (level != "character(0)") {
      pathway_id <- filter(reactome_db, child == level)[1] %>% as.character
      path_name <- filter(reactome_names, pathway_id == level)[2] %>% as.character
      order <- c(order, path_name)
      level <- pathway_id
    }
    # reverse it so the top pathway is first
    order <- rev(order)
    top_path <- order[1]
    order <- paste(order, collapse="; ")

    # add it to the data frame
    reactome_hierarchy <- reactome_hierarchy %>%
      append(list(c(original_id, original_name, top_path, order)))
  }
}

reactome_hierarchy_df <- reactome_hierarchy %>% as.data.frame() %>% t() %>% as.data.frame()
rownames(reactome_hierarchy_df) <- NULL
colnames(reactome_hierarchy_df) <- c('pathway_id', 'pathway_name', 'top_pathway_name', 'hierarchy')
reactome_hierarchy_df <- reactome_hierarchy_df %>%
  left_join(reactome_names %>% transmute(top_pathway = pathway_id,
                                         top_pathway_name = pathway_name))

# some of these will have character(0) as top pathway because some belong to
# multiple top pathways. Annotate them manually

reactome_dupe <- reactome_hierarchy_df %>% filter(top_pathway_name == 'character(0)')
reactome_dupe
load("data/manual_dupe_annotation.rda")
reactome_dupe_annot <- reactome_dupe_annot %>%
  left_join(reactome_names) %>%
  left_join(reactome_names %>% transmute(top_pathway = pathway_id,
                                         top_pathway_name = pathway_name))

# add the missing ones back in

reactome_all_annotated <-
  plyr::rbind.fill(
    reactome_levels %>%
      filter(!is.na(top_pathway)),
    reactome_hierarchy_df %>%
      filter(top_pathway_name != 'character(0)') %>%
      select(!hierarchy),
    reactome_dupe_annot
  )

# to check
# are all the child pathways in the final annotated dataframe?
all(reactome_db$child %in% reactome_all_annotated$pathway_id)

# are there any pathways that are not in the final annotated dataframe?
# These should only be the top pathways
not_in_df <- reactome_names$pathway_name[!reactome_names$pathway_id %in% reactome_all_annotated$pathway_id]
all(not_in_df %in% reactome_all_annotated$top_pathway_name)

# Now, what are the pathways that belong to multiple top pathways?
# There are 34 pathways mapped to multiple top pathways
multiple_top <- reactome_all_annotated$pathway_name[duplicated(reactome_all_annotated$pathway_id)]
reactome_all_annotated %>% filter(pathway_name %in% multiple_top)
# when plotting, should join it to the top pathway with more enriched pathways

# lastly, shrink top pathway names that are too long
reactome_all_annotated <- reactome_all_annotated %>%
  mutate(
    top_pathway_name_original = top_pathway_name,
    top_pathway_name = case_when(
      top_pathway_name == 'Gene expression (Transcription' ~ 'Gene expression',
      top_pathway_name == 'Transport of small molecules' ~ 'Transport small molecules',
      top_pathway_name == 'Extracellular matrix organization' ~ 'ECM organization',
      top_pathway_name == 'Cellular responses to stimuli' ~ 'Cell responses to stimuli',
      top_pathway_name == 'Organelle biogenesis and maintenance' ~ 'Organelle biogenesis',
      TRUE ~ top_pathway_name_original
  ))

# For mSigDB Hallmark pathways

hallmark <- msigdbr(category = 'H')
hallmark <- hallmark %>% separate(gs_name, sep = 'HALLMARK_', into = c('blank', 'gs_name'))
hallmark$gs_name <- gsub('_', ' ', hallmark$gs_name)
