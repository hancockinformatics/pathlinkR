library(tidyverse)

identify_table <- sigora_database %>%
  select(Ensembl.Gene.ID, pathway_id) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols     = pathway_id,
    names_from  = Ensembl.Gene.ID,
    values_from = present
  ) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("pathway_id") %>%
  as.matrix()

pathway_distances_jaccard <- identify_table %>%
  vegan::vegdist(
    method = "jaccard",
    binary = TRUE,
    diag = TRUE
  ) %>%
  as.matrix()

usethis::use_data(pathway_distances_jaccard, overwrite = TRUE)
