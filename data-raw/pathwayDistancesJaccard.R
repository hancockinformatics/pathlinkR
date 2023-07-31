# Load packages -----------------------------------------------------------

library(dplyr)


# Use existing sigoraDatabase to start ------------------------------------

load("data/sigoraDatabase.rda")


# Create identity matrix of all pathways ----------------------------------

identityTable <- sigoraDatabase %>%
  select(ensemblGeneId, pathwayId) %>%
  distinct() %>%
  mutate(present=1) %>%
  tidyr::pivot_wider(
    id_cols=pathwayId,
    names_from=ensemblGeneId,
    values_from=present
  ) %>%
  replace(is.na(.), 0) %>%
  tibble::column_to_rownames("pathwayId") %>%
  as.matrix()


# Calculate Jaccard distances ---------------------------------------------

pathwayDistancesJaccard <- identityTable %>%
  vegan::vegdist(
    method="jaccard",
    binary=TRUE,
    diag=TRUE
  ) %>%
  as.matrix()


# Save the data -----------------------------------------------------------

usethis::use_data(pathwayDistancesJaccard, overwrite=TRUE)
