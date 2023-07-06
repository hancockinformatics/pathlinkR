library(dplyr)

load("data/sigoraDatabase.rda")

identityTable <- sigoraDatabase %>%
  select(ensemblGeneId, pathwayId) %>%
  distinct() %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(
    id_cols     = pathwayId,
    names_from  = ensemblGeneId,
    values_from = present
  ) %>%
  replace(is.na(.), 0) %>%
  tibble::column_to_rownames("pathwayId") %>%
  as.matrix()

pathwayDistancesJaccard <- identityTable %>%
  vegan::vegdist(
    method = "jaccard",
    binary = TRUE,
    diag = TRUE
  ) %>%
  as.matrix()

usethis::use_data(pathwayDistancesJaccard, overwrite = TRUE, compress = "bzip2")
