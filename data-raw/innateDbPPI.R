# Load packages -----------------------------------------------------------

library(dplyr)


# Load and tidy InnateDB data ---------------------------------------------

innatedbInit <-
  readr::read_tsv("https://innatedb.com/download/interactions/all.mitab.gz")

# We just want human interactions
innatedbFiltered <- innatedbInit %>% filter(
  ncbi_taxid_A == "taxid:9606(Human)" & ncbi_taxid_B == "taxid:9606(Human)"
)

innatedbTrimmed <- innatedbFiltered %>%
  select(
    "ensemblGeneA"=alt_identifier_A,
    "ensemblGeneB"=alt_identifier_B
  ) %>%
  mutate(across(
    everything(),
    ~stringr::str_remove(.x, pattern="ensembl\\:")
  )) %>%
  distinct(ensemblGeneA, ensemblGeneB)


# Remove duplicate reversed interactions ----------------------------------

innatedbNoDups <- innatedbTrimmed[
  !duplicated(data.frame(t(apply(innatedbTrimmed, 1, sort)))),
]


# Remove promiscuous interactors ------------------------------------------

innateDbPPI <- innatedbNoDups %>%
  group_by(ensemblGeneA) %>%
  filter(n() < 1000) %>%
  ungroup() %>%
  group_by(ensemblGeneB) %>%
  filter(n() < 1000) %>%
  ungroup()


# Save the data -----------------------------------------------------------

usethis::use_data(innateDbPPI, overwrite=TRUE)
