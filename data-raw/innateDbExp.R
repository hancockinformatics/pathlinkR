# Load packages -----------------------------------------------------------

library(biomaRt)
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
    "ensemblGeneA" = alt_identifier_A,
    "ensemblGeneB" = alt_identifier_B
  ) %>%
  mutate(across(
    everything(),
    ~stringr::str_remove(.x, pattern = "ensembl\\:")
  )) %>%
  distinct(ensemblGeneA, ensemblGeneB)


# Remove duplicate reversed interactions ----------------------------------

innatedbNoDups <- innatedbTrimmed[
  !duplicated(data.frame(t(apply(innatedbTrimmed, 1, sort)))),
]


# Get gene mapping from biomaRt -------------------------------------------

biomartMapping <- getBM(
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
  attributes = c("ensembl_gene_id", "hgnc_symbol")
)


# Map the interactors -----------------------------------------------------

innatedbMapped <- innatedbNoDups %>%
  left_join(
    biomartMapping,
    by = c("ensemblGeneA" = "ensembl_gene_id"),
    multiple = "all"
  ) %>%
  rename("hgncSymbolA" = hgnc_symbol) %>%
  left_join(
    biomartMapping,
    by = c("ensemblGeneB" = "ensembl_gene_id"),
    multiple = "all"
  ) %>%
  rename("hgncSymbolB" = hgnc_symbol) %>%
  relocate(ends_with("A"))


# Remove promiscuous interactors ------------------------------------------

innateDbExp <- innatedbMapped %>%
  group_by(ensemblGeneA) %>%
  filter(n() < 1000) %>%
  ungroup() %>%
  group_by(ensemblGeneB) %>%
  filter(n() < 1000) %>%
  ungroup()


# Save the data -----------------------------------------------------------

usethis::use_data(innateDbExp, overwrite = TRUE)
