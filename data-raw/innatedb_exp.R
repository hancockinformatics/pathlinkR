# Load packages
library(biomaRt)
library(dplyr)


# Load and filter InnateDB data
innatedb_init <-
  readr::read_tsv("https://innatedb.com/download/interactions/all.mitab.gz")


# We just want human interactions
innatedb_filtered <- innatedb_init %>% filter(
  ncbi_taxid_A == "taxid:9606(Human)" & ncbi_taxid_B == "taxid:9606(Human)"
)

innatedb_trimmed <- innatedb_filtered %>%
  select(
    "ensembl_gene_A" = alt_identifier_A,
    "ensembl_gene_B" = alt_identifier_B
  ) %>%
  mutate(across(
    everything(),
    ~stringr::str_remove(.x, pattern = "ensembl\\:")
  )) %>%
  distinct(ensembl_gene_A, ensembl_gene_B)


# Remove interactions that are the same, but reversed between the two columns
innatedb_no_dups <- innatedb_trimmed[
  !duplicated(data.frame(t(apply(innatedb_trimmed, 1, sort)))),
]


# Get gene mapping from biomaRt
biomart_mapping <- getBM(
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
  attributes = c("ensembl_gene_id", "hgnc_symbol")
)

innatedb_mapped <- innatedb_no_dups %>%
  left_join(
    biomart_mapping,
    by = c("ensembl_gene_A" = "ensembl_gene_id"),
    multiple = "all"
  ) %>%
  rename("hgnc_symbol_A" = hgnc_symbol) %>%
  left_join(
    biomart_mapping,
    by = c("ensembl_gene_B" = "ensembl_gene_id"),
    multiple = "all"
  ) %>%
  rename("hgnc_symbol_B" = hgnc_symbol) %>%
  relocate(ends_with("A"))


# Remove proteins/genes with more than 1000 interactions (e.g. UBC)
innatedb_exp <- innatedb_mapped %>%
  group_by(ensembl_gene_A) %>%
  filter(n() < 1000) %>%
  ungroup() %>%
  group_by(ensembl_gene_B) %>%
  filter(n() < 1000) %>%
  ungroup()


# Save the data, being sure to enable compression to reduce the size
usethis::use_data(innatedb_exp, overwrite = TRUE, compress = "bzip2")
