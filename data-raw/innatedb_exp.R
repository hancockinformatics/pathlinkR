
# Load packages
library(dplyr)


# Load and filter InnateDB data
innatedb_init <- readr::read_tsv(
  file = "https://innatedb.com/download/interactions/all.mitab.gz",
  col_types = cols()
)

innatedb_filtered <- innatedb_init %>%
  filter(
    # We just want human interactions
    ncbi_taxid_A == "taxid:9606(Human)" & ncbi_taxid_B == "taxid:9606(Human)"
  )

innatedb_trimmed <- innatedb_filtered %>%
  select(
    "ensembl_gene_A" = alt_identifier_A,
    "ensembl_gene_B" = alt_identifier_B
  ) %>%
  mutate(across(everything(), stringr::str_remove, pattern = "ensembl\\:")) %>%
  distinct(ensembl_gene_A, ensembl_gene_B)


# Remove interactions that are the same, but reversed between the two columns
innatedb_no_dups <- innatedb_trimmed[
  !duplicated(data.frame(t(apply(innatedb_trimmed, 1, sort)))),
]


# Get gene mapping from biomaRt
biomart_mapping <- biomaRt::getBM(
  mart       = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  verbose    = TRUE
)

innatedb_mapped <- innatedb_no_dups %>%
  left_join(biomart_mapping, by = c("ensembl_gene_A" = "ensembl_gene_id")) %>%
  rename("hgnc_symbol_A" = hgnc_symbol) %>%
  left_join(biomart_mapping, by = c("ensembl_gene_B" = "ensembl_gene_id")) %>%
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


# Get a quick summary of the top interactors
innatedb_exp %>% count(hgnc_symbol_A) %>% arrange(desc(n))
innatedb_exp %>% count(hgnc_symbol_B) %>% arrange(desc(n))


usethis::use_data(innatedb_exp, overwrite = TRUE)
