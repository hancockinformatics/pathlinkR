# Load the required packages
library(biomaRt)
library(tidyverse)

# Use `biomaRt::getBM()` to create the conversion table, with the three most
# common human ID types.
biomart_id_mapping_human_0 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
)

# Tidy up the data
biomart_id_mapping_human <- biomart_id_mapping_human_0 %>%
  as_tibble() %>%
  rename("entrez_gene_id" = entrezgene_id) %>%
  mutate(
    across(everything(), as.character),
    across(everything(), ~str_replace(.x, "^$", NA_character_))
  ) %>%
  arrange(ensembl_gene_id, hgnc_symbol, entrez_gene_id) %>%
  distinct()

# Save data for use in the package
usethis::use_data(biomart_id_mapping_human, overwrite = TRUE)
