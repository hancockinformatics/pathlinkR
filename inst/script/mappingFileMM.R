# Load packages -----------------------------------------------------------

library(org.Mm.eg.db)
library(tidyverse)


# Using biomaRt to create most of the mapping -----------------------------

biomartData <- biomaRt::getBM(
  attributes=c(
    "ensembl_gene_id", "entrezgene_id", "mgi_symbol", "gene_biotype"
  ),
  filters="chromosome_name",
  values=c(1:22, "X", "Y", "MT"),
  mart=biomaRt::useDataset(
    dataset="mmusculus_gene_ensembl",
    mart=biomaRt::useEnsembl(biomart="genes")
  )
) %>%
  as_tibble() %>%
  mutate(
    across(everything(), as.character),
    across(everything(), ~na_if(.x, ""))
  )

## Let's remove any Ensembl IDs that have neither a gene name nor an Entrez ID
biomartDataNoNA <- biomartData %>%
  filter(!(is.na(entrezgene_id) & is.na(mgi_symbol)))

biomartColumns <- colnames(biomartDataNoNA)


# Use AnnotationDbi for some missing IDs ----------------------------------

annotationDbiEnsembl <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys=keys(org.Mm.eg.db, keytype="ENSEMBL"),
  columns=c("SYMBOL", "ENTREZID"),
  keytype="ENSEMBL"
) %>%
  as_tibble() %>%
  relocate(ENSEMBL, ENTREZID, SYMBOL)

annotationDbiEntrez <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys=keys(org.Mm.eg.db, keytype="ENTREZID"),
  columns=c("SYMBOL", "ENSEMBL"),
  keytype="ENTREZID"
) %>%
  as_tibble() %>%
  relocate(ENSEMBL, ENTREZID, SYMBOL)

annotationDbiResults <-
  distinct(bind_rows(annotationDbiEnsembl, annotationDbiEntrez))


# 1. Missing Entrez IDs with duplicated MGI symbols ----------------------

x1_biomartMissingEntrez <- biomartDataNoNA %>%
  group_by(mgi_symbol) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  filter(is.na(entrezgene_id)) %>%
  dplyr::select(!entrezgene_id)

x1_biomartAddedEntrez <- biomartDataNoNA %>%
  filter(
    mgi_symbol %in% x1_biomartMissingEntrez$mgi_symbol,
    !is.na(entrezgene_id)
  ) %>%
  distinct(entrezgene_id, mgi_symbol)

x1_fixedEntrez <- inner_join(
  x1_biomartMissingEntrez,
  x1_biomartAddedEntrez,
  by = "mgi_symbol"
) %>%
  relocate(all_of(biomartColumns))

mappingFile1 <- biomartDataNoNA %>%
  filter(!ensembl_gene_id %in% x1_fixedEntrez$ensembl_gene_id) %>%
  bind_rows(x1_fixedEntrez)


# 2. Remaining missing Entrez IDs -----------------------------------------

x2_biomartNoEntrez <- mappingFile1 %>%
  filter(is.na(entrezgene_id)) %>%
  dplyr::select(!entrezgene_id)

x2_annotationDbiEntrez <- annotationDbiEntrez %>%
  filter(ENSEMBL %in% x2_biomartNoEntrez$ensembl_gene_id) %>%
  dplyr::select("ensembl_gene_id" = ENSEMBL, "entrezgene_id" = ENTREZID)

x2_fixedEntrez <- inner_join(
  x2_biomartNoEntrez,
  x2_annotationDbiEntrez,
  by = "ensembl_gene_id"
) %>%
  relocate(all_of(biomartColumns))

mappingFile2 <- mappingFile1 %>%
  filter(!ensembl_gene_id %in% x2_fixedEntrez$ensembl_gene_id) %>%
  bind_rows(x2_fixedEntrez)


# 3. Add remainders using AnnotationDbi based on MGI symbols --------------

x3_biomartNoEntrez <- mappingFile2 %>%
  filter(is.na(entrezgene_id)) %>%
  dplyr::select(!entrezgene_id)

x3_annotationDbiEntrez <- annotationDbiResults %>%
  filter(SYMBOL %in% x3_biomartNoEntrez$mgi_symbol) %>%
  dplyr::select("mgi_symbol" = SYMBOL, "entrezgene_id" = ENTREZID)

x3_fixedEntrez <- inner_join(
  x3_annotationDbiEntrez,
  x3_biomartNoEntrez,
  by = "mgi_symbol"
) %>%
  relocate(all_of(biomartColumns)) |> 
  distinct()

mappingFile3 <- mappingFile2 %>%
  filter(!ensembl_gene_id %in% x3_fixedEntrez$ensembl_gene_id) %>%
  bind_rows(x3_fixedEntrez)


# 4. Ensembl IDs with no names --------------------------------------------

x4_biomartNoSymbol <- mappingFile3 %>%
  filter(is.na(mgi_symbol)) %>%
  dplyr::select(-mgi_symbol)

x4_annotationDbiSymbols <- annotationDbiResults %>%
  filter(
    ENTREZID %in% x4_biomartNoSymbol$entrezgene_id |
      ENSEMBL %in% x4_biomartNoSymbol$ensembl_gene_id
  ) %>%
  dplyr::rename(
    "ensembl_gene_id" = ENSEMBL,
    "entrezgene_id" = ENTREZID,
    "mgi_symbol" = SYMBOL
  )

x4_fixedSymbols <-
  inner_join(
    x4_biomartNoSymbol,
    x4_annotationDbiSymbols,
    by = c("ensembl_gene_id", "entrezgene_id")
  ) %>%
  relocate(all_of(biomartColumns))

mappingFile4 <- mappingFile3 %>%
  filter(!ensembl_gene_id %in% x4_fixedSymbols$ensembl_gene_id) %>%
  bind_rows(x4_fixedSymbols)


# 5. Ensembl IDs without Entrez IDs or MGI symbols -----------------------

x5_biomartMissingTwo <- biomartData %>%
  filter(is.na(entrezgene_id) & is.na(mgi_symbol))

x5_annotationDbiAddTwo <- annotationDbiEnsembl %>%
  filter(ENSEMBL %in% x5_biomartMissingTwo$ensembl_gene_id) %>%
  mutate(ensembl_gene_id = ENSEMBL)

x5_fixedTwo <- inner_join(
  x5_biomartMissingTwo,
  x5_annotationDbiAddTwo,
) %>%
  mutate(mgi_symbol = SYMBOL, entrezgene_id = ENTREZID) %>%
  dplyr::select(-c(ENSEMBL, ENTREZID, SYMBOL))

mappingFile5 <- mappingFile4 %>%
  filter(!ensembl_gene_id %in% x5_fixedTwo$ensembl_gene_id) %>%
  bind_rows(x5_fixedTwo)


# Clean up the mapping file -----------------------------------------------

duplicateEnsembl <- mappingFile5 %>%
  distinct() %>%
  group_by(ensembl_gene_id) %>%
  filter(n() > 1) %>%
  ungroup()

cleanEnsembl <- inner_join(
  duplicateEnsembl,
  mutate(
    annotationDbiEnsembl,
    entrezgene_id = ENTREZID,
    ensembl_gene_id = ENSEMBL,
    mgi_symbol = SYMBOL
  )
)

duplicateEnsembl2 <-
  cleanEnsembl$ensembl_gene_id[duplicated(cleanEnsembl$ensembl_gene_id)]

cleanEnsembl2 <- cleanEnsembl %>% filter(
  !ensembl_gene_id %in% duplicateEnsembl2
)

stillDuplicateEnsembl <- duplicateEnsembl %>% filter(
  !ensembl_gene_id %in% cleanEnsembl2$ensembl_gene_id,
  gene_biotype == "protein_coding"
)

cleanEnsembl3 <- bind_rows(cleanEnsembl2, stillDuplicateEnsembl)

mappingFileMM <- mappingFile5 %>%
  filter(!ensembl_gene_id %in% duplicateEnsembl$ensembl_gene_id) %>%
  bind_rows(cleanEnsembl3) %>%
  distinct() %>%
  select(
    "ensemblGeneId" = ensembl_gene_id,
    "mgiSymbol" = mgi_symbol,
    "entrezGeneId" = entrezgene_id
  ) %>%
  arrange(ensemblGeneId)


# Save the data -----------------------------------------------------------

usethis::use_data(mappingFileMM, overwrite=TRUE)
