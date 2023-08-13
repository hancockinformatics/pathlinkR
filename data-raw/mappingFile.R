# Generate the mapping file used
# Use both biomaRt and AnnotDBI to get a map for Ensembl gene ID to Entrez ID
# and HGNC names

library(org.Hs.eg.db)
library(biomaRt)
library(AnnotationDbi)
library(tidyverse)

# Using biomaRt to create most of the mapping -------------------------------------------------

# Look for ENSG IDs from biomaRt and fill in any gaps
ensembl <- useEnsembl(biomart = 'genes')
ensembl <- useDataset(dataset = 'hsapiens_gene_ensembl', mart = ensembl)

# 64,945 ENSG IDs mapped from all 22 autosomes, X, Y, mitochondrial chromosomes
biomartChr <- getBM(filters = 'chromosome_name',
                    attributes = c('ensembl_gene_id',
                                   'entrezgene_id', 
                                   'hgnc_symbol',
                                   'description',
                                   'gene_biotype',
                                   'chromosome_name'),
                    values = c(1:22, 'X', 'Y', 'MT'),
                    mart = ensembl)

# 3,245 ENSG IDs are duplicated
biomartChr$ensembl_gene_id %>% duplicated() %>% sum()

# 35,998 ENSG IDs do not have an Entrez ID
is.na(biomartChr$entrezgene_id) %>% sum()

# Make the '' into NA for missing gene names for easier filtering
biomartChr[biomartChr == ''] <- NA

# Let's also remove any ENSG IDs that have neither a gene name nor an Entrez ID
# This leaves 45,868 ENSG IDs with either a gene name or an Entrez ID

biomartChr0 <- biomartChr %>%
    filter(!(is.na(entrezgene_id) & is.na(hgnc_symbol)))



# Fill in gaps in biomaRt using AnnotationDBI ----------------------------------

## Load AnnotationDBI mapping data ---------------------------------------------

# All the Ensembl gene IDs
ensemblKeys <- keys(org.Hs.eg.db, keytype = 'ENSEMBL')
annotdbiENSG <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = ensemblKeys, 
                                      columns = c('SYMBOL', 'ENTREZID'), 
                                      keytype = 'ENSEMBL') 
# Total 39,008 ENSG IDs mapping to gene name, Entrez IDs
# All ENSG IDs map to gene name and Entrez ID
anyNA(annotdbiENSG$SYMBOL)
anyNA(annotdbiENSG$ENTREZID)
# 3,478 Entrez IDs are duplicated
annotdbiENSG$ENTREZID[duplicated(annotdbiENSG$ENTREZID)] %>% length()

# All the Entrez IDs
entrezKeys <- keys(org.Hs.eg.db, keytype = 'ENTREZID')
annotdbiEntrez <- AnnotationDbi::select(org.Hs.eg.db, 
                                        keys = entrezKeys, 
                                        columns = c('SYMBOL', 'ENSEMBL'), 
                                        keytype = 'ENTREZID') 
# Total 69,580 Entrez ids mapping to gene name, ENSG IDs
# All Entrez IDs map to a gene name but only 39,008 map to ENSG ID
anyNA(annotdbiEntrez$SYMBOL)
anyNA(annotdbiEntrez$ENSEMBL)
sum(!is.na(annotdbiEntrez$ENSEMBL)) #39,008 are not NA

# Also, there are 372 duplicate ENSG IDs for different Entrez IDs
duplicatedENSG <- annotdbiEntrez$ENSEMBL[duplicated(annotdbiEntrez$ENSEMBL)]
duplicatedENSG[!is.na(duplicatedENSG)] %>% length() 

# Let's combine these two datasets and remove duplicates
annotdbiResults <- rbind(annotdbiENSG, annotdbiEntrez)
annotdbiResults <- annotdbiResults[!duplicated(annotdbiResults),]

## Fill in ENSG IDs with missing Entrez IDs ------------------------------------
biomartNoEntrez <- biomartChr0 %>% 
    filter(is.na(entrezgene_id)) # 15,921 without Entrez IDs
biomartNoEntrez %>% filter(is.na(hgnc_symbol))

# First, if any are duplicated, use the Entrez IDs from the duplicated entries
# For example, POLR2J3 has two ENSG IDs: ENSG00000168255 and ENSG00000285437
# Only one (ENSG00000285437) has an annotated Entrez ID, so apply that to the
# missing duplicate entry

# Gene names that are duplicated
dupGene <- biomartChr0$hgnc_symbol[duplicated(biomartChr0$hgnc_symbol)]
# Entries with a duplicated name without Entrez ID
biomartNoEntrez1 <- biomartNoEntrez %>% filter(hgnc_symbol %in% dupGene)

# Populate these with Entrez IDs from other entries
biomartAddEntrez <- biomartChr0 %>% 
    filter(hgnc_symbol %in% biomartNoEntrez1$hgnc_symbol,
           !is.na(entrezgene_id)) %>% 
    transmute(entrez_add = entrezgene_id, hgnc_symbol)

biomartAddEntrez <- biomartAddEntrez[!duplicated(biomartAddEntrez),]
biomartNoEntrez1 <- inner_join(biomartNoEntrez1, biomartAddEntrez) %>% 
    mutate(entrezgene_id = entrez_add) %>% 
    select(!entrez_add)

# add these back in
biomartChr1 <- biomartChr0 %>% 
    filter(!ensembl_gene_id %in% biomartNoEntrez1$ensembl_gene_id) %>% 
    rbind(biomartNoEntrez1)

# Second, how many are still missing Entrez IDs? Let's fill these in with 
# AnnotDBI based off of ENSG IDs
biomartNoEntrez2 <- biomartChr1 %>% 
    filter(is.na(entrezgene_id))
annotdbiEntrezAdd <- annotdbiEntrez %>% 
    filter(ENSEMBL %in% biomartNoEntrez2$ensembl_gene_id) %>% 
    transmute(ensembl_gene_id = ENSEMBL, ENTREZID)

# 8,619 can be added in
# Since they aren't missing names, we'll keep the names from biomaRt because
# they're probably more up to date
biomartAddEntrez2 <- inner_join(biomartNoEntrez2, annotdbiEntrezAdd) %>% 
    mutate(entrezgene_id = ENTREZID) %>% 
    select(!ENTREZID)

# Add these in
biomartChr2 <- biomartChr1 %>% 
    filter(!ensembl_gene_id %in% biomartAddEntrez2$ensembl_gene_id) %>% 
    rbind(biomartAddEntrez2)

# Third, let's add any remainder using AnnotDBI based on HGNC symbols
biomartNoEntrez3 <- biomartChr2 %>% 
    filter(is.na(entrezgene_id))

annotdbiEntrezAdd2 <- annotdbiResults %>% 
    filter(SYMBOL %in% biomartNoEntrez3$hgnc_symbol) %>% 
    transmute(hgnc_symbol = SYMBOL, ENTREZID) # 6,842 more

annotdbiEntrezAdd2 <- inner_join(annotdbiEntrezAdd2, biomartNoEntrez3) %>% 
    mutate(entrezgene_id = ENTREZID) %>% 
    select(!ENTREZID)

biomartChr3 <- biomartChr2 %>% 
    filter(!ensembl_gene_id %in% annotdbiEntrezAdd2$ensembl_gene_id) %>% 
    rbind(annotdbiEntrezAdd2)

# Fourth, are there any ENSG IDs with no name
biomartNoName <- biomartChr3 %>% 
    filter(is.na(hgnc_symbol)) # 2,026 without a gene name
# They all have an ENSG ID and Entrez ID though

# Do these have any annotations in AnnotDBI
annotdbiNamesEntrez <- annotdbiEntrez %>% 
    filter(ENTREZID %in% biomartNoName$entrezgene_id) # 1220 entries
annotdbiNamesENSG <- annotdbiEntrez %>% 
    filter(ENSEMBL %in% biomartNoName$ensembl_gene_id) # 746 entries

# Weirdly, these don't fully overlap. It looks like the AnnotDBI names are 
# sometimes outdated. For example, ENSG00000249624 maps to IFNAR2 in AnnotDBI
# but it is actually IFNAR2-IL10RB readthrough on biomaRt. So we are only
# keeping names that map to both the Entrez and ENSG IDs

annotdbiNames <- rbind(annotdbiNamesEntrez, annotdbiNamesENSG) %>% 
    mutate(ensembl_gene_id = ENSEMBL,
           entrezgene_id = ENTREZID) 
annotdbiNames <- annotdbiNames[!duplicated(annotdbiNames),] #1234 entries

biomartNoNameAdd <- inner_join(biomartNoName, annotdbiNames) %>% 
    mutate(hgnc_symbol = SYMBOL) %>% 
    select(!ENSEMBL) %>% 
    select(!ENTREZID) %>% 
    select(!SYMBOL) #732 entries

# add these in
biomartChr4 <- biomartChr3 %>% 
    filter(!ensembl_gene_id %in% biomartNoNameAdd$ensembl_gene_id) %>% 
    rbind(biomartNoNameAdd)

# Fifth, let's look at the ENSG IDs without Entrez or HGNC and see if we can
# add any new entries for those

biomartMissing <- biomartChr %>%
    filter(is.na(entrezgene_id) & is.na(hgnc_symbol))

annotdbiAllAdd <- annotdbiENSG %>% 
    filter(ENSEMBL %in% biomartMissing$ensembl_gene_id) %>% 
    mutate(ensembl_gene_id = ENSEMBL)

annotdbiAllAdd <- inner_join(biomartMissing, annotdbiAllAdd) %>% 
    mutate(hgnc_symbol = SYMBOL,
           entrezgene_id = ENTREZID) %>% 
    select(!ENSEMBL) %>% 
    select(!ENTREZID) %>% 
    select(!SYMBOL)

biomartChr5 <- biomartChr4 %>% 
    filter(!ensembl_gene_id %in% annotdbiAllAdd$ensembl_gene_id) %>% 
    rbind(annotdbiAllAdd)

# Clean up the mapping file ----------------------------------------------------
# Since the main point is to convert ENSG IDs into Entrez IDs or gene names, we
# need to make there are no duplicate ENSG IDs

biomartChr5 <- biomartChr5[!duplicated(biomartChr5),]
dupENSG <- biomartChr5$ensembl_gene_id[duplicated(biomartChr5$ensembl_gene_id)]
# 3251 ENSG IDs are duplicated
# What are these? Looks like a lot are pseudo and ribosome genes with multiple
# Entrez IDs to one ENSG ID
dupENSG <- biomartChr5 %>% filter(ensembl_gene_id %in% dupENSG)

# Decision is to keep the ENSG/HGNC/Entrez triplets in AnnotDBI as they
# appear to be the "original" version (usually the lower number); for example,
# ARHGAP29 has both 9411 and 105378858 as Entrez IDs. However, only 9411 is 
# documented online as its Entrez ID, the other is an LOC. On AnnotDBI 9411 is
# the only Entrez ID it is matched to.

cleanENSG <- inner_join(
    dupENSG,
    annotdbiENSG %>% mutate(entrezgene_id = ENTREZID,
                            ensembl_gene_id = ENSEMBL,
                            hgnc_symbol = SYMBOL))

# Deal with any duplicates manually by looking up ENSG ID on Ensembl
dupENSG2 <- cleanENSG$ensembl_gene_id[duplicated(cleanENSG$ensembl_gene_id)]
cleanENSG %>% filter(ensembl_gene_id %in% dupENSG2)
keepEntrez <- c('414243', '414062', '645166', '102723360', '101927745',
                 '105372935')
dupENSGdf <- cleanENSG %>% filter(entrezgene_id %in% keepEntrez)
cleanENSG <- cleanENSG %>% filter(!ensembl_gene_id %in% dupENSG2)
cleanENSG2 <- rbind(cleanENSG, dupENSGdf) %>% 
    select(!ENSEMBL) %>% 
    select(!SYMBOL) %>% 
    select(!ENTREZID)

# We lost some that didn't go through the first pass
# There's too many to manually go through. So just going to process the protein
# coding genes and ignore the pseudogenes/snRNAs/lncRNAs etc.
stillDupENSG <- dupENSG %>% 
    filter(!ensembl_gene_id %in% cleanENSG2$ensembl_gene_id,
           gene_biotype == 'protein_coding')

# Let's manually add them in again
keepEntrez2 <- c('645202', '646066', '391747', '112488736', '112488737',
                 '112488738', '122394733', '28227')
stillDupENSG <- stillDupENSG %>% filter(entrezgene_id %in% keepEntrez2)
cleanENSG3 <- rbind(cleanENSG2, stillDupENSG)

biomartChr6 <- biomartChr5 %>% 
    filter(!ensembl_gene_id %in% dupENSG$ensembl_gene_id) %>% 
    rbind(cleanENSG3)

# Summary ----------------------------------------------------------------------
# In total, there are 43,783 non-duplicated ENSG IDs
nrow(biomartChr6)
any(duplicated(biomartChr6$ensembl_gene_id))

# 359 ENSG IDs missing Entrez IDs and 401 non-NA duplicate Entrez IDs
is.na(biomartChr6$entrezgene_id) %>% sum() 
sum(duplicated(biomartChr6$entrezgene_id[!is.na(biomartChr6$entrezgene_id)]))

# 926 ENSG IDs missing gene symbol and 212 non-NA duplicate gene names
is.na(biomartChr6$hgnc_symbol) %>% sum() 
sum(duplicated(biomartChr6$hgnc_symbol[!is.na(biomartChr6$hgnc_symbol)]))

# Create the final mapping file
mappingFile <- biomartChr6 %>% 
    transmute(ensemblGeneId = ensembl_gene_id,
              hgncSymbol = hgnc_symbol,
              entrezGeneId = entrezgene_id) %>% 
    as.tibble()

usethis::use_data(mappingFile, overwrite=TRUE)

