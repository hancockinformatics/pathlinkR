# Generate the mapping file used -----------------------------------------------
# Use both biomaRt and AnnotDBI to get a map for Ensembl gene ID to Entrez ID
# and HGNC names

library(org.Hs.eg.db)
library(biomaRt)
library(AnnotationDbi)
library(tidyverse)

# Creating a gene annotation map with AnnotDBI ---------------------------------

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


## biomaRt ---------------------------------------------------------------

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

# Let's also remove any ENSG IDs that have neither a gene name nor an Entrez ID
# This leaves 45,868 ENSG IDs with either a gene name or an Entrez ID
biomartChr[biomartChr == ''] <- NA
biomartChrClean <- biomartChr %>% 
    filter(!(is.na(entrezgene_id) & is.na(hgnc_symbol)))

# Fill in the gaps of missing Entrez IDs using AnnotDBI database of Entrez IDs
biomartNoEntrez <- biomartChrClean %>% 
    filter(is.na(entrezgene_id)) # 15,921 without Entrez IDs
biomartNoEntrez %>% filter(is.na(hgnc_symbol)) # all of them have names though

# First, if any are duplicated, use the Entrez IDs from the duplicated entries
# For example, POLR2J3 has two ENSG IDs: ENSG00000168255 and ENSG00000285437
# Only one (ENSG00000285437) has an annotated Entrez ID, so apply that to the
# missing duplicate entry

dupGene <- biomartChrClean$hgnc_symbol[duplicated(biomartChrClean$hgnc_symbol)]
biomartNoEntrez %>% filter(hgnc_symbol %in% dupGene)

# Can add in 8,730 Entrez IDs based from AnnotDBI
# Since they aren't missing names, we'll keep the names from biomaRt because
# they're probably more up to date
annotdbiEntrezAdd <- annotdbiEntrez %>% 
    filter(ENSEMBL %in% biomartNoEntrez$ensembl_gene_id) %>% 
    transmute(ensembl_gene_id = ENSEMBL, ENTREZID)

biomartAddEntrez <- inner_join(biomartNoEntrez, annotdbiEntrezAdd) %>% 
    mutate(entrezgene_id = ENTREZID) %>% 
    select(!ENTREZID)

# Add these in
biomartChr2 <- biomartChrClean %>% 
    filter(!ensembl_gene_id %in% biomartAddEntrez$ensembl_gene_id) %>% 
    rbind(biomartAddEntrez)

# Are there any with no name
biomartNoName <- biomartChr2 %>% 
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

biomartNoNameAdd <- inner_join(biomartNoName, annotdbiNames)


# 30,164 of 69,580 Entrez ids from AnnotDBI have annotatons
biomartEntrez <- getBM(filters = 'entrezgene_id',
                       attributes = c('ensembl_gene_id','entrezgene_id', 'hgnc_symbol','description','gene_biotype','chromosome_name'),
                       values = annotdbiResults$ENTREZID,
                       mart = ensembl)

# 43,126 of 39,008 ENSG IDs from AnnotDBI have annotations (duplicates)
biomartENSG <- getBM(filters = 'ensembl_gene_id',
                     attributes = c('ensembl_gene_id','entrezgene_id', 'hgnc_symbol','description','gene_biotype','chromosome_name'),
                     values = annotdbiResults$ENSEMBL,
                     mart = ensembl)

biomartResults <- rbind(biomartEntrez, biomartENSG)
biomartResults <- biomartResults[!duplicated(biomartResults),] #44,330 unique rows

# How many ensg_ids are duplicated?
biomartResults$ensembl_gene_id %>% duplicated() %>% sum() #4,631 duplicated

# What is wrong with these?
biomartResults %>% 
    filter(ensembl_gene_id %in% biomartResults$ensembl_gene_id[biomartResults$ensembl_gene_id %>% duplicated()]) %>% View()

# These multiple ensg ids, some have an odd chromosome location, not a number but a string e.g. CHR_HG1_PATCH
# Looking at a few manually, it looks like many of them are older ensg_ids
# Are any of these relevant, e.g. are they in the anntodbi_ensg_id keys or sequenced ensg_ids
chr_label <- biomartResults %>% filter(grepl('CHR', chromosome_name)) #2973 ensg_ids

any(chr_label$ensembl_gene_id %in% sequenced_ensg_ids) # none are in sequenced ensg_ids
chr_label$ensembl_gene_id %in% ensembl_keys %>% sum() # almost all (3350/3390) are in ensembl_keys from annotdbi. But they may be older versions, so remove for now

# removed those CHR entries, now only 208 duplicated entries
biomartResults_filtered <- biomartResults %>% filter(!grepl('CHR', chromosome_name))
biomartResults_filtered$ensembl_gene_id %>% duplicated() %>% sum() 

# What is wrong with the remainder duplicated ensg_ids? 385 rows with the same ensg_ids
duped_results <- biomartResults_filtered %>% filter(ensembl_gene_id %in% biomartResults_filtered$ensembl_gene_id[biomartResults_filtered$ensembl_gene_id %>% duplicated()])
# Many of these have the same gene symbol and ensg_id, but different entrez id.
# Which entrez_id is 'correct'?
# Can cross-reference to annotdbi results

duped_results$ensembl_gene_id %>% unique() %>% length() #177 unique ensg_ids with 385 entries

duped_combined <- left_join(duped_results, annotdbi_results %>% transmute(ensembl_gene_id = ENSEMBL, entrezgene_id = as.integer(ENTREZID), SYMBOL)) %>% 
    mutate(same = case_when(SYMBOL == hgnc_symbol ~ 'yes',
                            TRUE ~ 'no'))

# if the symbols from both annotdbi and biomart match, keep them
duped_consistent <- duped_combined %>% filter(same == 'yes') #167 ensg_ids
# are there any ensg_ids that are still duplicated
duped_consistent %>% filter(ensembl_gene_id %in% duped_consistent$ensembl_gene_id[duped_consistent$ensembl_gene_id %>% duplicated()])
# yes, CCL3L1/CCL3L3 and LINC00856/LINC00595
# Okay manually look at which ensg_ids they map to
# ENSG00000276085: CCL3L1 (based off ensembl.org)
# ENSG00000230417: LINC00595 (based off ensembl.org)
duped_consistent_edited <- duped_consistent %>% filter(!hgnc_symbol %in% c('CCL3L3', 'LINC00856')) #165 ensg_ids left

# Okay for the inconsistent ones, other than the ones that were replaced, are there still duplicates
duped_inconsistent <- duped_combined %>% filter(same != 'yes') #218 ensg_ids
duped_inconsistent %>% filter(ensembl_gene_id %in% duped_inconsistent$ensembl_gene_id[duped_inconsistent$ensembl_gene_id %>% duplicated()],
                              !ensembl_gene_id %in% duped_consistent$ensembl_gene_id) %>% write_csv('../duped_inconsistent.csv')

# No other way but to write this out as a csv file, manually curate, then return
write_csv(duped_inconsistent)
duped_inconsistent_edited <- read_csv('../duped_inconsistent.csv') #12 ensg_ids

duped_edited <- rbind(duped_consistent_edited, duped_inconsistent_edited) #177 ensg_ids
any(duped_edited$ensembl_gene_id %>% duplicated()) # ok no more duplicates

# Add these back (remove duplicated ensg_ids first)
biomart_final <- biomartResults_filtered %>% 
    filter(!ensembl_gene_id %in% c(duped_edited$ensembl_gene_id)) %>% 
    rbind(duped_edited %>% select(!c(SYMBOL, same)))

any(biomart_final$ensembl_gene_id %>% duplicated()) # no more duplicated ensg_ids, 60,619 ensg_ids

# Fill in empty entrez_ids or symbols ####
# There are a bunch that have no entrez_ids or hgnc symbols in this biomart_final dataframe
# Looks like annotdbi can be used to fill some in

### no entrez_id ####
no_entrez <- biomart_final %>% filter(is.na(entrezgene_id), hgnc_symbol != '') #15,223

# match by ensg_id and gene name to fill in any empty entrez_ids
no_entrez_joined1 <- left_join(no_entrez, annotdbi_results %>% transmute(ensembl_gene_id = ENSEMBL, ENTREZID, hgnc_symbol = SYMBOL)) %>% 
    mutate(entrezgene_id = ENTREZID) %>% 
    select(!ENTREZID) %>% 
    filter(!is.na(entrezgene_id)) #8256 that have entrez_ids filled in

# any duplicated? No
any(no_entrez_joined1$ensembl_gene_id %>% duplicated())

# For remainder, looks like some of the annotdbi data does not have ensg_id but has entrez_id and gene symbol, which is why it couldn't match by ensg_id
# Match based on solely symbol
no_entrez_joined2 <- left_join(
    no_entrez %>% filter(!ensembl_gene_id %in% no_entrez_joined1$ensembl_gene_id), 
    annotdbi_results %>% transmute(ENSEMBL, ENTREZID, hgnc_symbol = SYMBOL)) #6978 rows

# Some do not have any ensg_id, so doesn't really matter, can just fill in the entrez_id if it exists
# For those that do not have entrez_id, looks like they actually do online, but too much of a hassle to find. Most are lncRNA or pseudogenes anyways
no_entrez_joined3 <- no_entrez_joined2 %>% 
    filter(is.na(ENSEMBL)) %>% 
    mutate(entrezgene_id = ENTREZID) %>% 
    select(!c(ENSEMBL, ENTREZID)) #6907 that have entrez_ids filled in
any(no_entrez_joined3$ensembl_gene_id %>% duplicated()) # and no duplicates

# Okay, there are some that have conflicting ensg ids
no_entrez_joined4 <- no_entrez_joined2 %>% 
    filter(!is.na(ENSEMBL)) %>% 
    mutate(entrezgene_id = ENTREZID) #71 rows
# 6 of the conflicting ensg_ids are not already seen in the biomart annotation
no_entrez_joined4[!no_entrez_joined4$ENSEMBL %in% biomart_final$ensembl_gene_id,] 
# These 6 are seen in the weird CHR section
no_entrez_joined4[no_entrez_joined4$ENSEMBL %in% chr_label$ensembl_gene_id,]
# The conclusion is that:
# For conflicting ensg_ids - this is because multiple ensg_ids can map to a gene/entrez_id
# And those conflicting ensg_ids are already accounted for in the full biomart_final database
# The ones that are not are found in the CHR database, which we are not including anyways
# So adding the entrez_ids to these solely based on gene symbol should be okay
no_entrez_joined4 <-no_entrez_joined4 %>% 
    filter(!ENSEMBL %in% chr_label$ensembl_gene_id) %>% 
    select(!c(ENSEMBL, ENTREZID)) %>% #65 rows
    unique() #60 rows
any(no_entrez_joined4$ensembl_gene_id %>% duplicated()) # no duplicates

no_entrez_final <- rbind(no_entrez_joined1, no_entrez_joined3, no_entrez_joined4) #15223 rows that are annotated as best as one can

### no symbol ####
# There is an entrez ID but no symbol, give these genes a name
no_symbol <- biomart_final %>% filter(!is.na(entrezgene_id), hgnc_symbol == '') #855 

no_symbol_joined1 <- left_join(no_symbol, annotdbi_results %>% transmute(ensembl_gene_id = ENSEMBL, entrezgene_id = as.integer(ENTREZID), SYMBOL)) %>% 
    mutate(hgnc_symbol = SYMBOL) %>% 
    select(!SYMBOL) %>% 
    filter(!is.na(hgnc_symbol)) %>% 
    unique() #769 that have symbols filled in
any(no_symbol_joined1$ensembl_gene_id %>% duplicated()) # no duplicates

no_symbol_joined2 <- no_symbol %>% filter(!ensembl_gene_id %in% no_symbol_joined1$ensembl_gene_id) #86 that weren't annotated with above method
# probably because ensembl id did not match

no_symbol_joined3 <- no_symbol_joined2 %>% 
    left_join(annotdbi_results %>% transmute(ENSEMBL, entrezgene_id = as.integer(ENTREZID), SYMBOL)) %>% 
    mutate(hgnc_symbol = SYMBOL) %>% 
    select(!SYMBOL) #87
# again, are the ensg_ids already seen in either the biomart_final or CHR? Yes, so ignore ensg_id discrepancy
all(no_symbol_joined3$ensembl_gene_id %in% c(biomart_final$ensembl_gene_id, chr_label$ensembl_gene_id))

no_symbol_joined3 <- no_symbol_joined3 %>% 
    select(!ENSEMBL) %>% 
    unique() # 86, removed a duplicate
any(no_symbol_joined3$ensembl_gene_id %>% duplicated()) #no duplicates

no_symbol_final <- rbind(no_symbol_joined1, no_symbol_joined3) #855 genes annotated best as possible


### only ensg_id available ####
no_both <- biomart_final %>% filter(is.na(entrezgene_id), hgnc_symbol == '') #19,658

no_both_joined <- left_join(no_both, annotdbi_results %>% transmute(ensembl_gene_id = ENSEMBL, ENTREZID, SYMBOL)) %>% 
    mutate(hgnc_symbol = SYMBOL,
           entrezgene_id = as.integer(ENTREZID)) %>% 
    select(!c(ENTREZID, SYMBOL)) %>% 
    unique() #19,659 that have symbols and entrez_id filled in

no_both_joined[no_both_joined$ensembl_gene_id %>% duplicated(),] # one ensg_id duplicated
# ENSG00000249642, LOC107983963 vs LOC105377603 (the correct one)

no_both_final <- no_both_joined %>% filter(hgnc_symbol != 'LOC107983963' | is.na(hgnc_symbol))
any(no_both_final$ensembl_gene_id %>% duplicated()) #no duplicates, 19,658 entries

additional_annotated <- rbind(
    no_entrez_final,
    no_symbol_final,
    no_both_final
)

biomart_all_annotated <- rbind(
    biomart_final %>% filter(!ensembl_gene_id %in% additional_annotated$ensembl_gene_id),
    additional_annotated
) #60,619 ensembl gene ids annotated to the best of one's ability right now

# final check, are all sequenced ensg_ids in this dataframe? No...
all(sequenced_ensg_ids %in% biomart_all_annotated$ensembl_gene_id)
sequenced_ensg_ids[!sequenced_ensg_ids %in% biomart_all_annotated$ensembl_gene_id] #553 ensg_ids not in this
# quick looks appears that these ensg_ids are deprecated
# is this specific to a source (toronto vs quebec)
(!toronto_cts$ensembl_gene_id %in% biomart_all_annotated$ensembl_gene_id) %>% sum() #541 ensg_ids 
(!quebec_cts$ensg_id %in% biomart_all_annotated$ensembl_gene_id) %>% sum() #100 ensg_ids
# no, both datasets have deprecated ones
# decide to not include these in the final dataset, no point anyways since they don't annotate to anything

## write out file ####

# Reorganize columns
final_mapping <- biomart_all_annotated %>% transmute(
    ensg_id = ensembl_gene_id,
    gene_name = hgnc_symbol,
    entrez_id = entrezgene_id,
    gene_type = gene_biotype,
    desc_biomart = description,
    chromosome = chromosome_name) %>% arrange(ensg_id)

final_map_condense <- final_mapping %>% dplyr::select(ensg_id, gene_name, entrez_id)
mapping_readme <- 'This is a mapping guide for annotating ensembl ids to entrez ids and gene names.
Created Oct 12, 2022. R script titled 'gene_annotation.R'
It is a combination of using biomaRt and AnnotationDBI. 
Use the join() function to annotate DE gene lists, counts tables, etc.'
save(final_mapping, final_map_condense, mapping_readme, file = 'annotated_mapping_221012.rds')





# OLD CODE to make annotation file ####
# counts table of 58400 entries (5 of which are not annotated)
df <- read.csv('~/Desktop/COVID_data/counts_june2021.csv')
ensg_ids <- df$ensembl_gene_id %>% unlist()

# BioMart output is 58157 genes with duplicates
genes <- getBM(filters = 'hgnc_symbol',
               attributes = c('ensembl_gene_id','entrezgene_id', 'hgnc_symbol','description','gene_biotype'),
               values = 'ZNF285CP', 
               mart = ensembl)

# 453 ensembl ids not annotated by biomart
df %>% filter(!(ensembl_gene_id %in% genes$ensembl_gene_id)) %>% dplyr::select(ensembl_gene_id) -> not_in_biomart #453 ids

# AnnotationDBI output is 58639 genes with duplicates, all genes have annotations or are NA (some are duplicated)
annots <- AnnotationDbi::select(org.Hs.eg.db, # database
                                keys = ensg_ids,  # data to use for retrieval
                                columns = c('SYMBOL', 'ENTREZID','GENENAME'), # information to retreive for given data
                                keytype = 'ENSEMBL')

df %>% filter(!(ensembl_gene_id %in% annots$ENSEMBL)) %>% dplyr::select(ensembl_gene_id) -> not_in_annots #0 ids

# Create a dataframe 'total' with all genes in count table, biomart, and annotDBI annotations
combined <- data.frame(ensg_ids = ensg_ids[0:58395]) #remove the last five entries that are unmapped
names(genes)[1] <- 'ensg_ids'
names(annots)[1] <- 'ensg_ids'
total <- merge(combined, genes, by = 'ensg_ids', all.x = TRUE)
total <- merge(total, annots, by = 'ensg_ids', all.x = TRUE)

# Find the duplicated and unique ensg ids
ensg <- total$ensg_ids
unique_ensg <- total %>% filter(!ensg_ids %in% ensg[duplicated(ensg)]) #58173 ensg ids that only have one mapping for both annotDBI and biomart

dup_ensg <- total %>% filter(ensg_ids %in% ensg[duplicated(ensg)])
unique(dup_ensg$ensg_ids) %>% length() #222 ensg ids that are duplicated, 996 total entries

# Many of the duplicated entries are uncharacterized transcripts, pseudogenes, etc. filter those out
unchar <- dup_ensg %>% filter(grepl('uncharacterized', GENENAME))
pseudo <- dup_ensg %>% filter(grepl('pseudogene', GENENAME))
readthru <- dup_ensg %>% filter(grepl('readthrough', GENENAME))
linc <- dup_ensg %>% filter(grepl('LINC', SYMBOL))
loc <- dup_ensg %>% filter(grepl('LOC', SYMBOL))

cleaned_dup <- dup_ensg %>% filter(!(ENTREZID %in% c(unchar$ENTREZID, pseudo$ENTREZID, readthru$ENTREZID, linc$ENTREZID, loc$ENTREZID))) #Shrank to 597 entries with 180 genes (missing 42 ensg ids that aren't relevant anyways)

final_dup <- cleaned_dup %>% filter(entrezgene_id == ENTREZID & hgnc_symbol == SYMBOL) #Another filtering step looking for congruence between annotDBI and biomart, shrank to 167 genes with 166 ensg ids

# The 'lost' 56 ensg_ids (222 duplicated ids --> 166 ids after filtering
lost_ids <- dup_ensg %>% filter(!(ensg_ids %in% final_dup$ensg_ids)) %>% dplyr::select(ensg_ids) %>% unique() %>% unlist()
lost_genes <- filter(dup_ensg, ensg_ids %in% lost_ids) # 208 entries of 56 ensg ids
final_lost <- lost_genes %>% filter(entrezgene_id == ENTREZID & hgnc_symbol == SYMBOL) # 32 entries of 31 genes that are congruent, the last 25 genes are uncharacterized or pseudogenes

# The last 25 genes: 13 have consistent gene names but have incorrect entrez ids annotated by biomart, the rest are garbage genes anyways
remainders <- lost_genes %>% filter(!(ensg_ids %in% final_lost$ensg_ids)) %>% 
    filter(hgnc_symbol == SYMBOL) %>% transmute(ensg_ids, entrezgene_id = ENTREZID, hgnc_symbol, description, gene_biotype, SYMBOL, ENTREZID, GENENAME) 

# Combine all together
final_annot <- rbind(unique_ensg, final_dup, final_lost, remainders) #58385 entries of 58383 genes (out of 58395 genes from the start...good enough)

filter(final_annot, ensg_ids %in% final_annot$ensg_ids[duplicated(final_annot$ensg_ids)]) #two duplicates
#ENSG00000276085: CCL3L3 (414062) and CCL3L1 (6349) --> according to ENSEMBL, it is CCL3L1 (CCL3L3 is an alternative gene)
#ENSG00000230417: LINC00595 (414243) and LINC00856 (100132987) --> according to ENSEMBL, it is LINC00856
#Remove the duplicate ensg ids 
final_annot <- final_annot %>% filter(!(SYMBOL %in% c('CCL3L3','LINC00595')))
final_annot <- merge(combined, final_annot, by = 'ensg_ids', all.x = TRUE) # Add the last 12 genes as 'NA' at the end
final_annot[final_annot == ''] <- NA # Replace spaces with NA, weird processing kink from biomart

# Process each section of genes based on the amount of annotation they have

# Part 1: 16911 ids
# Have no annotations
section1 <- final_annot %>% filter(is.na(entrezgene_id) & is.na(hgnc_symbol) & is.na(ENTREZID) & is.na(SYMBOL)) %>% 
    mutate(final_symbol = NA,
           final_entrez = NA)

# Part 2: 6564 ids
# ensg IDs that have a gene symbol from biomart but nothing else
# choose only the protein coding ones and assign entrezids using annotdbi
# the remainder are non-coding etc and shouldn't really matter at this point, and will confuse as there are duplicate gene names (e.g. the same gene name and entrez id for protein coding and lnc entries)
section2 <- final_annot %>% filter(is.na(entrezgene_id) & !is.na(hgnc_symbol) & is.na(ENTREZID) & is.na(SYMBOL)) 

symbols <- section2 %>% filter(gene_biotype == 'protein_coding') %>% dplyr::select(hgnc_symbol) %>% unlist()
annots <- AnnotationDbi::select(org.Hs.eg.db, # database
                                keys = symbols,  # data to use for retrieval
                                columns = c('SYMBOL', 'ENTREZID','GENENAME'), # information to retreive for given data
                                keytype = 'SYMBOL')

names(annots)[1] <- 'hgnc_symbol'
section2 <- merge(section2, annots, by = 'hgnc_symbol', all.x = TRUE)
section2 <- section2 %>% transmute(ensg_ids,
                                   entrezgene_id,
                                   hgnc_symbol,
                                   description,
                                   gene_biotype,
                                   SYMBOL,
                                   ENTREZID = ENTREZID.x,
                                   GENENAME = GENENAME.x,
                                   final_symbol = hgnc_symbol,
                                   final_entrez = ENTREZID.y)

# Part 3: 42 ids
# Only has an entrez_id from biomart (with or without a gene symbol)
# Only one only had an entrez_id and not symbol: is a new gene (LOC389895)
section3 <- final_annot %>% filter(!is.na(entrezgene_id) & is.na(ENTREZID) & is.na(SYMBOL)) %>% 
    mutate(final_symbol = hgnc_symbol,
           final_entrez = entrezgene_id)

# Part 4: 1439 ids
# Only has an entrez id from annotdbi (with or without a gene symbol)
# None have only entrez id and no symbol
section4 <- final_annot %>% filter(is.na(entrezgene_id) & is.na(hgnc_symbol) & !is.na(ENTREZID)) %>% 
    mutate(final_symbol = SYMBOL,
           final_entrez = ENTREZID)

# Part 5: 8106 ids
# ensg ids that only have a gene symbol from biomart but both gene symbol and entrez id from annotdbi
# most have equivalent gene symbols, merge
# some have different gene symbols - some have LOC from annotdbi that haven't been updated, so keep the biomart gene symbol but the annotdbi entrez id
# for actually different ones, it seems that the annotdbi names are more 'updated' - so use annotdbi names (with biomart symbol in brackets) unless they are LOC or LINC gene names. Most of these are non-coding so again, not that important

section5 <- final_annot %>% filter(is.na(entrezgene_id) & !is.na(hgnc_symbol) & !is.na(ENTREZID) & !is.na(SYMBOL)) %>% 
    mutate(final_symbol = case_when(SYMBOL == hgnc_symbol ~ SYMBOL,
                                    SYMBOL != hgnc_symbol & (startsWith(SYMBOL, 'LOC') | startsWith(SYMBOL, 'LINC')) ~ hgnc_symbol,
                                    SYMBOL != hgnc_symbol & !(startsWith(SYMBOL, 'LOC') | startsWith(SYMBOL, 'LINC')) ~ paste0(SYMBOL, ' (', hgnc_symbol, ')')),
           final_entrez = ENTREZID)

# Part 6: 845 ids
# ensg ids that have an entrez id only from biomart but both gene symbol and entrez id from annotdbi
# use gene symbols and entrez ids from annotdbi (all the entrez ids are consistent)

section6 <- final_annot %>% filter(!is.na(entrezgene_id) & is.na(hgnc_symbol) & !is.na(ENTREZID) & !is.na(SYMBOL)) %>% 
    mutate(final_symbol = SYMBOL,
           final_entrez = ENTREZID)

# Part 7: 24488 ids
# ensg ids that are fully annotated
# all entrez ids are consistent, but some have different gene symbols. Use the annotdbi names (with biomart symbol in brackets) unless has LOC or LINC, then use biomart symbol

section7 <- final_annot %>% filter(!is.na(entrezgene_id) & !is.na(hgnc_symbol) & !is.na(ENTREZID) & !is.na(SYMBOL)) %>% 
    mutate(final_symbol = case_when(SYMBOL == hgnc_symbol ~ SYMBOL,
                                    SYMBOL != hgnc_symbol & (startsWith(SYMBOL, 'LOC') | startsWith(SYMBOL, 'LINC')) ~ hgnc_symbol,
                                    SYMBOL != hgnc_symbol & !(startsWith(SYMBOL, 'LOC') | startsWith(SYMBOL, 'LINC')) ~ paste0(SYMBOL, ' (', hgnc_symbol, ')')),
           final_entrez = ENTREZID)

# Combine all together and check no genes were lost:
final_mapping <- rbind(section1, section2, section3, section4, section5, section6, section7)
length(row.names(final_annot)) == length(row.names(final_mapping)) # Should be TRUE

# Reorganize columns
final_mapping <- final_mapping %>% transmute(
    ensg_id = ensg_ids,
    gene_name = final_symbol,
    entrez_id = final_entrez,
    gene_type = gene_biotype,
    name_biomart = hgnc_symbol,
    entrez_biomart = entrezgene_id,
    desc_biomart = description,
    name_annotdbi = SYMBOL,
    entrez_annotdbi = ENTREZID,
    desc_annotdbi = GENENAME) %>% arrange(ensg_id)

final_map_condense <- final_mapping %>% dplyr::select(ensg_id, gene_name, entrez_id)


# There are definitely a lot of duplicates in terms of gene names and entrez ids. How many of them are actually relevant, i.e. in the 'filtered' counts table that have low-counts ensg ids removed?

fcts <- df %>% filter((grepl('ENSG',ensembl_gene_id))) %>% #removes the non-ensg id entries like 'unmapped' etc.
    column_to_rownames(var = 'ensembl_gene_id') %>% 
    as.matrix()
keep <- rowMedians(fcts) >= 10
fcts <- fcts[keep,] 

filtered_df <- data.frame(ensg_id = row.names(fcts)) #14204 genes

filtered_map <- merge(filtered_df, final_map_condense)
with_names <- filtered_map %>% filter(!is.na(gene_name)) #13192 genes with annotations
with_entrez_id <- filtered_map %>% filter(!is.na(entrez_id)) #12816 genes with entrez ids

# Entries with duplicate entrez_ids
dup_entrez <- filter(filtered_map, entrez_id %in% filtered_map$entrez_id[duplicated(filtered_map$entrez_id)]) %>% 
    filter(!is.na(entrez_id))

# Entries with duplicate gene symbols
dup_name <- filter(filtered_map, gene_name %in% filtered_map$gene_name[duplicated(filtered_map$gene_name)]) %>%  
    filter(!is.na(gene_name))

# Check what these genes are in the mapping file
dups <- final_mapping %>% filter(ensg_id %in% c(dup_name$ensg_id, dup_entrez$ensg_id))
# The five genes are: c('BAZ2B', 'SLC35D2', 'POLR2J3', 'POLR2J4', 'TBCE')
# Some of the duplicates are for a protein coding and a non-coding entry
# Keep an eye out for these genes when doing DE analysis

mapping_readme <- 'This is a mapping guide for annotating ensembl ids to entrez ids. It is a combination of using biomaRt and AnnotationDBI. Use the merge() function to annotate DE gene lists, counts tables, etc.'

save(final_mapping, final_map_condense, mapping_readme, file = 'annotated_mapping.rds')