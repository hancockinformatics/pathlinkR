# Load packages
library(msigdbr)
library(dplyr)

# Get Hallmark info
hallmark <- msigdbr(category = "H")

hallmark <- hallmark %>%
  tidyr::separate(gs_name, sep = "HALLMARK_", into = c("blank", "gs_name"))

hallmark$gs_name <- gsub("_", " ", hallmark$gs_name)

# Annotate the gene sets to their top gene set
top_terms <- list(
  cell_comp = c("APICAL SURFACE", "APICAL JUNCTION", "PEROXISOME"),
  development = c(
    "MYOGENESIS",
    "ANGIOGENESIS",
    "SPERMATOGENESIS",
    "ADIPOGENESIS",
    "EPITHELIAL MESENCHYMAL TRANSITION",
    "PANCREAS BETA CELLS"
  ),
  dna_damage = c("DNA REPAIR", "UV RESPONSE DN", "UV RESPONSE UP"),
  immune = c(
    "ALLOGRAFT REJECTION",
    "COMPLEMENT",
    "COAGULATION",
    "INFLAMMATORY RESPONSE",
    "INTERFERON ALPHA RESPONSE",
    "INTERFERON GAMMA RESPONSE"
  ),
  metabolism = c(
    "BILE ACID METABOLISM",
    "CHOLESTEROL HOMEOSTASIS",
    "GLYCOLYSIS",
    "XENOBIOTIC METABOLISM",
    "FATTY ACID METABOLISM",
    "HEME METABOLISM",
    "OXIDATIVE PHOSPHORYLATION"
  ),
  stress = c(
    "HYPOXIA",
    "APOPTOSIS",
    "UNFOLDED PROTEIN RESPONSE",
    "PROTEIN SECRETION",
    "REACTIVE OXYGEN SPECIES PATHWAY"
  ),
  proliferation = c(
    "E2F TARGETS",
    "G2M CHECKPOINT",
    "MITOTIC SPINDLE",
    "P53 PATHWAY",
    "MYC TARGETS V1",
    "MYC TARGETS V2"
  ),
  signaling = c(
    "ANDROGEN RESPONSE",
    "ESTROGEN RESPONSE EARLY",
    "ESTROGEN RESPONSE LATE",
    "KRAS SIGNALING DN",
    "KRAS SIGNALING UP",
    "MTORC1 SIGNALING",
    "NOTCH SIGNALING",
    "HEDGEHOG SIGNALING",
    "PI3K AKT MTOR SIGNALING",
    "TGF BETA SIGNALING",
    "IL6 JAK STAT3 SIGNALING",
    "IL2 STAT5 SIGNALING",
    "TNFA SIGNALING VIA NFKB",
    "WNT BETA CATENIN SIGNALING"
  )
)

hallmark_annotated <- hallmark %>% mutate(top_pathways = case_when(
  gs_name %in% top_terms$cell_comp ~ "Cellular",
  gs_name %in% top_terms$development ~ "Development",
  gs_name %in% top_terms$dna_damage ~ "DNA Damage",
  gs_name %in% top_terms$immune ~ "Immune",
  gs_name %in% top_terms$metabolism ~ "Metabolic",
  gs_name %in% top_terms$stress ~ "Stress",
  gs_name %in% top_terms$proliferation ~ "Growth",
  gs_name %in% top_terms$signaling ~ "Signaling"
))

hallmark_db <- hallmark_annotated %>%
  transmute(
    pathway_id = gs_name,
    top_pathways,
    pathway_name = gs_name,
    ensg_id = ensembl_gene,
    grouped_pathway = NA,
    top_pathways_original = NA
  )

msigdbr_t2g <- hallmark_db %>%
  distinct(pathway_id, ensg_id) %>%
  as_tibble()

# Save this for gene set enrichment for hallmark
usethis::use_data(msigdbr_t2g, overwrite = TRUE, compress = "bzip2")
