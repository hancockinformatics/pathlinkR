#' Human gene ID mappings from biomaRt
#'
#' A tibble containing gene ID mapping information for Ensembl, Entrez, and
#' HGNC gene identifiers.
#'
#' @format A data frame (tibble) with 75130 rows and 3 columns:
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene IDs}
#'   \item{hgnc_symbol}{HGNC symbols}
#'   \item{entrez_gene_id}{Entrez (NCBI) gene IDs}
#' }
"biomart_id_mapping_human"

#' List of example results from DESeq2
#'
#' @format A six-element list of data frames, each with 6 columns, containing
#'   example results from DESeq2. Each data frame has the following columns.
#' \describe{
#' \item{baseMean}{A combined score for the gene}
#' \item{log2FoldChange}{Fold change value for the gene}
#' \item{lfcSE}{Standard error for the fold change value}
#' \item{stat}{The statistic value}
#' \item{pvalue}{The nominal p value for the gene}
#' \item{padj}{The adjusted p value for the gene}
#' }
"deseq_example_list"

#' Experimentally verified human PPI data downloaded from InnateDB
#'
#' A data frame containing human PPI data from InnateDB, from the entry
#' "All Experimentally Validated Interactions (updated weekly)" at
#' <https://innatedb.com/redirect.do?go=downloadImported>. A few important steps
#' have been taken to filter the data, namely the removal of duplicate
#' interactions, and removing interactions that have the same components but are
#' swapped between A and B
#'
#' @format A data frame with 152259 rows and 4 columns:
#' \describe{
#'   \item{ensembl_gene_A}{Ensembl gene ID for the first gene/protein in the
#'     interaction}
#'   \item{hgnc_symbol_A}{HGNC symbol for the first gene/protein in the
#'     interaction}
#'   \item{ensembl_gene_B}{Ensembl gene ID for the second gene/protein in the
#'     interaction}
#'   \item{hgnc_symbol_B}{HGNC symbol for the second gene/protein in the
#'     interaction}
#' }
"innatedb_exp"

#' Table of human gene ID mappings
#'
#' @format A data frame with 42,099 rows and 3 columns
#' \describe{
#' \item{ensg_id}{Ensembl IDs}
#' \item{gene_name}{HGNC symbols}
#' \item{entrez_id}{NCBI Entrez IDs}
#' }
"mapping_file"

#' Matrix of pairwise Jaccard indices for all human Reactome pathways
#'
#' @format A matrix with 1001 rows and columns
#' \describe{
#' \item{Rows}{Reactome pathway identifiers}
#' \item{Columns}{Reactome pathway identifiers}
#' }
"pathway_distances_jaccard"

#' Table of all Sigora pathways and their constituent genes
#'
#' @format A tibble with 60775 rows and 4 columns
#' \describe{
#' \item{pathway_id}{Reactome pathway identifier}
#' \item{EnsemblGene.ID}{Ensembl gene identifier}
#' \item{Symbol}{HGNC gene symbol}
#' \item{pathway_name}{Reactome pathway description}
#' }
"sigora_database"

#' Data frame of Sigora results to be used for testing, as an example
#'
#' @format A data frame (tibble) with 44 rows and 7 columns
#' \describe{
#' \item{direction}{Type of dysregulation, based on fold change values from input DE genes}
#' \item{pathway_id}{Reactome pathway identifier}
#' \item{description}{Pathway name}
#' \item{pvalue}{Nominal p-value for enrichment}
#' \item{bonferroni}{Adjusted p-value for enrichment}
#' \item{level_1}{Top level Reactome term for the pathway}
#' \item{level_2}{Second level Reactome term for the pathway}
#' }
"sigora_example_1"

#' Another data frame of Sigora results to be used for testing, as an example
#'
#' @format A data frame (tibble) with 23 rows and 8 columns
#' \describe{
#' \item{pathway_id}{Reactome pathway identifier}
#' \item{description}{Pathway name}
#' \item{direction}{Type of dysregulation, based on fold change values from input DE genes}
#' \item{pvalue}{Nominal p-value for enrichment}
#' \item{bonferroni}{Adjusted p-value for enrichment}
#' \item{level_1}{Top level Reactome term for the pathway}
#' \item{level_2}{Second level Reactome term for the pathway}
#' \item{genes}{Overlapping genes between the pathway and input list}
#' }
"sigora_example_2"

#' Yet another data frame of Sigora results to be used for testing, as an example
#'
#' @format A data frame (tibble) with 23 rows and 8 columns
#' \describe{
#' \item{direction}{Type of dysregulation, based on fold change values from input DE genes}
#' \item{pathway_id}{Reactome pathway identifier}
#' \item{description}{Pathway name}
#' \item{pvalue}{Nominal p-value for enrichment}
#' \item{bonferroni}{Adjusted p-value for enrichment}
#' \item{level_1}{Top level Reactome term for the pathway}
#' }
"sigora_example_3"

#' Colours to use for the manually-assigned "top_pathways"
#'
#' @format A named vector of pathway types and their corresponding hex colour,
#'   from the RColorBrewer "Set2" palette
"top_pathway_colours"

#' Manually-curated list of Reactome pathways with a simple category assignment
#'
#' @format A data frame (tibble) with 1298 rows and 3 columns
#' \describe{
#' \item{pathway_id}{Reactome pathway identifier}
#' \item{pathway_name}{Pathway name}
#' \item{grouped_pathway}{Manually-curated pathway type; seven possible values}
#' }
"top_pathways_slim"

#' Manually-curated list of Reactome pathways with a simple category assignment
#'
#' @format A data frame (tibble) with 1298 rows and 6 columns
#' \describe{
#' \item{pathway_id}{Reactome pathway identifier}
#' \item{top_pathways}{Top hierarchy pathway term}
#' \item{second_top_pathway}{Seond-highest hierarchy term}
#' \item{comments}{Comments about the pathway}
#' \item{pathway_name}{Pathway name}
#' \item{grouped_pathway}{Manually-curated pathway type; seven possible values}
#' }
"top_pathways"
