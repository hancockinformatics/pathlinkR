#' List of example results from DESeq2
#'
#' @format A list of two data frames, each with 6 columns, containing
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

#' Colour assignments for grouped pathways
#' @format A length 8 named vector of hex colour values
"grouped_pathway_colours"

#' Human gene ID mapping from Sigora
#'
#' A data frame of ID mappings for different human ID types. From the Sigora
#' package.
#'
#' @format A data frame with 47037 rows and 3 columns
#' \describe{
#'   \item{Ensembl.Gene.ID}{Ensembl gene ID}
#'   \item{EntrezGene.ID}{Entrez/NCBI gene ID}
#'   \item{Symbol}{HGNC symbol}
#' }
"idmap"

#' Experimentally verified human PPI data downloaded from InnateDB
#'
#' A data frame containing human PPI data from InnateDB, from the entry
#' "All Experimentally Validated Interactions (updated weekly)" at
#' <https://innatedb.com/redirect.do?go=downloadImported>.
#' A few important steps have been taken to filter the data, namely the removal
#' of duplicate interactions, and removing interactions that have the same
#' components but are swapped between A and B.
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

#' Table of Reactome pathways that had multiple top pathways originally, and
#' manually curated to only have one top pathway in find_top_pathways.R
#'
#' @format A tibble with 25 rows and 2 columns
#' \describe{
#' \item{pathway_id}{Reactome pathway ID}
#' \item{pathway_name}{Name of the pathway}
#' \item{top_pathway_name}{Name of the top pathway}
#' \item{species}{species}
#' \item{top_pathway}{ID of the top pathway}
#' }
"manual_dupe_annotation"

#' Table of human gene ID mappings
#'
#' @format A data frame with 42,099 rows and 3 columns
#' \describe{
#' \item{ensg_id}{Ensembl IDs}
#' \item{gene_name}{HGNC symbols}
#' \item{entrez_id}{NCBI Entrez IDs}
#' }
"mapping_file"

#' Table of Hallmark gene sets and their genes
#'
#' @format A data frame with 8,209 rows and 2 columns
#' \describe{
#' \item{pathway_id}{Name of the Hallmark Gene Set}
#' \item{ensg_id}{Ensembl IDs}
#' }
"msigdbr_t2g"

#' Matrix of pairwise Jaccard indices for all human Reactome pathways
#'
#' @format A matrix with 1001 rows and columns
#' \describe{
#' \item{Rows}{Reactome pathway identifiers}
#' \item{Columns}{Reactome pathway identifiers}
#' }
"pathway_distances_jaccard"

#' Table of all Reactome pathways and genes
#'
#' @format A data frame with 123519 rows and 3 columns
#' \describe{
#'   \item{pathway_id}{Reactome pathway ID}
#'   \item{entrez_id}{Entrez gene ID}
#'   \item{pathway_name}{Name of the Reactome pathway}
#' }
"reactome_database"

#' The default Sigora Gene Pair Signature (GPS) object
#'
#' @format A list of length nine with the following top-level elements
#' \describe{
#'   \item{origRepo}{The original pathway data}
#'   \item{L1}{Level 1 pathways}
#'   \item{L2}{Level 2 pathways}
#'   \item{L3}{Level 3 pathways}
#'   \item{L4}{Level 4 pathways}
#'   \item{L5}{Level 5 pathways}
#'   \item{repoName}{Name of the repository used to build the GPS object}
#'   \item{pathwaydescriptions}{Pathway names}
#'   \item{call}{Function call to create the GPS object}
#' }
"reaH"

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

#' Example outut from `enrich_pathways` run on "deseq_example_list"
#'
#' @format A data frame (tibble) with 66 rows and 12 columns
#' \describe{
#' \item{pathway_id}{Reactome pathway identifier}
#' \item{pathway_description}{Description of the pathway}
#' \item{direction}{Was the pathway enriched in up or down regulated genes}
#' \item{p_value}{Nominal p value for the enrichment}
#' \item{p_value_adjusted}{p value adjusted for multiple testing}
#' \item{genes}{Genes in the pathway/input}
#' \item{num_candidate_genes}{Analyzed genes found in the pathway of interest}
#' \item{num_bg_genes}{All genes from the pathway database}
#' \item{gene_ratio}{Quotient of the number of candidate and background genes}
#' \item{top_pathways}{Pathway category}
#' \item{comparison}{Comparison from which results are dervied; name of input list}
#' \item{total_genes}{Total number of input genes}
#' }
"sigora_examples"

#' Manually-curated list of Reactome and Hallmark pathways and their top
#' pathways and grouped pathways
#'
#' @format A data frame (tibble) with 2671 rows and 5 columns
#' \describe{
#' \item{pathway_id}{Reactome or Hallmark pathway identifier}
#' \item{top_pathways}{Top hierarchy pathway term, shortened in some cases}
#' \item{pathway_name}{Pathway name}
#' \item{grouped_pathway}{Top grouped pathway, 8 for Reactome}
#' \item{top_pathways_original}{Original top pathway name}
#' }
"top_pathways_more"

#' Manually-curated Reactome pathways with a simple category assignment
#'
#' @format A data frame (tibble) with 1298 rows and 6 columns
#' \describe{
#' \item{pathway_id}{Reactome pathway identifier}
#' \item{top_pathways}{Top hierarchy pathway term}
#' \item{second_top_pathway}{Seond-highest hierarchy term}
#' \item{comments}{Comments about the pathway}
#' \item{pathway_name}{Pathway name}
#' \item{grouped_pathway}{Manually-curated pathway type; seven possible values}
#' \item{top_pathways_original}{Original top pathway}
#' }
"top_pathways"
