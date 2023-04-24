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

#' Manually-curated list of Reactome pathways with a simple category assignment
#'
#' @format A data frame (tibble) with 1298 rows and 3 columns
#' \describe{
#' \item{pathway_id}{Reactome pathway identifier}
#' \item{pathway_name}{Pathway name}
#' \item{grouped_pathway}{Manually-curated pathway type; seven possible values}
#' }
"top_pathways"

#' Colours to use for the manually-assigned "top_pathways"
#'
#' @format A named vector of pathway types and their corresponding hex colour,
#'   from the RColorBrewer "Set2" palette
"top_pathway_colours"


