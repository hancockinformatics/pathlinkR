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

#' Matrix of pairwise Jaccard indices for all human Reactome pathways
#'
#' @format A matrix with 873 rows and columns
#' \describe{
#' \item{Rows}{Reactome pathway identifiers}
#' \item{Columns}{Reactome pathway identifiers}
#' }
"pathway_distances_jaccard"

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
"sigora_result_eg"
