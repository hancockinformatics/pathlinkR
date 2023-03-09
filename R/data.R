#' Table of all Sigora pathways and their constituent genes
#'
#' @format A tibble with 71450 rows and 5 columns
#' \describe{
#' \item{EntrezGene.ID}{Entrez gene identifier}
#' \item{EnsemblGene.ID}{Ensembl gene identifier}
#' \item{Symbol}{HGNC gene symbol}
#' \item{pathway_id}{Reactome pathway identifier}
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
