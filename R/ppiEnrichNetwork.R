#' Test a PPI network for enriched pathways
#'
#' @param network A tidygraph network object, with Ensembl IDs in the first
#'   column of the node table
#' @param analysis Default is "sigora", but can also be "reactomepa" or
#'   "hallmark"
#' @param filterResults Should the output be filtered for significance? Use
#'   `1` to return the unfiltered results, or any number less than 1 for a
#'   custom p-value cutoff. If left as `default`, the significance cutoff
#'   for Sigora is 0.001, or 0.05 for ReactomePA and Hallmark.
#' @param gpsRepo Only applies to `analysis="sigora"`. Gene Pair Signature
#'   object for Sigora to use to test for enriched pathways. Leaving this set
#'   as "default" will use the "reaH" GPS object from `Sigora`, or you can
#'   provide your own custom GPS repository.
#' @param geneUniverse Only applies when `analysis` is "reactomepa" or
#'   "hallmark". The set of background genes to use when testing with ReactomePA
#'   or Hallmark gene sets. For ReactomePA this must be a character vector of
#'   Entrez genes. For Hallmark, it must be Ensembl IDs.
#'
#' @return A data frame (tibble) of enriched pathways
#' @export
#'
#' @import dplyr
#'
#' @importFrom tibble column_to_rownames
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' exNetwork <- ppiBuildNetwork(
#'     deseqResults=deseqExampleList[[1]],
#'     filterInput=TRUE,
#'     order="zero"
#' )
#'
#' ppiEnrichNetwork(
#'     network=exNetwork,
#'     analysis="sigora",
#'     gpsRepo="default"
#' )
#'
ppiEnrichNetwork <- function(
        network,
        analysis="sigora",
        filterResults="default",
        gpsRepo="default",
        geneUniverse=NULL
) {

    networkTable <- network %>%
        tibble::as_tibble() %>%
        column_to_rownames("name")

    newList <- list("network" = networkTable)

    pathwayEnrichment(
        inputList=newList,
        analysis=analysis,
        filterInput=FALSE,
        split=FALSE,
        filterResults=filterResults,
        gpsRepo=gpsRepo,
        geneUniverse=geneUniverse
    )
}
