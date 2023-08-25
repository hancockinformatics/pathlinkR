#' Create the foundation for pathway networks using pathway distances
#'
#' @param mat Matrix of distances between pathways, i.e. 0 means two pathways
#'   are identical. Should match the output from `getPathwayDistances`.
#' @param maxDistance Numeric distance cutoff (less than or equal) used to
#'   determine if two pathways should share an edge. Pathway pairs with a
#'   distance of 0 are always removed. One of `maxDistance` or
#'   `propToKeep` must be provided.
#' @param propToKeep Top proportion of pathway pairs to keep as edges, ranked
#'   based distance. One of `maxDistance` or `propToKeep` must be
#'   provided.
#'
#' @return A tibble of interacting pathway pairs
#' @export
#'
#' @import dplyr
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#'
#' @description From a "n by n" distance matrix, generate a table of interacting
#'   pathways to use in constructing a pathway network. The cutoff can be
#'   adjusted to have more or fewer edges in the final network, depending on the
#'   number of pathways involved, i.e. the number of enriched pathways you're
#'   trying to visualize.
#'
#'  The desired cutoff will also vary based on the distance measure used, so
#'  some trial-and-error may be needed to find an appropriate value.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' startingPathways <- pathnetFoundation(
#'     mat=pathwayDistancesJaccard,
#'     maxDistance=0.8
#' )
#'
pathnetFoundation <- function(mat, maxDistance=NA, propToKeep=NA) {
    stopifnot(all(rownames(mat) == colnames(mat)))

    matTibble <- as.data.frame(mat) %>%
        rownames_to_column("pathway1") %>%
        pivot_longer(
            -pathway1,
            names_to="pathway2",
            values_to="distance"
        ) %>%
        distinct() %>%
        filter(distance != 0) %>%
        arrange(distance) %>%
        mutate(across(where(is.factor), as.character))

    if (!is.na(maxDistance)) {
        edgeTable <- filter(matTibble, distance <= maxDistance)
    } else if (!is.na(propToKeep)) {
        edgeTable <- slice_head(matTibble, prop=propToKeep)
    }

    annoEdgeTable <- edgeTable %>%
        left_join(
            distinct(sigoraDatabase, pathwayId, pathwayName),
            by=c("pathway1" = "pathwayId"),
            multiple="all"
        ) %>%
        left_join(
            distinct(sigoraDatabase, pathwayId, pathwayName),
            by=c("pathway2" = "pathwayId"),
            suffix=c("1", "2"),
            multiple="all"
        ) %>%
        relocate(contains("name"), distance) %>%
        mutate(
            across(where(is.factor), as.character),
            across(where(is.character), trimws)
        )

    return(annoEdgeTable)
}
