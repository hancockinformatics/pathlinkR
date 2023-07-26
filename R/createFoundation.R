#' Create the foundation for pathway networks
#'
#' @param mat Matrix of distances between pathways (i.e. 0 means two pathways
#'   are identical). Should match output from `getPathwayDistances`
#' @param maxDistance Distance cutoff (less than or equal) used to determine if
#'   two pathways should share an edge. Pairs with a distance of 0 are always
#'   removed. One of `maxDistance` or `propToKeep` must be provided.
#' @param propToKeep Top proportion of pathway pairs to keep as edges. One of
#'   `maxDistance` or `propToKeep` must be provided.
#'
#' @return A tibble of interacting pathway pairs
#' @export
#'
#' @import dplyr
#' @import stringr
#' @import tidyr
#' @importFrom tibble rownames_to_column
#'
#' @description From a n by n distance matrix, generate a tibble of pathway
#'   interactions to use in constructing a pathway network. The cutoff can be
#'   adjusted to have more or fewer edges in the final network, depending on the
#'   number of pathways involved.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' startingPathways <- createFoundation(
#'     mat=pathwayDistancesJaccard,
#'     maxDistance=0.8
#' )
#'
createFoundation <- function(mat, maxDistance=NA, propToKeep=NA) {

    ## Input checks
    stopifnot(all(rownames(mat) == colnames(mat)))

    matTibble <- mat %>%
        as.data.frame() %>%
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
        message(
            "Defining interactions with a distance cutoff of ",
            maxDistance,
            "..."
        )
        edgeTable <- filter(matTibble, distance <= maxDistance)
    } else if (!is.na(propToKeep)) {
        message(
            "Defining edges using the top ",
            signif(propToKeep * 100),
            "% of pathway interactions..."
        )
        edgeTable <- slice_head(matTibble, prop=propToKeep)
    }

    annoEdgeTable <- edgeTable %>%
        left_join(
            distinct(select(sigoraDatabase, pathwayId, pathwayName)),
            by=c("pathway1" = "pathwayId"),
            multiple="all"
        ) %>%
        left_join(
            distinct(select(sigoraDatabase, pathwayId, pathwayName)),
            by=c("pathway2" = "pathwayId"),
            suffix=c("1", "2"),
            multiple="all"
        ) %>%
        relocate(contains("name"), distance) %>%
        mutate(
            across(where(is.factor), as.character),
            across(where(is.character), str_trim)
        )

    message("Done! Foundation contains ", nrow(edgeTable), " interactions.")
    return(annoEdgeTable)
}
