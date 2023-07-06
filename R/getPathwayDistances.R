#' getPathwayDistances
#'
#' @param pathwayData Three column-data frame of pathways and their constituent
#'   genes. Defaults to the provided `sigoraDatabase` object. Must contain
#'   Ensembl gene IDs in the first column, pathway IDs in the second, and
#'   pathway descriptions in the third.
#' @param distMethod Character; method used to determine pairwise pathway
#'   distances. Can be any option supported by `vegan::vegdist()`.
#'
#' @return Matrix of the pairwise pathway distances (dissimilarity) based on
#'   overlap of their constituent genes.
#'
#' @export
#'
#' @import dplyr
#' @import stringr
#' @import tibble
#' @import tidyr
#' @importFrom purrr map
#'
#' @description Given a data frame of pathways and their member genes, calculate
#'   the pairwise distances using a constructed identity matrix.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' \donttest{
#'     getPathwayDistances(
#'         pathwayData = sigoraDatabase,
#'         distMethod = "jaccard"
#'     )
#' }
#'
getPathwayDistances <- function(
        pathwayData = sigoraDatabase,
        distMethod = "jaccard"
) {

    stopifnot(is(sigoraDatabase, "data.frame"))

    geneIdCol <- colnames(pathwayData)[
        unlist(map(pathwayData[1, ], ~str_detect(.x, "ENSG")))
    ]

    pathwayIdCol <- colnames(pathwayData)[
        unlist(map(
            pathwayData[1, ],
            ~str_detect(.x, "R-[A-Z]{3}-[0-9]{1,10}")
        ))
    ]

    message(
        "Using '", geneIdCol,"' for gene IDs and '",
        pathwayIdCol, "' for pathway IDs..."
    )

    message("Creating identity matrix...")
    identityTable <- pathwayData %>%
        select(all_of(c(geneIdCol, pathwayIdCol))) %>%
        distinct() %>%
        mutate(present = 1) %>%
        pivot_wider(
            id_cols     = all_of(pathwayIdCol),
            names_from  = all_of(geneIdCol),
            values_from = "present"
        ) %>%
        replace(is.na(.), 0) %>%
        column_to_rownames(all_of(pathwayIdCol)) %>%
        as.matrix()

    if (length(unique(pathwayData[[pathwayIdCol]])) > 500) {
        message("Running distance calculations (this may take a while)...")
    } else {
        message("Running distance calculations...")
    }

    distanceMatrix <- as.matrix(vegan::vegdist(
        identityTable,
        method = "jaccard",
        binary = TRUE,
        diag = TRUE
    ))

    message("Done!\n")
    return(distanceMatrix)
}