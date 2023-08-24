#' Calculate pairwise distances from a table of pathways and genes
#'
#' @param pathwayData Three column data frame of pathways and their constituent
#'   genes. Defaults to the provided `sigoraDatabase` object, but can be any set
#'   of Reactome pathways. Must contain Ensembl gene IDs in the first column,
#'   human Reactome pathway IDs in the second, and pathway descriptions in the
#'   third.
#' @param distMethod Character; method used to determine pairwise pathway
#'   distances. Can be any option supported by `vegan::vegdist()`.
#'
#' @return Matrix of the pairwise pathway distances (dissimilarity) based on
#'   overlap of their constituent genes.
#' @export
#'
#' @importFrom dplyr %>% all_of distinct mutate select
#' @importFrom purrr map
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider
#' @importFrom vegan vegdist
#'
#' @description Given a data frame of pathways and their member genes, calculate
#'   the pairwise distances using a constructed identity matrix. Zero means two
#'   pathways are identical, while one means two pathways share no genes in
#'   common.
#'
#'   `pathlinkR` includes an example distance object
#'   (`pathwayDistancesJaccard`), created using the included
#'   "sigoraDatabase" object and the Jaccard distance measure.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' # Here we'll use a subset of pathways to save time
#' getPathwayDistances(
#'     pathwayData=dplyr::slice_head(
#'         dplyr::arrange(sigoraDatabase, pathwayId),
#'         prop=0.25
#'     ),
#'     distMethod="jaccard"
#' )
#'
getPathwayDistances <- function(
        pathwayData=sigoraDatabase,
        distMethod="jaccard"
) {
    stopifnot(is(pathwayData, "data.frame"))

    ## Identify which columns have Ensembl and pathway IDs
    geneIdCol <- colnames(pathwayData)[
        unlist(map(pathwayData[1, ], ~grepl(x=.x, pattern="ENSG")))
    ]

    pathwayIdCol <- colnames(pathwayData)[
        unlist(map(
            pathwayData[1, ],
            ~grepl(x=.x, pattern="R-[A-Z]{3}-[0-9]{1,10}")
        ))
    ]

    stopifnot(
        "Couldn't find a column of Ensembl gene IDs" = length(geneIdCol) > 0
    )
    stopifnot(
        "Couldn't find a column of pathway IDs" = length(pathwayIdCol) > 0
    )

    message(
        "Using '", geneIdCol,"' for gene IDs and '",
        pathwayIdCol, "' for pathway IDs..."
    )

    identityTable <- pathwayData %>%
        select(all_of(c(geneIdCol, pathwayIdCol))) %>%
        distinct() %>%
        mutate(present=1) %>%
        pivot_wider(
            id_cols=all_of(pathwayIdCol),
            names_from =all_of(geneIdCol),
            values_from="present"
        ) %>%
        replace(is.na(.), 0) %>%
        column_to_rownames(pathwayIdCol) %>%
        as.matrix()

    if (length(unique(pathwayData[[pathwayIdCol]])) > 500) {
        message("Running distance calculations (this may take a while)...")
    }

    distanceMatrix <- as.matrix(vegdist(
        identityTable,
        method="jaccard",
        binary=TRUE,
        diag=TRUE
    ))

    message("Done!\n")
    return(distanceMatrix)
}
