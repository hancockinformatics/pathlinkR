#' create_foundation
#'
#' @param mat Matrix of distances between pathways (i.e. 0 means two pathways
#'   are identical).
#' @param max_distance Distance cutoff (less than or equal) used to determine if
#'   two pathways should share an edge. Pairs with a distance of 0 are always
#'   removed.
#'
#' @return A tibble
#' @export
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#'
#' @description From a n by n distance matrix, generate a three column tibble to
#'   use in constructing a pathway network. The cutoff can be adjusted to have
#'   more or fewer edges in the final network, depending on the number of
#'   pathways involved.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
create_foundation <- function(mat, max_distance = 0.1) {

  stopifnot(all(rownames(mat) == colnames(mat)))

  mat_tibble <- mat %>%
    as.data.frame() %>%
    rownames_to_column("pathway_1") %>%
    pivot_longer(-pathway_1, names_to = "pathway_2", values_to = "distance")

  edge_table <- mat_tibble %>%
    filter(distance != 0, distance <= max_distance)

  return(edge_table)
}
