#' create_foundation
#'
#' @param mat Matrix of distances between pathways (i.e. 0 means two pathways
#'   are identical).
#' @param cutoff Distance cutoff used to determine if two pathways should share
#'   an edge.
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
create_foundation <- function(mat, cutoff = 0.1) {

  stopifnot(all(rownames(mat) == colnames(mat)))

  mat_tibble <- mat %>%
    as.data.frame() %>%
    rownames_to_column("pathway_1") %>%
    pivot_longer(-pathway_1, names_to = "pathway_2", values_to = "distance")

  edge_table <- mat_tibble %>%
    filter(distance < cutoff)

  return(edge_table)
}
