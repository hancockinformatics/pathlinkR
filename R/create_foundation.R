#' create_foundation
#'
#' @param mat Matrix of distances between pathways (i.e. 0 means two pathways
#'   are identical).
#' @param max_distance Distance cutoff (less than or equal) used to determine if
#'   two pathways should share an edge. Pairs with a distance of 0 are always
#'   removed. One of `max_distance` or `prop_to_keep` must be provided.
#' @param prop_to_keep Top proportion of pathway pairs to keep as edges. One of
#'   `max_distance` or `prop_to_keep` must be provided.
#'
#' @return A tibble
#' @export
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#'
#' @description From a n by n distance matrix, generate a tibble to use in
#'   constructing a pathway network. The cutoff can be adjusted to have more or
#'   fewer edges in the final network, depending on the number of pathways
#'   involved.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
create_foundation <- function(mat, max_distance = NA, prop_to_keep = NA) {

  stopifnot(all(rownames(mat) == colnames(mat)))

  mat_tibble <- mat %>%
    as.data.frame() %>%
    rownames_to_column("pathway_1") %>%
    pivot_longer(-pathway_1, names_to = "pathway_2", values_to = "distance") %>%
    distinct() %>%
    filter(distance != 0) %>%
    arrange(distance)

  if (!is.na(max_distance)) {
    edge_table <- filter(mat_tibble, distance <= 0.5)
  } else if (!is.na(prop_to_keep)) {
    edge_table <- slice_head(mat_tibble, prop = prop_to_keep)
  }

  anno_edge_table <- edge_table %>%
    left_join(
      distinct(select(sigora_database, pathway_id, pathway_name)),
      by = c("pathway_1" = "pathway_id")
    ) %>%
    left_join(
      distinct(select(sigora_database, pathway_id, pathway_name)),
      by = c("pathway_2" = "pathway_id"),
      suffix = c("_1", "_2")
    ) %>%
    relocate(contains("name"), distance)

  return(anno_edge_table)
}
