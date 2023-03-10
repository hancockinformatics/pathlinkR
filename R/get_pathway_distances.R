#' get_pathway_distances
#'
#' @param pathway_data Three column-data frame of pathways and their constituent
#'   genes. Defaults to the provided `sigora_database` object. Must contain
#'   Ensembl gene IDs in the first column, pathway IDs in the second, and
#'   pathway descriptions in the third.
#' @param dist_method Character; method used to determine pairwise pathway
#'   distances. Can be any option supported by `vegan::vegdist()`.
#'
#' @return Matrix of the pairwise pathway distances (dissimilarity) based on
#'   overlap of their constituent genes.
#'
#' @export
#'
#' @import dplyr
#' @import purrr
#' @import stringr
#' @import tibble
#' @import tidyr
#'
#' @description Given a data frame of pathways and their member genes, calculate
#'   the pairwise distances using a constructed identity matrix.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
get_pathway_distances <- function(pathway_data = sigora_database,
                                  dist_method = "jaccard") {

  gene_id_col <- colnames(pathway_data)[
    unlist(map(pathway_data[1, ], ~str_detect(.x, "ENSG")))
  ]

  pathway_id_col <- colnames(pathway_data)[
    unlist(map(pathway_data[1, ], ~str_detect(.x, "R-[A-Z]{3}-[0-9]{1,10}")))
  ]

  message(glue::glue(
    "Using '{gene_id_col}' for gene IDs and '{pathway_id_col}' ",
    "for pathway IDs..."
  ))

  message("Creating identity matrix...")
  identify_table <- pathway_data %>%
    select(all_of(c(gene_id_col, pathway_id_col))) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(
      id_cols     = all_of(pathway_id_col),
      names_from  = all_of(gene_id_col),
      values_from = "present"
    ) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames(all_of(pathway_id_col)) %>%
    as.matrix()

  message("Running distance calculations (this may take a while)...")
  distance_matrix <- identify_table %>%
    vegan::vegdist(
      method = "jaccard",
      binary = TRUE,
      diag = TRUE
    ) %>%
    as.matrix()

  message("Done!\n")
  return(distance_matrix)
}
