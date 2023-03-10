#' get_distances
#'
#' @param dist_method Character; method used to determine pairwise pathway
#'   distances. Can be any option supported by `vegan::vegdist()`.
#' @param pathway_data Data frame of pathways and their constituent genes.
#'   Defaults to the provided `sigora_database` object.
#'
#' @return Matrix of the pairwise pathway distances based on their constituent
#'   genes.
#'
#' @export
#'
#' @import cli
#' @import dplyr
#' @import purrr
#' @import readr
#' @import stringr
#' @import tibble
#' @import tidyr
#'
#' @description Given a data frame of all Reactome pathways and their genes,
#'   calculate the pairwise distances using a constructed identity matrix.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
get_distances <- function(dist_method, pathway_data = sigora_database) {


  # Helper functions --------------------------------------------------------

  path_steps_local <- function(pathway_id) {

    # Looks for the top pathway
    if (filter(HSA_react, Child == pathway_id)[1] %>% as.character != "character(0)") {

      order <- c()
      level <- pathway_id

      while (level != "character(0)") {
        hsa_id <- filter(HSA_react, Child == level)[1] %>% as.character
        path_name <- filter(reactome_names, hsa_id == level)[2] %>% as.character
        order <- c(order, path_name)
        level <- hsa_id
      }

      order <- rev(order)
      order <- paste(order, collapse = "; ")

      # Returns the top pathway name as well as the hierarchy of pathways
      # (arranged from highest to lowest)
      return(list(path_name, order))
    }

    if (pathway_id == "R-HSA-194840") {
      path_name <- "Signal Transduction"
      order <- paste0("Signal Transduction; Signaling by Rho GTPases, ",
                      "Miro GTPases and RHOBTB3; Signaling by Rho GTPases; ",
                      "RHO GTPase cycle")

      return(list(path_name, order))

    } else {
      # If already a top pathway, just return top pathway info
      path_name <- filter(reactome_names, hsa_id == pathway_id)[2] %>% as.character
      return(list(path_name, path_name))
    }
  }

  get_identity_matrix <- function(list) {
    list %>%
      map(~data.frame(id = .x)) %>%
      bind_rows(.id = "name") %>%
      mutate(present = 1) %>%
      distinct(name, id, .keep_all = TRUE) %>%
      pivot_wider(
        id_cols     = "id",
        names_from  = "name",
        values_from = "present"
      ) %>%
      replace(is.na(.), 0) %>%
      column_to_rownames(var = "id")
  }


  # First get Reactome relationships ----------------------------------------

  cli_text("Downloading Reactome data")
  reactome_db <- read_tsv(
    "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
    col_names = c("Parent", "Child"),
    col_types = cols()
  )

  reactome_names <- read_tsv(
    "https://reactome.org/download/current/ReactomePathways.txt",
    col_names = c("hsa_id", "name", "species"),
    col_types = cols()
  )


  # Add missing pathway -----------------------------------------------------

  # R-HSA-194840 is missing (possibly updated), so we'll add it manually; it's
  # "RHO GTPase cycle"
  reactome_names <- rbind(
    reactome_names,
    tibble(hsa_id  = "R-HSA-194840",
           name    = "RHO GTPase cycle",
           species = "Homo sapiens")
  )

  HSA_react <- reactome_db %>% filter(grepl("HSA", Parent))


  # Create the gene list for `get_identity_matrix()` usage --------------------------

  all_pathways <- sigora_database$pathway_id %>% unique()

  all_pathways_df <- data.frame(
    ID = NA,
    Top = NA,
    level = NA,
    hierarchy = NA
  )

  cli_progress_bar(
    "Getting pathway hierarchy information:",
    total = length(all_pathways),
    clear = FALSE
  )
  for (i in all_pathways) {
    results <- path_steps_local(i)
    hierarchy <- str_split(results[[2]], "; ") %>% unlist()
    df <- data.frame(
      ID = i,
      Top = results[[1]][1],
      level = length(hierarchy),
      hierarchy = results[[2]][1]
    )
    all_pathways_df <- rbind(all_pathways_df, df)
    cli_progress_update()
  }
  cli_progress_done()


  level_four_pathways <- all_pathways_df %>% filter(level <= 4, level != 1)
  # 873 pathways total

  pathway_description <- sigora_database %>%
    filter(pathway_id %in% level_four_pathways$ID)

  pathway_list <- list()

  cli_progress_bar(
    "Getting pathway genes:",
    total = length(level_four_pathways$ID),
    clear = FALSE
  )
  for (id in level_four_pathways$ID) {
    genes <- pathway_description %>% filter(pathway_id == id) %>%
      .$Ensembl.Gene.ID %>%
      as.character()
    list_genes <- list(id = genes)
    names(list_genes) <- id
    pathway_list <- append(pathway_list, list_genes)
    cli_progress_update()
  }
  cli_progress_done()


  # Matrix of distances between pathways ------------------------------------

  # 0 = all overlap, 1 = no overlap
  cli_text(glue::glue("Calculating distances with {dist_method} method"))
  pathway_distances_out <- pathway_list %>%
    get_identity_matrix() %>%
    t() %>%
    vegan::vegdist(
      method = dist_method,
      binary = TRUE,
      diag   = TRUE
    ) %>%
    as.matrix()

  cli_alert_success("All done!")
  return(pathway_distances_out)
}
