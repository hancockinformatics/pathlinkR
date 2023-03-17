#' create_pathnet
#'
#' @param sigora_result Data frame of results from Sigora. Must contain columns
#' "pathway_id" and "bonferroni".
#' @param foundation List of pathway pairs to use in constructing a network.
#'   Output from `create_foundation`.
#' @param trim Remove subgraphs which don't contain any enriched pathways.
#' @param trim_order Order to use when removing subgraphs.
#'
#' @return A tidygraph network object
#' @export
#'
#' @import dplyr
#' @import stringr
#' @importFrom igraph neighborhood as.igraph
#' @importFrom tidygraph tbl_graph
#' @importFrom purrr map map_chr
#'
#' @description Description will go here
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
create_pathnet <- function(sigora_result, foundation, trim = TRUE, trim_order = 1) {

  tr_trunc_neatly <- function(x, l = 60) {

    if (str_length(x) <= l) {
      return(x)
    } else {
      shortened <- x %>%
        as.character() %>%
        str_sub(., start = 1, end = l) %>%
        str_replace(., pattern = "\\s([^\\s]*)$", replacement = "...")
      return(shortened)
    }
  }

  stopifnot(all(c("pathway_id", "bonferroni") %in% colnames(sigora_result)))

  starting_nodes <- foundation %>%
    select(pathway_1, pathway_name_1) %>%
    distinct() %>%
    left_join(sigora_result, by = c("pathway_1" = "pathway_id")) %>%
    replace_na(list(bonferroni = 1))

  starting_edges <- foundation %>%
    mutate(similarity = 1 / distance) %>%
    select(pathway_1, pathway_2, similarity, distance)

  pathways_as_network <- tbl_graph(
    nodes = starting_nodes,
    edges = starting_edges,
    directed = FALSE
  ) %>%
    mutate(
      rn = row_number(),
      node_label = map_chr(
        as.character(pathway_name_1),
        ~str_wrap(.x, width = 20)
      ),
      node_label = str_replace(node_label, "^$", NA_character_)
    )

  if (trim) {
    x1 <- pathways_as_network %>%
      filter(bonferroni < 1) %>%
      pull(rn)

    valid_nodes <- map(x1, ~neighborhood(
      graph = as.igraph(pathways_as_network),
      order = trim_order,
      nodes = .x
    )) %>%
      unlist() %>%
      unique()

    pathways_as_network_trimmed <- pathways_as_network %>%
      filter(rn %in% c(x1, valid_nodes))

    return(pathways_as_network_trimmed)
  }

  return(pathways_as_network)
}
