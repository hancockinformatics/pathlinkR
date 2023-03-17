create_pathnet <- function(sigora_result, foundation) {

  stopifnot(all(c("pathway_id", "bonferroni") %in% colnames(sigora_result)))

  starting_nodes <- foundation %>%
    select(pathway_1, pathway_name_1) %>%
    distinct() %>%
    left_join(sigora_result, by = c("pathway_1" = "pathway_id")) %>%
    replace_na(list(bonferroni = 1))

  starting_edges <- foundation %>%
    mutate(similarity = 1 / distance) %>%
    select(pathway_1, pathway_2, similarity, distance)

  pathways_as_network <- tidygraph::tbl_graph(
    nodes = starting_nodes,
    edges = starting_edges,
    directed = FALSE
  ) %>%
    mutate(
      rn = row_number(),
      node_label = if_else(is.na(direction), "", pathway_name_1),
      node_label = map_chr(node_label, ~tRavis::tr_trunc_neatly(.x, l = 50)),
      node_label = str_replace(node_label, "^$", NA_character_)
    )

  return(pathways_as_network)
}
