#' pathnet_visNetwork
#'
#' @param network Tidygraph network object as output by `create_pathnet`
#'
#' @return Interactive visNetwork plot
#' @export
#'
#' @import dplyr
#' @import visNetwork
#' @importFrom igraph as.igraph
#'
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
pathnet_visNetwork <- function(network) {
  visGroups_all <- function(x) {
    x %>%
      visGroups(groupname = names(top_pathway_colours)[1], color = top_pathway_colours[[1]]) %>%
      visGroups(groupname = names(top_pathway_colours)[2], color = top_pathway_colours[[2]]) %>%
      visGroups(groupname = names(top_pathway_colours)[3], color = top_pathway_colours[[3]]) %>%
      visGroups(groupname = names(top_pathway_colours)[4], color = top_pathway_colours[[4]]) %>%
      visGroups(groupname = names(top_pathway_colours)[5], color = top_pathway_colours[[5]]) %>%
      visGroups(groupname = names(top_pathway_colours)[6], color = top_pathway_colours[[6]]) %>%
      visGroups(groupname = names(top_pathway_colours)[7], color = top_pathway_colours[[7]])
  }

  my_igraph <- network %>%
    mutate(
      title = pathway_name_1,
      shape = if_else(is.na(bonferroni), "dot", "diamond")
    ) %>%
    rename("group" = grouped_pathway) %>%
    as.igraph()

  visIgraph(my_igraph) %>%
    visEdges(width = 6, color = "#848484") %>%
    visNodes(size = 35) %>%
    visOptions(highlightNearest = TRUE) %>%
    visLegend(position = "right", main = "Grouped pathway") %>%
    visGroups_all() %>%
    visExport(float = "left")
}
