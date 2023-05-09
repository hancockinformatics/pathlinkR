#' Find and return the largest subnetwork
#'
#' @param input Graph object
#'
#' @return Largest subnetwork from the input network list as an igraph object
#'
#' @importFrom igraph components induced_subgraph
#'
#' @seealso <https://github.com/hancockinformatics/pathnet/>
#'
remove_subnetworks <- function(input) {
  igraph::V(input)$comp <- components(input)$membership

  max_subnet_id <- igraph::V(input)$comp %>%
    table() %>%
    tibble::enframe() %>%
    arrange(desc(value)) %>%
    magrittr::extract2(1, 1) %>%
    as.numeric()

  induced_subgraph(input, igraph::V(input)$comp == max_subnet_id)
}
