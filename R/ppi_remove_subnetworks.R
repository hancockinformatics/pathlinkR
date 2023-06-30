#' INTERNAL Find and return the largest subnetwork
#'
#' @param network Graph object
#'
#' @return Largest subnetwork from the input network list as an igraph object
#'
#' @importFrom igraph components induced_subgraph
#'
#' @seealso <https://github.com/hancockinformatics/pathnet/>
#'
ppi_remove_subnetworks <- function(network) {
    igraph::V(network)$comp <- components(network)$membership

    max_subnet_id <- igraph::V(network)$comp %>%
        table() %>%
        tibble::enframe() %>%
        arrange(desc(value)) %>%
        .[[1, 1]] %>%
        as.numeric()

    induced_subgraph(network, igraph::V(network)$comp == max_subnet_id)
}
