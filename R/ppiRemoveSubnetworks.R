#' INTERNAL Find and return the largest subnetwork
#'
#' @param network Graph object
#'
#' @return Largest subnetwork from the input network list as an igraph object
#'
#' @importFrom igraph V components induced_subgraph
#' @importFrom tibble enframe
#'
#' @seealso <https://github.com/hancockinformatics/pathnet/>
#'
ppiRemoveSubnetworks <- function(network) {
    igraph::V(network)$comp <- components(network)$membership

    maxSubnetId <- igraph::V(network)$comp %>%
        table() %>%
        enframe() %>%
        arrange(desc(value)) %>%
        .[[1, 1]] %>%
        as.numeric()

    induced_subgraph(network, igraph::V(network)$comp == maxSubnetId)
}
