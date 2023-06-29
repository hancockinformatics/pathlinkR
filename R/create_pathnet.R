#' create_pathnet
#'
#' @param sigora_result Data frame of results from Sigora. Must contain columns
#' "pathway_id" and "p_value_adjusted".
#' @param foundation List of pathway pairs to use in constructing a network,
#'   output from `create_foundation`.
#' @param trim Remove subgraphs which don't contain any enriched pathways
#'   (default is `TRUE`).
#' @param trim_order Order to use when removing subgraphs; Higher values will
#'   keep more non-enriched pathway nodes. Defaults to `1`.
#'
#' @return A tidygraph network object
#' @export
#'
#' @import dplyr
#' @import stringr
#' @importFrom igraph neighborhood as.igraph
#' @importFrom tidygraph tbl_graph activate
#' @importFrom purrr map map_chr
#'
#' @description Creates a tidygraph network object from the pathway information,
#'   ready to be visualized with `plot_pathnet`.
#'
#' @details With the "trim" option enabled, nodes (pathways) and subgraphs which
#'   are not sufficiently connected to enriched pathways will be removed. How
#'   aggressively this is done can be controlled via the `trim_order` argument.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
#' @examples
#' ex_starting_pathways <- create_foundation(
#'     mat = pathway_distances_jaccard,
#'     max_distance = 0.8
#' )
#'
#' ex_pathnet <- create_pathnet(
#'     sigora_result = sigora_examples,
#'     foundation = ex_starting_pathways,
#'     trim = TRUE,
#'     trim_order = 1
#' )
#'
create_pathnet <- function(
        sigora_result,
        foundation,
        trim = TRUE,
        trim_order = 1
) {

    stopifnot(
        all(c("pathway_id", "p_value_adjusted") %in% colnames(sigora_result))
    )

    starting_nodes <- foundation %>%
        select(pathway_1, pathway_name_1) %>%
        distinct() %>%
        left_join(
            sigora_result,
            by = c("pathway_1" = "pathway_id"),
            multiple = "all"
        )

    starting_edges <- foundation %>%
        mutate(similarity = 1 / distance) %>%
        select(pathway_1, pathway_2, similarity, distance)

    pathways_as_network <- tbl_graph(
        nodes = starting_nodes,
        edges = starting_edges,
        directed = FALSE
    ) %>%
        mutate(rn = row_number())

    pathways_as_network_2 <-
        if (trim) {
            x1 <- pathways_as_network %>%
                filter(!is.na(p_value_adjusted)) %>%
                pull(rn)

            valid_nodes <- map(x1, ~neighborhood(
                graph = as.igraph(pathways_as_network),
                order = trim_order,
                nodes = .x
            )) %>%
                unlist() %>%
                unique()

            pathways_as_network %>%
                filter(rn %in% c(x1, valid_nodes)) %>%
                select(-rn)
        } else {
            select(pathways_as_network, -rn)
        }

    pathways_as_network_3 <- pathways_as_network_2 %>%
        activate("edges") %>%
        distinct() %>%
        activate("nodes")

    pathways_as_network_4 <- pathways_as_network_3 %>%
        left_join(
            x  = .,
            y  = select(
                top_pathways_more, pathway_id,
                pathway_name, grouped_pathway
            ),
            by = c("pathway_1" = "pathway_id"),
            multiple = "all"
        ) %>%
        select(
            pathway_1,
            pathway_name_1,
            everything(),
            -any_of(c("pathway_name", "level_1", "level_2"))
        )

    return(pathways_as_network_4)
}
