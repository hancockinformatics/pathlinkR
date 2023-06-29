#' Construct a PPI network from input genes and InnateDB's database
#'
#' @param df Input data frame containing genes of interest
#' @param col Column of input genes as Ensembl IDs (character)
#' @param order Desired network order. Possible options are "zero" (default),
#'   "first," "min_simple," or "min_steiner."
#' @param hub_measure Character denoting what measure should be used in
#'   determining which nodes to highlight as hubs when plotting the network.
#'   Options include "betweenness" (default), "degree", and "hubscore". These
#'   represent network statistics calculated by their respective
#'   `tidygraph::centrality_x`, functions, specifically `degree`,
#'    `betweenness`, and `hubscore`.
#' @param ppi_data Data frame of PPI data; must contain rows of interactions as
#'   pairs of Ensembl gene IDs, with columns named "ensembl_gene_A" and
#'   "ensembl_gene_B". Defaults to pre-packaged InnateDB PPI data.
#'
#' @return `tidygraph` object for plotting or further analysis
#' @export
#'
#' @import dplyr
#' @import tidygraph
#'
#' @details The "min_steiner" method is implemented with the `SteinerNet`
#'   package.
#'
#' The "hub_measure" argument determines how `build_network` assesses
#' connectedness of nodes in the network, which will be used to highlight nodes
#' when visualizing with `plot_network`. The options are "degree",
#' "betweenness", or "hubscore". This last option uses the igraph implementation
#' of the Kleinburg hub centrality score - details on this method can be found
#' at `?igraph::hub_score`.
#'
#' @references See
#'   <https://cran.r-project.org/web/packages/SteinerNet/index.html> for details
#'   on the Steiner network trimming.
#'
#' @seealso <https://github.com/hancockinformatics/pathnet/>
#'
#' @examples
#' library(dplyr)
#'
#' ex_de_genes <- deseq_example_list[[1]] %>%
#'     tibble::rownames_to_column("gene") %>%
#'     as_tibble() %>%
#'     filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
#'
#' ppi_build_network(
#'     df = ex_de_genes,
#'     col = "gene",
#'     order = "zero"
#' )
#'
ppi_build_network <- function(
        df,
        col,
        order,
        hub_measure = "betweenness",
        ppi_data = innatedb_exp
) {

    stopifnot(is(df, "data.frame"))
    stopifnot(col %in% colnames(df))
    stopifnot(order %in% c("zero", "first", "min_simple", "min_steiner"))
    stopifnot(hub_measure %in% c("betweenness", "degree", "hubscore"))
    stopifnot(
        "'ppi_data' must have columns 'ensembl_gene_A', 'ensembl_gene_B'" = all(
            c("ensembl_gene_A", "ensembl_gene_B") %in% colnames(ppi_data)
        )
    )

    # Check for and remove any duplicate IDs
    message("Cleaning input data...")
    df_clean <- distinct(df, !!sym(col), .keep_all = TRUE)
    gene_vector <- unique(df_clean[[col]])

    stopifnot(
        "Input genes must be human Ensembl IDs" = grepl(
            x = gene_vector[1],
            pattern = "^ENSG"
        )
    )

    lost_ids <- df[[col]][duplicated(df[[col]])]

    if (length(gene_vector) < nrow(df)) {
        num_dups <- nrow(df) - length(gene_vector)

        message(
            "INFO: Found ", num_dups,
            "duplicate IDs in the input column, which have been removed:"
        )

        if (num_dups <= 10) {
            message(stringr::str_wrap(
                paste(lost_ids, collapse = ", "),
                indent = 2,
                exdent = 2
            ))
        } else {
            message(stringr::str_wrap(
                paste0(paste(lost_ids[seq_len(10)], collapse = ", "), "..."),
                indent = 2,
                exdent = 2
            ))
        }
    }

    ppi_data_ensembl <- select(ppi_data, starts_with("ensembl"))

    message("Finding interactions...")
    if (order == "zero") {
        edge_table <- ppi_data_ensembl %>% filter(
            ensembl_gene_A %in% gene_vector & ensembl_gene_B %in% gene_vector
        )
    } else {
        edge_table <- ppi_data_ensembl %>% filter(
            ensembl_gene_A %in% gene_vector | ensembl_gene_B %in% gene_vector
        )
    }

    message("Creating network...")
    network_init <- edge_table %>%
        as_tbl_graph(directed = FALSE) %>%
        ppi_remove_subnetworks() %>%
        as_tbl_graph() %>%
        mutate(
            degree = centrality_degree(),
            betweenness = centrality_betweenness(),
            seed = (name %in% gene_vector)
        ) %>%
        select(-comp)

    # Perform node filtering/trimming for minimum order networks, and
    # recalculate degree and betweenness
    if (order == "min_simple") {
        message("Performing 'simple' minimum network trimming...")

        network_out_1 <- network_init %>%
            filter(
                !(degree == 1 & !seed),
                !(betweenness == 0 & !seed)
            ) %>%
            mutate(
                degree = centrality_degree(),
                betweenness = centrality_betweenness()
            )

    } else if (order == "min_steiner") {
        message("Performing 'Steiner' minimum network trimming...")

        terminals <- network_init %>%
            activate(nodes) %>%
            pull(name) %>%
            intersect(gene_vector)

        network_out_1 <- SteinerNet::steinertree(
            type      = "SP",
            terminals = terminals,
            graph     = network_init,
            color     = FALSE
        ) %>%
            .[[1]] %>%
            as_tbl_graph(directed = FALSE) %>%
            mutate(
                degree = centrality_degree(),
                betweenness = centrality_betweenness()
            )

    } else {
        network_out_1 <- network_init
    }

    network_out_2 <-
        if (hub_measure == "betweenness") {
            network_out_1 %>% mutate(hub_score_btw = betweenness)
        } else if (hub_measure == "degree") {
            network_out_1 %>% mutate(hub_score_deg = degree)
        } else if (hub_measure == "hubscore") {
            network_out_1 %>% mutate(hub_score_hub = centrality_hub())
        }

    if (nrow(as_tibble(network_out_2)) > 2000) {
        message(
            "Your network contains more than 2000 nodes, and will likely be ",
            "difficult to interpret when plotted."
        )
    }

    message("Mapping input Ensembl IDs to HGNC symbols...")
    network_mapped <- left_join(
        network_out_2,
        select(mapping_file, "name" = ensg_id, gene_name),
        by = "name",
        multiple = "all"
    )

    network_final <- left_join(
        network_mapped,
        df_clean,
        by = c("name" = col),
        multiple = "all"
    )

    attr(network_final, "order") <- order

    message("Done.\n")
    return(network_final)
}
