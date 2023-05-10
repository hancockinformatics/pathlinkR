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
#'   `tidygraph::centrality_x`, functions, specifically `degree`, `betweenness`,
#'   and `hub_score`.
#' @param ppi_data Data frame of PPI data; must contain rows of interactions as
#'   pairs of Ensembl gene IDs, with columns named "ensembl_gene_A" and
#'   "ensembl_gene_B". Defaults to pre-packaged InnateDB PPI data.
#' @param seed Number used in call to `set.seed()` to allow for reproducible
#'   network generation
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
build_network <- function(df,
                          col,
                          order,
                          hub_measure = "betweenness",
                          ppi_data = innatedb_exp,
                          seed = 1) {

  # Check for and remove any duplicate IDs, which will cause problems later.
  # Make sure to warn the user about this.
  message("Cleaning input data...")
  df_clean <- distinct(df, !!sym(col), .keep_all = TRUE)
  gene_vector <- unique(df_clean[[col]])
  lost_ids <- df[[col]][duplicated(df[[col]])]

  if (length(gene_vector) < nrow(df)) {

    num_dups <- nrow(df) - length(gene_vector)

    message(paste0(
      "  INFO: Found ",
      num_dups,
      " duplicate IDs in the input column, which have been removed:"
    ))

    if (num_dups <= 10) {
      message(stringr::str_wrap(
        paste(lost_ids, collapse = ", "),
        indent = 2,
        exdent = 2
      ))
    } else {
      message(stringr::str_wrap(
        paste0(paste(lost_ids[1:10], collapse = ", "), "..."),
        indent = 2,
        exdent = 2
      ))
    }
  }

  if (!grepl(x = gene_vector[1], pattern = "^ENSG")) {
    stop("Input genes must be human Ensembl IDs")
  }

  if ( !all(c("ensembl_gene_A", "ensembl_gene_B") %in% colnames(ppi_data)) ) {
    stop("Argument 'ppi_data' must be a data frame containing columns ",
         "'ensembl_gene_A' and 'ensembl_gene_B' (case sensitive)")
  }

  ppi_data_ensembl <- ppi_data %>%
    dplyr::select(starts_with("ensembl"))

  message("Finding interactions...")
  if (order == "zero") {
    edge_table <- ppi_data_ensembl %>%
      filter(ensembl_gene_A %in% gene_vector & ensembl_gene_B %in% gene_vector)
  } else if (order %in% c("first", "min_simple", "min_steiner")) {
    edge_table <- ppi_data_ensembl %>%
      filter(ensembl_gene_A %in% gene_vector | ensembl_gene_B %in% gene_vector)
  } else {
    stop(
      "Argument 'order' must be one of: ",
      "'zero', 'first', 'min_simple', or 'min_steiner'"
    )
  }

  if (hub_measure == "betweenness") {
    hub_fn <- centrality_betweenness
  } else if (hub_measure == "degree") {
    hub_fn <- centrality_degree
  } else if (hub_measure == "hubscore") {
    hub_fn <- centrality_hub
  } else {
    stop(paste0(
      "Argument 'hub_measure' must be one of 'betweenness', 'degree', or ",
      "'hubscore'"
    ))
  }

  message("Creating network...")
  network_init_1 <- edge_table %>%
    as_tbl_graph(directed = FALSE) %>%
    remove_subnetworks() %>%
    as_tbl_graph() %>%
    mutate(
      degree      = centrality_degree(),
      betweenness = centrality_betweenness(),
      seed        = (name %in% gene_vector)
    ) %>%
    dplyr::select(-comp)

  network_init_2 <-
    if (hub_measure == "betweenness") {
      network_init_1 %>% mutate(hub_score_btw = betweenness)
    } else if (hub_measure == "degree") {
      network_init_1 %>% mutate(hub_score_deg = degree)
    } else if (hub_measure == "hubscore") {
      network_init_1 %>% mutate(hub_score_hub = centrality_hub())
    }


  # Perform node filtering/trimming for minimum order networks, and recalculate
  # the network statistics
  if (order == "min_simple") {

    message("Performing 'simple' minimum network trimming...")
    network_out_1 <- network_init_2 %>%
      filter(
        !(degree == 1 & !seed),
        !(betweenness == 0 & !seed)
      ) %>%
      mutate(
        degree = centrality_degree(),
        betweenness = centrality_betweenness()
      )

    network_out_2 <-
      if (hub_measure == "betweenness") {
        network_out_1 %>% mutate(hub_score_btw = betweenness)
      } else if (hub_measure == "degree") {
        network_out_1 %>% mutate(hub_score_deg = degree)
      } else if (hub_measure == "hubscore") {
        network_out_1 %>% mutate(hub_score_hub = centrality_hub())
      }

  } else if (order == "min_steiner") {

    message("Performing 'Steiner' minimum network trimming...")
    set.seed(seed)

    terminals <- network_init_2 %>%
      activate(nodes) %>%
      pull(name) %>%
      intersect(gene_vector)

    network_out_1 <- SteinerNet::steinertree(
      type      = "SP",
      terminals = terminals,
      graph     = network_init_2,
      color     = FALSE
    ) %>%
      magrittr::extract2(1) %>%
      as_tbl_graph(directed = FALSE) %>%
      mutate(
        degree = centrality_degree(),
        betweenness = centrality_betweenness()
      )

    network_out_2 <-
      if (hub_measure == "betweenness") {
        network_out_1 %>% mutate(hub_score_btw = betweenness)
      } else if (hub_measure == "degree") {
        network_out_1 %>% mutate(hub_score_deg = degree)
      } else if (hub_measure == "hubscore") {
        network_out_1 %>% mutate(hub_score_hub = centrality_hub())
      }

  } else {
    network_out_2 <- network_init_2
  }

  if (nrow(as_tibble(network_out_2)) > 2000) {
    message(
      "Warning:\nYour network contains more than 2000 nodes, and will likely ",
      "be difficult to interpret when plotted."
    )
  } else if (nrow(as_tibble(network_out_2)) > 2000 & order != "zero") {
    message(
      "\nWarning:\nYour network contains more than 2000 nodes, and will ",
      "likely be difficult to interpret when plotted. Consider switching to a ",
      "zero order network to improve legibility.\n"
    )
  }

  message("Mapping input Ensembl IDs to HGNC symbols...")
  ensembl_to_hgnc <- mapping_file %>%
    dplyr::select("name" = ensg_id, gene_name)

  network_mapped <- left_join(
    network_out_2,
    ensembl_to_hgnc,
    by = "name"
  )

  network_final <- left_join(
    network_mapped,
    df_clean,
    by = c("name" = col)
  )

  attr(network_final, "order") <- order

  message("Done.\n")
  return(network_final)
}
