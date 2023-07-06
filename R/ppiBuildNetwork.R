#' Construct a PPI network from input genes and InnateDB's database
#'
#' @param df Input data frame containing genes of interest
#' @param col Column of input genes as Ensembl IDs (character)
#' @param order Desired network order. Possible options are "zero" (default),
#'   "first," "minSimple," or "minSteiner."
#' @param hubMeasure Character denoting what measure should be used in
#'   determining which nodes to highlight as hubs when plotting the network.
#'   Options include "betweenness" (default), "degree", and "hubscore". These
#'   represent network statistics calculated by their respective
#'   `tidygraph::centrality_x`, functions, specifically `degree`,
#'    `betweenness`, and `hubscore`.
#' @param ppiData Data frame of PPI data; must contain rows of interactions as
#'   pairs of Ensembl gene IDs, with columns named "ensemblGeneA" and
#'   "ensemblGeneB". Defaults to pre-packaged InnateDB PPI data.
#'
#' @return `tidygraph` object for plotting or further analysis
#' @export
#'
#' @import dplyr
#' @import tidygraph
#'
#' @details The "minSteiner" method is implemented with the `SteinerNet`
#'   package.
#'
#' The "hubMeasure" argument determines how `ppiBuildNetwork` assesses
#' connectedness of nodes in the network, which will be used to highlight nodes
#' when visualizing with `ppiPlotNetwork`. The options are "degree",
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
#' exDEGenes <- dplyr::filter(
#'     dplyr::as_tibble(tibble::rownames_to_column(
#'         deseqExampleList[[1]],
#'         "gene"
#'     )),
#'     padj < 0.05,
#'     abs(log2FoldChange) > log2(1.5)
#' )
#'
#' ppiBuildNetwork(
#'     df = exDEGenes,
#'     col = "gene",
#'     order = "zero"
#' )
#'
ppiBuildNetwork <- function(
        df,
        col,
        order,
        hubMeasure = "betweenness",
        ppiData = innateDbExp
) {

    stopifnot(is(df, "data.frame"))
    stopifnot(col %in% colnames(df))
    stopifnot(order %in% c("zero", "first", "minSimple", "minSteiner"))
    stopifnot(hubMeasure %in% c("betweenness", "degree", "hubscore"))
    stopifnot(
        "'ppiData' must have columns 'ensemblGeneA', 'ensemblGeneB'" = all(
            c("ensemblGeneA", "ensemblGeneB") %in% colnames(ppiData)
        )
    )

    ## Check for and remove any duplicate IDs, warning the user when this occurs
    message("Cleaning input data...")
    dfClean <- distinct(df, !!sym(col), .keep_all = TRUE)
    geneVector <- unique(dfClean[[col]])

    stopifnot(
        "Input genes must be human Ensembl IDs" = grepl(
            x = geneVector[1],
            pattern = "^ENSG"
        )
    )

    lostIds <- df[[col]][duplicated(df[[col]])]

    if (length(geneVector) < nrow(df)) {
        numDups <- nrow(df) - length(geneVector)

        message(
            "INFO: Found ", numDups,
            "duplicate IDs in the input column, which have been removed:"
        )

        if (numDups <= 10) {
            message(stringr::str_wrap(
                paste(lostIds, collapse = ", "),
                indent = 2,
                exdent = 2
            ))
        } else {
            message(stringr::str_wrap(
                paste0(paste(lostIds[seq_len(10)], collapse = ", "), "..."),
                indent = 2,
                exdent = 2
            ))
        }
    }

    ppiDataEnsembl <- select(ppiData, starts_with("ensembl"))

    message("Finding interactions...")
    if (order == "zero") {
        edgeTable <- ppiDataEnsembl %>% filter(
            ensemblGeneA %in% geneVector & ensemblGeneB %in% geneVector
        )
    } else {
        edgeTable <- ppiDataEnsembl %>% filter(
            ensemblGeneA %in% geneVector | ensemblGeneB %in% geneVector
        )
    }

    message("Creating network...")
    networkInit <- edgeTable %>%
        as_tbl_graph(directed = FALSE) %>%
        ppiRemoveSubnetworks() %>%
        as_tbl_graph() %>%
        mutate(
            degree = centrality_degree(),
            betweenness = centrality_betweenness(),
            seed = (name %in% geneVector)
        ) %>%
        select(-comp)

    ## Perform node filtering/trimming for minimum order networks, and
    ## recalculate degree and betweenness
    if (order == "minSimple") {
        message("Performing 'simple' minimum network trimming...")

        networkOut1 <- networkInit %>%
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

        terminals <- networkInit %>%
            activate(nodes) %>%
            pull(name) %>%
            intersect(geneVector)

        networkOut1 <- SteinerNet::steinertree(
            type      = "SP",
            terminals = terminals,
            graph     = networkInit,
            color     = FALSE
        ) %>%
            .[[1]] %>%
            as_tbl_graph(directed = FALSE) %>%
            mutate(
                degree = centrality_degree(),
                betweenness = centrality_betweenness()
            )

    } else {
        networkOut1 <- networkInit
    }

    networkOut2 <-
        if (hubMeasure == "betweenness") {
            networkOut1 %>% mutate(hubScoreBtw = betweenness)
        } else if (hubMeasure == "degree") {
            networkOut1 %>% mutate(hubScoreDeg = degree)
        } else if (hubMeasure == "hubscore") {
            networkOut1 %>% mutate(hubScoreHub = centrality_hub())
        }

    if (nrow(as_tibble(networkOut2)) > 2000) {
        message(
            "Your network contains more than 2000 nodes, and will likely be ",
            "difficult to interpret when plotted."
        )
    }

    message("Mapping input Ensembl IDs to HGNC symbols...")
    networkMapped <- left_join(
        networkOut2,
        select(mappingFile, "name" = ensemblGeneId, hgncSymbol),
        by = "name",
        multiple = "all"
    )

    networkFinal <- left_join(
        networkMapped,
        dfClean,
        by = c("name" = col),
        multiple = "all"
    )

    attr(networkFinal, "order") <- order

    message("Done.\n")
    return(networkFinal)
}
