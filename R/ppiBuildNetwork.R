#' Construct a PPI network from input genes and InnateDB's database
#'
#' @param deseqResults Data frame of DESeq2 results with Ensembl gene IDs as
#'   rownames
#' @param filterInput If providing list of data frames containing the
#'   unfiltered output from `DESeq2::results()`, set this to TRUE to filter for
#'   DE genes using the thresholds set by the `pCutoff` and `fcCutoff`
#'   arguments. When FALSE it's assumed your passing the filtered
#'   results into `inputList` and no more filtering will be done.
#' @param pCutoff Adjusted p value cutoff, defaults to <0.05
#' @param fcCutoff Absolute fold change cutoff, defaults to an absolute value
#'   of >1.5
#' @param order Desired network order. Possible options are "zero" (default),
#'   "first," "minSimple," or "minSteiner."
#' @param hubMeasure Character denoting what measure should be used in
#'   determining which nodes to highlight as hubs when plotting the network.
#'   Options include "betweenness" (default), "degree", and "hubscore". These
#'   represent network statistics calculated by their respective
#'   `tidygraph::centrality_x`, functions.
#' @param ppiData Data frame of PPI data; must contain rows of interactions as
#'   pairs of Ensembl gene IDs, with columns named "ensemblGeneA" and
#'   "ensemblGeneB". Defaults to pre-packaged InnateDB PPI data.
#'
#' @return `tidygraph` object for plotting or further analysis
#' @export
#'
#' @import dplyr
#' @import stringr
#' @import tidygraph
#' @importFrom SteinerNet steinertree
#'
#' @description Creates a protein-protein interaction (PPI) network using
#'   data from InnateDB, with options for network order, and filtering input.
#'
#' @details The "minSteiner" method is implemented with the `SteinerNet`
#'   package.
#'
#'   The "hubMeasure" argument determines how `ppiBuildNetwork` assesses
#'   connectedness of nodes in the network, which will be used to highlight
#'   nodes when visualizing with `ppiPlotNetwork`. The options are "degree",
#'   "betweenness", or "hubscore". This last option uses the igraph
#'   implementation of the Kleinburg hub centrality score - details on this
#'   method can be found at `?igraph::hub_score`.
#'
#' @references See
#'   <https://cran.r-project.org/web/packages/SteinerNet/index.html> for details
#'   on the Steiner network trimming.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR/>
#'
#' @examples
#' ppiBuildNetwork(
#'     deseqResults=deseqExampleList[[1]],
#'     filterInput=TRUE,
#'     order="zero"
#' )
#'
ppiBuildNetwork <- function(
        deseqResults,
        filterInput=TRUE,
        pCutoff=0.05,
        fcCutoff=1.5,
        order="zero",
        hubMeasure="betweenness",
        ppiData=innateDbExp
) {

    stopifnot(is(deseqResults, "data.frame"))
    stopifnot(grepl(pattern="ENSG", x=rownames(deseqResults)[1]))
    stopifnot(all(c("padj", "log2FoldChange") %in% colnames(deseqResults)))
    stopifnot(order %in% c("zero", "first", "minSimple", "minSteiner"))
    stopifnot(hubMeasure %in% c("betweenness", "degree", "hubscore"))
    stopifnot(
        "'ppiData' must have columns 'ensemblGeneA', 'ensemblGeneB'"=all(
            c("ensemblGeneA", "ensemblGeneB") %in% colnames(ppiData)
        )
    )

    df <- as_tibble(rownames_to_column(deseqResults, "gene"))

    if (filterInput) {
        df <- filter(
            df,
            padj < pCutoff,
            abs(log2FoldChange) > log2(fcCutoff)
        )
    }

    message(
        "ppiBuildNetwork will use ", nrow(df),
        " genes for network construction..."
    )

    ## Check for and remove any duplicate IDs, warning the user when this occurs
    dfClean <- distinct(df, gene, .keep_all=TRUE)
    geneVector <- unique(dfClean[["gene"]])

    lostIds <- df[["gene"]][duplicated(df[["gene"]])]

    if (length(geneVector) < nrow(df)) {
        numDups <- nrow(df) - length(geneVector)

        message(
            "INFO: Found ", numDups,
            "duplicate IDs in the input column, which have been removed:"
        )

        if (numDups <= 10) {
            message(str_wrap(
                paste(lostIds, collapse=", "),
                indent=2,
                exdent=2
            ))
        } else {
            message(str_wrap(
                paste0(paste(lostIds[seq_len(10)], collapse=", "), "..."),
                indent=2,
                exdent=2
            ))
        }
    }

    ppiDataEnsembl <- select(ppiData, starts_with("ensembl"))

    if (order == "zero") {
        edgeTable <- ppiDataEnsembl %>% filter(
            ensemblGeneA %in% geneVector & ensemblGeneB %in% geneVector
        )
    } else {
        edgeTable <- ppiDataEnsembl %>% filter(
            ensemblGeneA %in% geneVector | ensemblGeneB %in% geneVector
        )
    }

    networkInit <- edgeTable %>%
        as_tbl_graph(directed=FALSE) %>%
        ppiRemoveSubnetworks() %>%
        as_tbl_graph() %>%
        mutate(
            degree=centrality_degree(),
            betweenness=centrality_betweenness(),
            seed=(name %in% geneVector)
        ) %>%
        select(-comp)

    ## Perform node filtering/trimming for minimum order networks, and
    ## recalculate degree and betweenness
    if (order == "minSimple") {
        message("Performing 'simple' minimum network trimming...")

        networkOut1 <- networkInit %>%
            filter(!(degree == 1 & !seed), !(betweenness == 0 & !seed)) %>%
            mutate(
                degree=centrality_degree(),
                betweenness=centrality_betweenness()
            )

    } else if (order == "minSteiner") {
        message("Performing 'Steiner' minimum network trimming...")

        terminals <- networkInit %>%
            activate(nodes) %>%
            pull(name) %>%
            intersect(geneVector)

        networkOut1 <- steinertree(
            type="SP",
            terminals=terminals,
            graph=networkInit,
            color=FALSE
        ) %>%
            .[[1]] %>%
            as_tbl_graph(directed=FALSE) %>%
            select(-color) %>%
            mutate(
                degree=centrality_degree(),
                betweenness=centrality_betweenness()
            )

    } else {
        networkOut1 <- networkInit
    }

    networkOut2 <-
        if (hubMeasure == "betweenness") {
            networkOut1 %>% mutate(hubScoreBtw=betweenness)
        } else if (hubMeasure == "degree") {
            networkOut1 %>% mutate(hubScoreDeg=degree)
        } else if (hubMeasure == "hubscore") {
            networkOut1 %>% mutate(hubScoreHub=centrality_hub())
        }

    if (nrow(as_tibble(networkOut2)) > 2000) {
        message(
            "Your network contains more than 2000 nodes, and will likely be ",
            "difficult to interpret when plotted."
        )
    }

    networkFinal <- networkOut2 %>%
        left_join(
            select(mappingFile, "name"=ensemblGeneId, hgncSymbol),
            by="name",
            multiple="all"
        ) %>%
        left_join(dfClean, by=c("name"="gene"), multiple="all")

    attr(networkFinal, "order") <- order
    return(networkFinal)
}
