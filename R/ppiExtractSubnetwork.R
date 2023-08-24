#' Extract a subnetwork based on pathway genes
#'
#' @param network Input network object; output from `ppiBuildNetwork`
#' @param genes List of Ensembl gene IDs to use as the starting point to extract
#'   a subnetwork from the initial network. You must provide either the
#'   `genes` or `pathwayEnrichmentResult` argument.
#' @param pathwayEnrichmentResult Pathway enrichment result, output from
#'   `ppiEnrichNetwork`. You must provide either `genes` or
#'   `pathwayEnrichmentResult` argument.
#' @param pathwayToExtract Name of the pathway determining what genes (nodes)
#'   are pulled from the input network. Must be present in
#'   `pathwayEnrichmentResults`
#'
#' @return A `tidygraph` object for plotting or further analysis
#'
#' @export
#'
#' @import dplyr
#' @importFrom igraph as.igraph V induced.subgraph decompose.graph simplify
#'   delete.vertices get.shortest.paths
#'
#' @details Uses functions from the igraph package to extract a minimally
#'   connected subnetwork from the starting network, using either a list of
#'   Ensembl genes or genes from an enriched pathway as the basis. To see what
#'   genes were pulled out for the pathway, check `attr(x, "starters")`.
#'
#' @references Code for network module (subnetwork) extraction was based off of
#' that used in "jboktor/NetworkAnalystR" on Github.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
#' @examples
#' exNetwork <- ppiBuildNetwork(
#'     deseqResults=deseqExampleList[[1]],
#'     filterInput=TRUE,
#'     order="zero"
#' )
#'
#' exPathways <- ppiEnrichNetwork(
#'     network=exNetwork,
#'     analysis="sigora"
#' )
#'
#' ppiExtractSubnetwork(
#'     network=exNetwork,
#'     pathwayEnrichmentResult=exPathways,
#'     pathwayToExtract="Interferon alpha/beta signaling"
#' )
#'
ppiExtractSubnetwork <- function(
        network,
        genes=NULL,
        pathwayEnrichmentResult=NULL,
        pathwayToExtract
) {

    stopifnot(
        "You must specify either 'genes' or 'pathwayEnrichmentResult' to
        provide genes to extract from the initial network" = {
            !all(is.null(genes), is.null(pathwayEnrichmentResult))
        }
    )

    if (!is.null(genes)) {
        stopifnot(
            "Argument 'genes' must be a list of Ensembl gene IDs" = {
                grepl(x = genes[[1]], pattern = "^ENSG[0-9]+$")
            }
        )
    }

    if (!is.null(pathwayEnrichmentResult)) {
        stopifnot(
            "Argument 'pathwayEnrichmentResult' must contain the columns
            'pathwayName' and 'genes'" = {
                all(c("pathwayName", "genes") %in%
                        colnames(pathwayEnrichmentResult))
            }
        )

        stopifnot(
            "Provided 'pathwayToExtract' must be present in
            'pathwayEnrichmentResult" = {
                pathwayToExtract %in% pathwayEnrichmentResult[["pathwayName"]]
            }
        )

        stopifnot(
            "The 'genes' column must contain HGNC symbols separated with a ';'" = {
                grepl(
                    x=pathwayEnrichmentResult[["genes"]][1],
                    pattern="\\w+;"
                )
            }
        )
    }


    if (!is.null(genes)) {
        message("Using provided list of Ensembl genes...", appendLF = FALSE)
        genesToExtract <- genes

    } else if (!is.null(pathwayEnrichmentResult)) {
        message("Pulling genes for given pathway...", appendLF = FALSE)

        pathwayGenesHGNC <- pathwayEnrichmentResult %>%
            filter(pathwayName == pathwayToExtract) %>%
            pull(genes) %>%
            strsplit(., split=";") %>%
            unlist()

        genesToExtract <- mappingFile %>%
            filter(hgncSymbol %in% pathwayGenesHGNC) %>%
            pull(ensemblGeneId) %>%
            unique()
    }

    message("found ", length(genesToExtract), " genes.")


    geneNodeIds <- as_tibble(network) %>%
        mutate(rn=row_number()) %>%
        filter(name %in% genesToExtract) %>%
        pull(rn)

    # Get subgraphs, which will only contain the specified nodes
    message("Calculating subgraphs from specified nodes...")
    allSubgraphs <- induced.subgraph(
        graph=as.igraph(network),
        vids=geneNodeIds
    )

    # Decompose each subgraph
    allComponents <- decompose.graph(graph=allSubgraphs, min.vertices=1)

    # If all the specified nodes form a single, connected network, pull that...
    if (length(allComponents) == 1) {
        message("All nodes form a single network...")
        moduleNetwork <- allComponents[[1]]

        # ...or we need to minimally connect each subgraph we've identified
    } else {
        message("Determining shortest paths between nodes...")

        moduleShortestPaths <- list()

        for (i in seq(length(geneNodeIds))) {
            moduleShortestPaths[[i]] <- get.shortest.paths(
                as.igraph(network),
                geneNodeIds[i],
                geneNodeIds[-(1:i)]
            )$vpath
        }

        moduleNodes <- unique(unlist(moduleShortestPaths))
        NodesToRemove <- igraph::V(network)$name[-moduleNodes]
        moduleNetwork <- simplify(delete.vertices(network, NodesToRemove))
    }

    moduleNetworkTidygraph <- as_tbl_graph(moduleNetwork)

    message(
        "Done, new subnetwork contains ",
        nrow(as_tibble(moduleNetworkTidygraph)),
        " nodes.\n"
    )

    attr(moduleNetworkTidygraph, "starters") <- genesToExtract
    return(moduleNetworkTidygraph)
}
