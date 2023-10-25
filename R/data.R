#' List of example results from DESeq2
#'
#' @usage data(exampleDESeqResults)
#' @return An object of class "list"
#' @format A list of two "DESeqResults" objects, each with 5000 rows and 6
#'   columns:
#' \describe{
#'   \item{baseMean}{A combined score for the gene}
#'   \item{log2FoldChange}{Fold change value for the gene}
#'   \item{lfcSE}{Standard error for the fold change value}
#'   \item{stat}{The statistic value}
#'   \item{pvalue}{The nominal p value for the gene}
#'   \item{padj}{The adjusted p value for the gene}
#' }
#'
#' @source For details on DESeq2 and its data structures/methods, please see
#' \url{https://bioconductor.org/packages/DESeq2/}
#'
"exampleDESeqResults"

#' Colour assignments for grouped pathways
#'
#' @usage data(groupedPathwayColours)
#' @return An object of class "character"
#' @format A length 8 named vector of hex colour values
#'
"groupedPathwayColours"

#' Table of Hallmark gene sets and their genes
#'
#' @usage data(hallmarkDatabase)
#' @return An object of class "tbl", "tbl.df", "data.frame"
#' @format A data frame (tibble) with 8,209 rows and 2 columns
#' \describe{
#'   \item{pathwayId}{Name of the Hallmark Gene Set}
#'   \item{ensemblGeneId}{Ensembl gene IDs}
#' }
#'
#' @source For more information on the MSigDB Hallmark gene sets, please see
#' \url{https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp}
#'
"hallmarkDatabase"

#' InnateDB PPI data
#'
#' @description A data frame containing human PPI data from InnateDB, from the
#'   entry "All Experimentally Validated Interactions (updated weekly)" at
#'   <https://innatedb.com/redirect.do?go=downloadImported>. A few important
#'   steps have been taken to filter the data, namely the removal of duplicate
#'   interactions, and removing interactions that have the same components but
#'   are swapped between A and B.
#'
#' @usage data(innateDbPPI)
#' @return An object of class "tbl", "tbl.df", "data.frame"
#' @format A data frame (tibble) with 152,256 rows and 2 columns:
#' \describe{
#'   \item{ensemblGeneA}{Ensembl gene ID for the first gene/protein in the
#'     interaction}
#'   \item{ensemblGeneB}{Ensembl gene ID for the second gene/protein in the
#'     interaction}
#' }
#'
#' @source For more details on the data sourced from InnateDB, please see
#' their website: \url{https://www.innatedb.com}
#'
"innateDbPPI"

#' Table of human gene ID mappings
#'
#' @description A data frame to aid in mapping human gene IDs between different
#' formats, inclusing Ensembl IDs, HGNC symbols, and Entrez IDs. Mapping
#' information was sourced using \code{biomaRt} and \code{AnnotationDbi}.
#'
#' @usage data(mappingFile)
#' @return An object of class "tbl", "tbl.df", "data.frame"
#' @format A data frame (tibble) with 43,993 rows and 3 columns
#' \describe{
#'   \item{ensemblGeneId}{Ensembl IDs}
#'   \item{hgncSymbol}{HGNC symbols}
#'   \item{entrezGeneId}{NCBI Entrez IDs}
#' }
#'
#' @source See \url{https://bioconductor.org/packages/biomaRt/} and
#'   \url{https://bioconductor.org/packages/AnnotationDbi/} for information on
#'   each of the utilized packages and functions.
#'
"mappingFile"

#' Top-level pathway categories
#'
#' @description A data frame containing all Reactome pathways and Hallmark
#' terms, along with a manually-curated top-level category for each entry.
#'
#' @usage data(pathwayCategories)
#' @return An object of class "tbl", "tbl.df", "data.frame"
#' @format A data frame (tibble) with 2685 rows and 5 columns
#' \describe{
#'   \item{pathwayId}{Reactome or Hallmark pathway identifier}
#'   \item{pathwayName}{Pathway name}
#'   \item{topLevelPathway}{Top hierarchy pathway term, shortened in some cases}
#'   \item{groupedPathway}{Top grouped pathway, 8 for Reactome}
#'   \item{topLevelOriginal}{Original top pathway name}
#' }
#'
#' @source See \url{https://reactome.org/} and
#'   \url{https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp} for
#'   information on each of these databases.
#'
"pathwayCategories"

#' Table of all Reactome pathways and genes
#'
#' @usage data(reactomeDatabase)
#' @return An object of class "tbl", "tbl.df", "data.frame"
#' @format A data frame (tibble) with 123574 rows and 3 columns
#' \describe{
#'   \item{pathwayId}{Reactome pathway ID}
#'   \item{entrezGeneId}{Entrez gene ID}
#'   \item{pathwayName}{Name of the Reactome pathway}
#' }
#'
#' @source See \url{https://reactome.org/} for information on each of this
#'   patwhay resource.
#'
"reactomeDatabase"

#' Table of all Sigora pathways and their constituent genes
#'
#' @usage data(sigoraDatabase)
#' @return An object of class "tbl", "tbl.df", "data.frame"
#' @format A data frame (tibble) with 60775 rows and 4 columns
#' \describe{
#'   \item{pathwayId}{Reactome pathway identifier}
#'   \item{pathwayName}{Reactome pathway description}
#'   \item{ensemblGeneId}{Ensembl gene identifier}
#'   \item{hgncSymbol}{HGNC gene symbol}
#' }
#'
#' @source Please refer to the Sigora package for more details:
#'   \url{https://cran.r-project.org/package=sigora}
#'
"sigoraDatabase"

#' Sigora enrichment example
#'
#' @description Example Sigora output from running `pathwayEnrichment()` on
#'   "exampleDESeqResults"
#'
#' @usage data(sigoraExamples)
#' @return An object of class "tbl", "tbl.df", "data.frame"
#' @format A data frame (tibble) with 66 rows and 12 columns
#' \describe{
#'   \item{comparison}{Comparison from which results are derived; names of the
#'     input list}
#'   \item{direction}{Was the pathway enriched in up or down regulated genes}
#'   \item{pathwayId}{Reactome pathway identifier}
#'   \item{pathwayName}{Description of the pathway}
#'   \item{pValue}{Nominal p value for the enrichment}
#'   \item{pValueAdjusted}{p value adjusted for multiple testing}
#'   \item{genes}{Genes in the pathway/input}
#'   \item{numCandidateGenes}{Analyzed genes found in the pathway of interest}
#'   \item{numBgGenes}{All genes from the pathway database}
#'   \item{geneRatio}{Quotient of the number of candidate and background genes}
#'   \item{totalGenes}{Total number of input genes}
#'   \item{topLevelPathway}{Pathway category}
#' }
#'
#' @source Please refer to the Sigora package for more details on that method:
#'   \url{https://cran.r-project.org/package=sigora}
#'
"sigoraExamples"
