#' List of example results from DESeq2
#'
#' @format A list of two data frames, each with 6 columns (and rownames),
#'   containing example results from DESeq2
#' \describe{
#'   \item{baseMean}{A combined score for the gene}
#'   \item{log2FoldChange}{Fold change value for the gene}
#'   \item{lfcSE}{Standard error for the fold change value}
#'   \item{stat}{The statistic value}
#'   \item{pvalue}{The nominal p value for the gene}
#'   \item{padj}{The adjusted p value for the gene}
#' }
"deseqExampleList"

#' Colour assignments for grouped pathways
#'
#' @format A length 8 named vector of hex colour values
"groupedPathwayColours"

#' A data frame (tibble) of ID mappings for different human ID types, from
#' the Sigora package.
#'
#' @format A data frame with 47037 rows and 3 columns
#' \describe{
#'   \item{Ensembl.Gene.ID}{Ensembl gene ID}
#'   \item{EntrezGene.ID}{Entrez/NCBI gene ID}
#'   \item{Symbol}{HGNC symbol}
#' }
"idmap"

#' A data frame containing human PPI data from InnateDB, from the entry
#' "All Experimentally Validated Interactions (updated weekly)" at
#' <https://innatedb.com/redirect.do?go=downloadImported>.
#' A few important steps have been taken to filter the data, namely the removal
#' of duplicate interactions, and removing interactions that have the same
#' components but are swapped between A and B.
#'
#' @format A data frame with 152259 rows and 4 columns:
#' \describe{
#'   \item{ensemblGeneA}{Ensembl gene ID for the first gene/protein in the
#'     interaction}
#'   \item{hgncSymbolA}{HGNC symbol for the first gene/protein in the
#'     interaction}
#'   \item{ensemblGeneB}{Ensembl gene ID for the second gene/protein in the
#'     interaction}
#'   \item{hgncSymbolB}{HGNC symbol for the second gene/protein in the
#'     interaction}
#' }
"innateDbExp"

#' Table of human gene ID mappings
#'
#' @format A data frame (tibble) with 43,993 rows and 3 columns
#' \describe{
#'   \item{ensemblGeneId}{Ensembl IDs}
#'   \item{hgncSymbol}{HGNC symbols}
#'   \item{entrezGeneId}{NCBI Entrez IDs}
#' }
"mappingFile"

#' Table of Hallmark gene sets and their genes
#'
#' @format A data frame (tibble) with 8,209 rows and 2 columns
#' \describe{
#'   \item{pathwayId}{Name of the Hallmark Gene Set}
#'   \item{ensemblGeneId}{Ensembl IDs}
#' }
"mSigDbTermToGene"

#' Manually-curated list of Reactome and Hallmark pathways and their top
#' pathways and grouped pathways
#'
#' @format A data frame (tibble) with 2685 rows and 5 columns
#' \describe{
#'   \item{pathwayId}{Reactome or Hallmark pathway identifier}
#'   \item{pathwayName}{Pathway name}
#'   \item{topLevelPathway}{Top hierarchy pathway term, shortened in some cases}
#'   \item{groupedPathway}{Top grouped pathway, 8 for Reactome}
#'   \item{topLevelOriginal}{Original top pathway name}
#' }
"pathwayCategories"

#' Matrix of pairwise Jaccard indices for all human Reactome pathways
#'
#' @format A matrix with 1001 rows and columns
#' \describe{
#'   \item{Rows}{Reactome pathway identifiers}
#'   \item{Columns}{Reactome pathway identifiers}
#' }
"pathwayDistancesJaccard"

#' Table of all Reactome pathways and genes
#'
#' @format A data frame (tibble) with 123574 rows and 3 columns
#' \describe{
#'   \item{pathwayId}{Reactome pathway ID}
#'   \item{entrezGeneId}{Entrez gene ID}
#'   \item{pathwayName}{Name of the Reactome pathway}
#' }
"reactomeDatabase"

#' The default Sigora Gene Pair Signature (GPS) object
#'
#' @format A list of length nine with the following top-level elements
#' \describe{
#'   \item{origRepo}{The original pathway data}
#'   \item{L1}{Level 1 pathways}
#'   \item{L2}{Level 2 pathways}
#'   \item{L3}{Level 3 pathways}
#'   \item{L4}{Level 4 pathways}
#'   \item{L5}{Level 5 pathways}
#'   \item{repoName}{Name of the repository used to build the GPS object}
#'   \item{pathwaydescriptions}{Pathway names}
#'   \item{call}{Function call to create the GPS object}
#' }
"reaH"

#' Table of all Sigora pathways and their constituent genes
#'
#' @format A data frame (tibble) with 60775 rows and 4 columns
#' \describe{
#'   \item{pathwayId}{Reactome pathway identifier}
#'   \item{pathwayName}{Reactome pathway description}
#'   \item{ensemblGeneId}{Ensembl gene identifier}
#'   \item{hgncSymbol}{HGNC gene symbol}
#' }
"sigoraDatabase"

#' Example Sigora output from running `pathwayEnrichment` on
#'   "deseqExampleList"
#'
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
"sigoraExamples"
