% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathwayEnrichment.R
\name{pathwayEnrichment}
\alias{pathwayEnrichment}
\title{Test significant DE genes for enriched pathways}
\usage{
pathwayEnrichment(
  inputList,
  columnFC = NA,
  columnP = NA,
  filterInput = TRUE,
  pCutoff = 0.05,
  fcCutoff = 1.5,
  split = TRUE,
  analysis = "sigora",
  filterResults = "default",
  gpsRepo = "reaH",
  gpsLevel = "default",
  geneUniverse = NULL,
  verbose = FALSE
)
}
\arguments{
\item{inputList}{A list, with each element containing RNA-Seq results as a
"DESeqResults", "TopTags", or "data.frame" object. Rownames of each table
must contain Ensembl Gene IDs. The list names are used as the comparison
name for each element (e.g. "COVID vs Healthy"). See Details for more
information on supported input types.}

\item{columnFC}{Character; Column to plot along the x-axis, typically log2
fold change values. Only required when \code{rnaseqResult} is a simple data
frame. Defaults to NA.}

\item{columnP}{Character; Column to plot along the y-axis, typically nominal
or adjusted p values. Only required when \code{rnaseqResult} is a simple data
frame. Defaults to NA.}

\item{filterInput}{When providing list of data frames containing the
unfiltered RNA-Seq results (i.e. not all genes are significant), set this
to \code{TRUE} to remove non-significant genes using the thresholds set by the
\code{pCutoff} and \code{fcCutoff}. When this argument is \code{FALSE}
its assumed your passing a pre-filtered data in \code{inputList}, and no
more filtering will be done.}

\item{pCutoff}{Adjusted p value cutoff when filtering. Defaults to < 0.05.}

\item{fcCutoff}{Minimum absolute fold change value when filtering. Defaults
to > 1.5}

\item{split}{Boolean (TRUE); Split into up- and down-regulated DE genes using
the fold change column, and do enrichment independently on each. Results
are combined at the end, with an added "direction" column.}

\item{analysis}{Method/database to use for enrichment analysis. The default
is "sigora", but can also be "reactome"/"reactomepa", "hallmark" or "kegg".}

\item{filterResults}{Should the output be filtered for significance? Use \code{1}
to return the unfiltered results, or any number less than 1 for a custom
p-value cutoff. If left as \code{default}, the significance cutoff for
\code{analysis="sigora"} is 0.001, or 0.05 for "reactome", "hallmark", and
"kegg".}

\item{gpsRepo}{Only applies to \code{analysis="sigora"}. Gene Pair Signature (GPS)
object for Sigora to use to test for enriched pathways. "reaH" (default)
will use the Reactome GPS object from \code{Sigora}; "kegH" will use the KEGG
GPS. One can also provide their own GPS object; see Sigora's documentation
for details.}

\item{gpsLevel}{Only applies to \code{analysis="sigora"}. If left as \code{default},
will be set to \code{4} for \code{gpsRepo="reaH"} or \code{2} for \code{gpeRepo="kegH"}. If
providing your own GPS object, can be set as desired; see Sigora's
documentation for details.}

\item{geneUniverse}{Only applies when \code{analysis} is "reactome"/"reactomepa",
"hallmark", or "kegg". The set of background genes to use when testing with
Reactome, Hallmark, or KEGG gene sets. For Reactome this must be a
character vector of Entrez genes. For Hallmark or KEGG, it must be Ensembl
IDs.}

\item{verbose}{Logical; If FALSE (the default), don't print info/progress
messages.}
}
\value{
A "data.frame" (tibble) of pathway enrichment results for all input
comparisons, with the following columns:
\item{comparison}{Source comparison from the names of \code{inputList}}
\item{direction}{Whether the pathway was enriched in all genes
(\code{split=FALSE}), or up- or down-regulated genes (\code{split=TRUE})}
\item{pathwayId}{Pathway identifier}
\item{pathwayName}{Pathway name}
\item{pValue}{Nominal p value for the pathway}
\item{pValueAdjusted}{p value, corrected for multiple testing}
\item{genes}{Candidate genes, which were DE for the comparison and also in
the pathway}
\item{numCandidateGenes}{Number of candidate genes}
\item{numBgGenes}{Number of background genes for the pathway}
\item{geneRatio}{Ratio of candidate and background genes}
\item{totalGenes}{Number of DE genes which were tested for enriched
pathways}
\item{topLevelPathway}{High level Reactome term which serves to group
similar pathways}
}
\description{
This function provides a simple and consistent interface to
three different pathway enrichment tools: Sigora and ReactomePA (which both
test for Reactome pathways), and MSigDB Hallmark gene set enrichment.
}
\details{
\code{inputList} must be a named list of RNA-Seq results, with each
element being of class "DESeqResults" from \code{DESeq2}, "TopTags" from
\code{edgeR}, or a simple data frame. For the first two cases, column names are
expected to be the standard defined by each class ("log2FoldChange" and
"padj" for "DESeqResults", and "logFC" and "FDR" for "TopTags"). Hence for
these two cases the arguments \code{columnFC} and \code{columnP}
can be left as \code{NA}.

In the last case (elements are "data.frame"), both \code{columnFC} and
\code{columnP} must be supplied when \code{filterInput=TRUE},
and \code{columnFC} must be given if \code{split=TRUE}.

Setting \code{analysis} to any of "reactome", "reactomepa", "hallmark", or
"kegg" will execute traditional over-representation analysis, the only
difference being the database used ("reactome" and "reactomepa" are treated
the same). Setting \code{analysis="sigora"} will use a gene pair-based approach,
which can be performed on either Reactome data when \code{gpsRepo="reaH"} or
KEGG data with \code{gpsRepo="kegH"}.
}
\examples{
data("exampleDESeqResults")

pathwayEnrichment(
    inputList=exampleDESeqResults[1],
    filterInput=TRUE,
    split=TRUE,
    analysis="hallmark",
    filterResults="default"
)

}
\references{
Sigora: \url{https://cran.r-project.org/package=sigora}
ReactomePA: \url{https://www.bioconductor.org/packages/ReactomePA/}
Reactome: \url{https://reactome.org/}
MSigDB/Hallmark: \url{https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp}
KEGG: \url{https://www.kegg.jp/}
}
\seealso{
\url{https://github.com/hancockinformatics/pathlinkR}
}
