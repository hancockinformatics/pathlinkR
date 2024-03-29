% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{mappingFile}
\alias{mappingFile}
\title{Table of human gene ID mappings}
\format{
A data frame (tibble) with 43,993 rows and 3 columns
\describe{
\item{ensemblGeneId}{Ensembl IDs}
\item{hgncSymbol}{HGNC symbols}
\item{entrezGeneId}{NCBI Entrez IDs}
}
}
\source{
See \url{https://bioconductor.org/packages/biomaRt/} and
\url{https://bioconductor.org/packages/AnnotationDbi/} for information on
each of the utilized packages and functions.
}
\usage{
data(mappingFile)
}
\value{
An object of class "tbl", "tbl.df", "data.frame"
}
\description{
A data frame to aid in mapping human gene IDs between different
formats, inclusing Ensembl IDs, HGNC symbols, and Entrez IDs. Mapping
information was sourced using \code{biomaRt} and \code{AnnotationDbi}.
}
\keyword{datasets}
