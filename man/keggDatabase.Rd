% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{keggDatabase}
\alias{keggDatabase}
\title{Table of KEGG pathways and genes}
\format{
A data frame (tibble) with 32883 rows and 4 columns
\describe{
\item{pathwayId}{KEGG pathway ID}
\item{pathwayName}{Name of the Reactome pathway}
\item{ensemblGeneId}{Ensembl gene ID}
\item{hgncSymbol}{HGNC gene symbol}
}
}
\source{
See \url{https://kegg.jp} for more information.
}
\usage{
data(keggDatabase)
}
\value{
An object of class "tbl", "tbl.df", "data.frame"
}
\description{
Table of KEGG pathways and genes
}
\keyword{datasets}
