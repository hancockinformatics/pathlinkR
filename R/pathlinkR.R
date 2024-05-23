#' pathlinkR a package for analyzing RNA-Seq data
#'
#' The pathlinkR package is a suite of functions designed to facilitate the
#' analysis and visualization of RNA-Seq results. The main functions are:
#'
#' \itemize{
#' \item \code{\link{eruption}} - Create volcano plots from RNA-Seq results
#' \item \code{\link{plotFoldChange}} - Heatmaps to visualize and compare gene
#' expression across multiple conditions
#' \item \code{\link{pathwayEnrichment}} - Test DE genes for enriched Reactome
#' pathways or Hallmark terms, with different methods supported. Results can be
#' visualized with \code{\link{pathwayPlots}}
#' \item \code{\link{ppiBuildNetwork}} - Construct PPI networks from DE genes,
#' using interaction data from InnateDB. Networks can be plotted with
#' \code{\link{ppiPlotNetwork}}, tested for enriched pathways with
#' \code{\link{ppiEnrichNetwork}}, or subnetworks extracted using
#' \code{\link{ppiExtractSubnetwork}}
#' \item \code{\link{pathnetCreate}} - Turn pathway enrichment results into a
#' network of connected pathways, and create static plots with
#' \code{\link{pathnetGGraph}} or interactive plots with
#' \code{\link{pathnetVisNetwork}}
#' }
#'
#' For more details, please see the package vignette by entering
#' \code{vignette("pathlinkR")} into the console. Another document with more
#' examples is linked near the top of the included vignette.
#'
#' Any software-related questions can be posted on the Bioconductor Support
#' site: \url{https://support.bioconductor.org}
#'
#' The code is made publicly available on our Github page:
#' \url{https://github.com/hancockinformatics/pathlinkR}
#'
#' @author Travis Blimkie, Andy An
#'
#' @docType _PACKAGE
#' @name pathlinkR-package
#' @title pathlinkR
#' @aliases pathlinkR
#' @keywords package
NULL
