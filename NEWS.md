CHANGES IN VERSION 1.1.17
-------------------------

* Fix for default parameters and theme settings in ppiPlotNetwork
* Support for GSEA (via package fgsea), both as an option in pathwayPlots and pathwayEnrichment
* Global fontSize argument added to pathwayPlots


CHANGES IN VERSION 1.1.10
-------------------------

* Support for KEGG pathways when running pathwayEnrichment, either through Sigora or traditional over-representation


CHANGES IN VERSION 1.1.8
-------------------------

* pathwayPlots now supports the output from fgsea on Reactome pathways (109a173)


CHANGES IN VERSION 1.1.7
-------------------------

* Size of asterisks in pathwayPlots is set based on tringle size (1f9b751)


CHANGES IN VERSION 1.1.6
-------------------------

* Updated installation information in README
* "eruption" gained argument "labelCutoffs" to add labels to p value and fold change cutoff lines
* "ppiPlotNetwork" checks for argument "labelColumn" when "label=TRUE"
* Fixed some tests


CHANGES IN VERSION 0.99.380
----------------------------

Lots of updates, including:

* Updated documentation for data and functions
* Cleaner vignette
* Better support for Bioconductor classes
* Better handling of data loading


CHANGES IN VERSION 0.99.352
----------------------------

* Added value tag to data documentation


CHANGES IN VERSION 0.99.350
----------------------------

* Rewrite of pathwayPlots to be more efficient


CHANGES IN VERSION 0.99.332
----------------------------

* Changed methodology for importing functions
* Added two new functions for PPI networks: "ppiEnrichNetwork" and "ppiExtractSubnetwork"


CHANGES IN VERSION 0.99.317
----------------------------

* Renamed a number of functions to improve consistency, and make it clear which 
  ones are related and part of the same "workflow"


CHANGES IN VERSION 0.99.310
----------------------------

* Big update to the vignette
* Update to the README and added in the package hex logo
* Various function tweaks and minor changes


CHANGES IN VERSION 0.99.300
----------------------------

* Renamed package to pathlinkR


CHANGES IN VERSION 0.99.1
----------------------------

* Prepared for bioconductor submission
