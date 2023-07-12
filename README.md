# pathnet

![GitHub R package version (branch)](https://img.shields.io/github/r-package/v/hancockinformatics/pathnet/master?label=pathnet%40master)

<!-- badges: start -->
[![R-CMD-check](https://github.com/hancockinformatics/pathnet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hancockinformatics/pathnet/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**pathnet** is an R package designed to facilitate analysis of RNA-Seq results.
Specifically, our aim with pathnet was to provide a number of tools which 
take a list of DE genes and perform different analyses on them, aiding with the
interpretation of results.

## Installation
```r
devtools::install_github("https://github.com/hancockinformatics/pathnet")
```

## Workflow & functions

The functions provided in `pathnet` can be grouped into a number of different 
general approaches: 

- Visualization of differential expression results with volcano plots and heatmaps
- Protein-Protein Interaction (PPI) network creation and visualization
- Pathway enrichment, and plots to summarize findings, including drawing results
as networks of connected pathways

See the vignette for an example workflow and more details on each of these
methods and the functions they use.
    
## Contributors
Andy An & Travis Blimkie at the CMDR REW Hancock Lab.

## Versioning
This package makes use of [SemVer](https://semver.org/).

## License
This project uses the GNU General Public License v3.0, available
[here](https://github.com/hancockinformatics/SeptiSearch/blob/master/LICENSE).

<br>

[<img src="man/figures/hancock-lab-logo.svg" height="40px">](http://cmdr.ubc.ca/bobh/)
