# pathnet
**pathnet** is an R package to create PPI-style networks from Reactome pathways, 
using overlap between constituent genes to define interactions/edges. This will
allow for visualization of similar pathways based on their underlying genes,
instead of simple grouping by category/hierarchy term. One can then overlay the
results of pathway enrichment analyses onto this "foundational" network, mapping
dysregulation and significance onto nodes (pathways) in the network.

## Installation
```r
devtools::install_github("https://github.com/hancockinformatics/pathnet")
```

## Workflow

- Use the included `sigora_database` table (containing all human
  Reactome pathways and their genes) as input to `get_pathway_distances` to 
  build a pairwise distance matrix, using any distance metric supported by 
  `vegan::vegdist()`
  - Alternatively, we provide a pre-made distance matrix (Jaccard) that can be 
    used as-is
- With `create_foundation()`, build a data frame of pathway nodes/edges, where
  pathways with a distance below a desired threshold are considered "connected"
  - One can also simply select a top proportion of all pathways to consider as
    connected, which can more easily handle a variety of distance measures
- Turn this node/edge data frame into a network object with `create_pathnet`,
  with support for trimming non-enriched pathways
- Visualize this pathway network using the `pathnet_ggraph` function, which
  supports a large number of layout options and visual tweaks
  - By default, nodes are coloured by a manually-curated top-level pathway type,
    with filled nodes denoting enriched pathways and black nodes "interactor 
    pathways"
  - Alternatively, you can use `pathnet_visNetwork` to create an interactive 
    network using the `visNetwork` package
    
## Contributors
Andy An & Travis Blimkie at the CMDR REW Hancock Lab.

## Versioning
This package makes use of [SemVer](https://semver.org/).

## License
This project uses the GNU General Public License v3.0, available
[here](https://github.com/hancockinformatics/SeptiSearch/blob/master/LICENSE).

<br>

[<img src="man/figures/hancock-lab-logo.svg" height="40px">](http://cmdr.ubc.ca/bobh/)
