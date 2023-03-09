# pathnet
**pathnet** is an R package to create PPI-style networks from pathways, using
overlap between their constituent genes to define interactions/edges. This will
allow for visualization of similar pathways based on their underlying genes,
instead of a category/hierarchy term. One can then overlay the results of
pathway enrichment analyses onto this "foundational" network, mapping
dysregulation and significance onto nodes (pathways) in the network.


## Workflow
The proposed workflow currently looks like this:

- Start with the provided `sigora_database` object, which contains all human 
  Reactome pathways and their constituent genes
- Use the included `sigora_database` as input to `unnamed_function` to create a
  pairwise distance matrix, using any distance metric supported by
  `vegan::vegdist()`
  - Alternatively, we provide a pairwise Jaccard matrix that can be used as-is
- With `create_foundation()`, create a data frame of pathway nodes/edges, where
  pathways with a distance below a desired threshold are considered "connected"
- Turn this node/edge data frame into a network object that can be plotted and 
  explored using `ggraph` functions
  - We will make our own wrappers here to simplify the plotting process
- Add information from enriched pathways generated with `sigora`, such as 
  p-values, to include when visualizing the pathway networks
