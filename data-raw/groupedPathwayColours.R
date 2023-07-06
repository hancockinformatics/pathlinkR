# Use the colours from one of the RColorBrewer sets, plus one custom colour
# as the eighth from the set is grey, which I'd like to avoid using here.

groupedPathwayColours <- c(
  RColorBrewer::brewer.pal(n = 7, name = "Set2"),
  "#A65628"
)

names(groupedPathwayColours) <- c(
  "Cell Process",
  "Cell Replication",
  "Signaling",
  "Tissue Function",
  "Immune/Hemostasis",
  "Metabolism",
  "Gene Expression",
  "Disease"
)

usethis::use_data(groupedPathwayColours, overwrite = TRUE)
