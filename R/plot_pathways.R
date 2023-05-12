#' plot_pathways
#'
#' @param enriched_results Tibble of results from the function `enrich_pathway`
#' @param columns Number of columns to split the pathways across, especially if
#'   there are many pathways. Can specify up to 3 columns, with a default of 1.
#' @param specific_top_pathways Only plot pathways from a specific vector of
#'   "top_pathways". Defaults to "any" which includes all pathway results.
#' @param specific_pathways Only plot specific pathways. Defaults to "any".
#' @param name_width How many characters to show for pathway name before
#'   truncating? Defaults to 35.
#' @param name_rows How much to rows to wrap across for the pathway name?
#'   Defaults to 1.
#' @param x_angle Angle of x axis labels, set to "angled" (45 degrees),
#'   "horizontal" (0 degrees), or "vertical" (90 degrees).
#' @param max_pval P values below "10 ^ -max_pval" will be set to that value
#' @param intercepts Add vertical lines to separate different groupings, is a
#'   vector of intercepts (e.g. c(1.5, 2.5)). Defaults to "NA".
#' @param include_gene_ratio Boolean; should be plotted? If so, then it is
#'   attributed to the size of the triangles.
#' @param size Size of points if not scaling to gene ratio. Defaults to 5.
#' @param legend_multiply Size of the legend, e.g. increase if there are a lot
#'   of pathways which makes the legend small and unreadable by comparison.
#'   Defaults to 1, i.e. no increase in legend size.
#' @param show_num_genes Boolean, defaults to FALSE. Show the number of genes
#'   for each comparison as brackets under the comparison's name.
#'
#' @return A ggplot object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import stringr
#' @importFrom ggforce facet_col
#' @importFrom ggpubr rremove ggarrange
#'
#' @seealso <https://github.com/hancockinformatics/pathnet>
#'
plot_pathways <- function(
    enriched_results,
    columns = 1,
    specific_top_pathways = "any",
    specific_pathways = "any",
    name_width = 35,
    name_rows = 1,
    x_angle = "angled",
    max_pval = 50,
    intercepts = NA,
    include_gene_ratio = FALSE,
    size = 5,
    legend_multiply = 1,
    show_num_genes = FALSE
) {

  # TODO 1. Change top pathway names to shorten some of them

  # Convert to -log10 p value and arbitrarily set to a max -log10 p value of 50
  # (i.e. adj pval = 10^-50), which some enrichment results surpass, especially
  # for SIGORA.
  enriched_results <- enriched_results %>% mutate(
    logmax = case_when(
      -log10(p_value_adjusted) > max_pval ~ max_pval,
      -log10(p_value_adjusted) <= max_pval ~ -log10(p_value_adjusted)
    )
  )

  # Order the directionality of results
  enriched_results$direction <- factor(
    enriched_results$direction,
    levels = c("Up", "Down", "All")
  )

  # Add in the number of genes for each comparison if indicated
  if (show_num_genes) {
    enriched_results <- enriched_results %>% mutate(
      comparison = paste0(comparison, "\n(", total_genes, ")")
    )
    enriched_results$comparison <- factor(
      enriched_results$comparison,
      levels = unique(enriched_results$comparison)
    )
  }

  # If did not specify to only plot specific pathways, otherwise filter them
  if (specific_top_pathways[1] == "any") {
    specific_top_pathways <- unique(enriched_results$top_pathways)
  }
  if (specific_pathways[1] == "any") {
    specific_pathways <- unique(enriched_results$pathway_description)
  }

  enriched_results_graph <- enriched_results %>% filter(
    top_pathways %in% specific_top_pathways,
    pathway_description %in% specific_pathways
  )

  # enriched_results_graph$description <- factor(enriched_results_graph$description, levels = rev(specific_pathways))
  # enriched_results_graph <- enriched_results %>% filter(group %in% group_name)
  # # Organize the groups (e.g. D1, D7, etc.)
  # enriched_results_graph$group <- factor(enriched_results_graph$group, levels = group_name)
  # # Organize the id/comparisons (the different column names, e.g. sepsis vs healthy) in the appropriate order
  # enriched_results_graph$id <- factor(enriched_results_graph$id, levels = unique(enriched_results_graph$id))

  # In certain cases, a pathway may be enriched by both up- and down-regulated
  # genes. Find duplicated pathways, and only show the one that is more enriched
  # (lower p-value)
  duplicates <- enriched_results_graph %>%
    group_by(pathway_description, comparison) %>%
    summarise(counts = n()) %>% filter(counts > 1) %>%
    mutate(unique_id = paste0(pathway_description, comparison))

  enriched_results_graph <- enriched_results_graph %>%
    mutate(unique_id = paste0(pathway_description, comparison))

  enriched_results_clean <- enriched_results_graph %>%
    filter(!unique_id %in% duplicates$unique_id)

  enriched_results_dupes <- enriched_results_clean[0, ]

  # This chooses the top enriched pathway of duplicates and also adds the other
  # pathway into enriched_results_dupes.
  if (nrow(duplicates) > 0) {
    message(
      "\nWARNING: The following pathways were enriched in both directions for ",
      "the given comparisons. These are indicated with an asterisk over the ",
      "triangle, which is only shown for the lower p value result."
    )
    print(as_tibble(duplicates[, c(2, 1)]))

    for (i in 1:nrow(duplicates)) {
      row <- duplicates[i, ]

      choices <- enriched_results_graph %>%
        filter(pathway_description == row$pathway_description & comparison == row$comparison) %>%
        arrange(p_value_adjusted) # Choose the one that has lowest p value

      # Add the lower p value to the clean dataframe
      enriched_results_clean <- rbind(enriched_results_clean, choices[1, ])

      # keep the other enrichment to the dupes dataframe
      enriched_results_dupes <- rbind(enriched_results_dupes, choices[2, ])
    }
  }

  # Now organize pathways into multiple columns. Maximum is 3 columns to graph,
  # if inputted larger, will be set to 3.
  if (columns > 3) {
    message("Maximum is three columns to graph. Plotting three columns.")
    columns <- 3
  }

  ## How many pathways per top pathway? Add 1 to account for the extra space the
  ## facet takes up.
  num_pathways <- enriched_results_clean %>%
    dplyr::select(top_pathways, pathway_description) %>%
    unique() %>%
    group_by(top_pathways) %>%
    summarise(pathways = n() + 1)

  # Arrange ascending
  num_pathways <- num_pathways %>% arrange(pathways)

  # The n largest top pathways are chosen to start off the columns
  column_splitting <- num_pathways %>% tail(columns)

  # For cases when you specify specific pathways, make sure that there are more
  # pathways than columns...
  if (nrow(column_splitting) < columns) {
    columns <- nrow(column_splitting)

  } else {
    for (i in (columns + 1):nrow(num_pathways)) {
      # The "n + 1" largest top pathway to add
      add <- num_pathways %>% tail(i) %>% .[1,]

      column_smallest <- column_splitting %>%
        arrange(pathways) %>%
        head(1) %>%
        .$top_pathways # The current column that has the least number of pathways, add to this

      # Now, add the name of the top pathway to one of the columns and increase
      # the number of pathways
      column_splitting <- column_splitting %>% mutate(
        pathways = case_when(
          top_pathways == column_smallest ~ pathways + add$pathways,
          TRUE ~ pathways
        ),
        top_pathways = case_when(
          top_pathways == column_smallest ~ paste0(top_pathways, "," ,add$top_pathways),
          TRUE ~ top_pathways
        )
      )
    }
  }

  # Now create a list of top pathways for each column
  column_list <- column_splitting$top_pathways %>%
    as.list() %>%
    str_split(",")

  # Plot pathways
  plotlist <- list()
  name_trunc = name_width * name_rows - 5 # set where the pathway name should be truncated

  # Can be set to angled (45 degrees), "horizontal" (0 degrees), or "vertical"
  # (90 degrees)
  if(x_angle == "angled") {
    angle <- 45
    hjust <- 1
    vjust <- 1
  } else if (x_angle == "flat") {
    angle <- 0
    hjust <- 0.5
    vjust <- 1
  } else if (x_angle == "vertical") {
    angle <- 90
    hjust <- 0.5
    vjust <- 0.5
  }


  for (n in 1:length(column_list)) {
    plot <- ggplot(
      enriched_results_clean %>% filter(top_pathways %in% column_list[n][[1]]),
      aes(
        x = comparison,
        y = pathway_description,
        fill = logmax,
        shape = direction
      )
    ) +
      {if (include_gene_ratio) geom_point(aes(size = gene_ratio))} +
      {if (!include_gene_ratio) geom_point(size = size)} +
      geom_point(
        data = enriched_results_dupes %>% filter(top_pathways %in% column_list[n][[1]]),
        aes(x = comparison, y = pathway_description),
        shape = 8,
        size = 2,
        colour = "white",
        show.legend = FALSE
      ) +
      # Wrap and truncate pathway names if necessary
      scale_y_discrete(
        labels = function(x) str_wrap(trunc_neatly(x, name_trunc), width = name_width),
        position = "right"
      ) +
      # Keeps comparisons even if they don"t enrich for any pathways
      scale_x_discrete(drop = FALSE) +
      facet_col(
        facets = ~top_pathways,
        scales = "free_y",
        space = "free"
      ) +
      theme_bw() +
      theme(
        strip.text.x = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 12 * legend_multiply),
        legend.title = element_text(size = 13 * legend_multiply),
        axis.text.y = element_text(face = "bold", colour = "black", size = 12),
        axis.text.x = element_text(
          face = "bold",
          colour = "black",
          size = 12,
          angle = angle,
          hjust = hjust,
          vjust = vjust
        )
      ) +
      rremove("xlab") +
      rremove("ylab") +
      scale_shape_manual(
        values = c("Down" = 25 , "Up" = 24, "All" = 21),
        name = "Regulation",
        na.value = NA
      ) +
      scale_fill_continuous(
        name = expression(P[adjusted]),
        limits = c(0, 50),
        breaks = c(10, 20, 30, 40, 50),
        labels = c(
          expression(10 ^ -10),
          expression(10 ^ -20),
          expression(10 ^ -30),
          expression(10 ^ -40),
          expression(10 ^ -50)
        ),
        low = "blue",
        high = "red",
        na.value = NA
      ) +
      guides(
        shape = guide_legend(override.aes = list(size = 5 * legend_multiply)),
        size  = guide_legend(override.aes = list(shape = 24, fill = "black"))
      ) +
      # Can also add lines to separate different groups
      {if (!is.na(intercepts[1])) geom_vline(xintercept = intercepts)}

    plotlist <- append(plotlist, list(plot))
  }

  if (columns > 1) {
    plot <- ggarrange(
      plotlist = plotlist,
      ncol = columns,
      common.legend = TRUE,
      legend = "right",
      align = "v"
    )
    return(plot)
  } else {
    return(plot)
  }
}
