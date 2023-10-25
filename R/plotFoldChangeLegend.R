#' INTERNAL Construct heatmap legend
#'
#' @param .matFC Matrix of fold change values
#' @param .log2FoldChange Boolean denoting if values will be in log2
#' @param .cellColours Colours for fold change values
#'
#' @return A list containing heatmap legend parameters and colour function
#'
#' @importFrom circlize colorRamp2
#'
#' @description Helper function to handle heatmap legends without clutteing up
#'   the main function.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
.plotFoldChangeLegend <- function(.matFC, .log2FoldChange, .cellColours) {

    parameters <- list(
        title=ifelse(.log2FoldChange, "Log2 fold\nchange", "Fold change")
    )

    limit <- ceiling(max(abs(.matFC)))

    if (.log2FoldChange) {
        both <- c(-limit, -limit / 2, 0, limit / 2, limit)
        parameters <- append(parameters, list(at=both, labels=both))

    } else {
        if ((limit %% 2) != 0) limit <- limit + 1

        parameters <- append(parameters, list(
            at=c(-limit, -limit / 2, 0, limit / 2, limit),
            labels=c(
                -2 ^ limit,
                -2 ^ (limit / 2),
                1,
                2 ^ (limit / 2),
                2 ^ limit
            )
        ))
    }

    if (min(.matFC) >= 0) {
        myColFun <- colorRamp2(
            breaks=c(0, max(.matFC)),
            colors=c(.cellColours[2], .cellColours[3])
        )

        parameters[["at"]] <- parameters[["at"]][parameters[["at"]] >= 0]
        parameters[["labels"]] <-
            parameters[["labels"]][parameters[["labels"]] >= 0]

    } else if (max(.matFC) <= 0) {
        myColFun <- colorRamp2(
            breaks=c(min(.matFC), 0),
            colors=c(.cellColours[1], .cellColours[2])
        )

        parameters[["at"]] <- parameters[["at"]][parameters[["at"]] <= 0]
        parameters[["labels"]] <-
            parameters[["labels"]][parameters[["labels"]] <= 0]

    } else {
        myColFun <- colorRamp2(
            breaks=c(min(.matFC), 0, max(.matFC)),
            colors=.cellColours
        )
    }

    return(list(
        parameters,
        myColFun
    ))
}
