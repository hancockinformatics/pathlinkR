#' INTERNAL Create manual breaks/labels for volcano plots
#'
#' @param x Length-two numeric vector to manually specify limits of the
#'   x-axis in log2 fold change; defaults to NA which lets ggplot2 determine the
#'   best values.
#'
#' @return ggplot scale object
#'
#' @import dplyr
#' @import ggplot2
#'
#' @description Internal function which is used to create even breaks for
#'   volcano plots produced by `eruption`.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
.eruptionBreaks <- function(x) {

    ## Up to log2FC=6
    if (between(max(abs(x)), 0, 6)) {
        scale_x_continuous(
            breaks=c(-6, -4, -2, 0, 2, 4, 6),
            labels=c(-32, -16, -4, 1, 4, 16, 32),
            limits=x
        )

    ## Up to log2FC=6.64
    } else if (between(max(abs(x)), 6.001, log2(100))) {
        scale_x_continuous(
            breaks=c(-log2(100), -log2(10), 0, log2(10), log2(100)),
            labels=c(
                expression(-10 ^ 2),
                expression(-10 ^ 1),
                1,
                expression(10 ^ 1),
                expression(10 ^ 2)
            ),
            limits=x
        )

    ## Up to log2FC=9.97
    } else if (between(max(abs(x)), log2(100), log2(1000))) {
        scale_x_continuous(
            breaks=c(
                -log2(1000),
                -log2(100),
                -log2(10),
                0,
                log2(10),
                log2(100),
                log2(1000)
            ),
            labels=c(
                expression(-10 ^ 3),
                expression(-10 ^ 2),
                expression(-10 ^ 1),
                1,
                expression(10 ^ 1),
                expression(10 ^ 2),
                expression(10 ^ 3)
            ),
            limits=x
        )

    ## Up to log2FC=13.29
    } else if (between(max(abs(x)), log2(1000), log2(10000))) {
        scale_x_continuous(
            breaks=c(
                -log2(10000),
                -log2(1000),
                -log2(100),
                -log2(10),
                0,
                log2(10),
                log2(100),
                log2(1000),
                log2(10000)
            ),
            labels=c(
                expression(-10 ^ 4),
                expression(-10 ^ 3),
                expression(-10 ^ 2),
                expression(-10 ^ 1),
                1,
                expression(10 ^ 1),
                expression(10 ^ 2),
                expression(10 ^ 3),
                expression(10 ^ 4)
            ),
            limits=x
        )

    ## Up to log2FC=16.61
    } else if (between(max(abs(x)), log2(10000), log2(100000))) {
        scale_x_continuous(
            breaks=c(
                -log2(10 ^ 5),
                -log2(10 ^ 3),
                -log2(10),
                0,
                log2(10),
                log2(10 ^ 3),
                log2(10 ^ 5)
            ),
            labels=c(
                expression(-10 ^ 5),
                expression(-10 ^ 3),
                expression(-10 ^ 1),
                1,
                expression(10 ^ 1),
                expression(10 ^ 3),
                expression(10 ^ 5)
            ),
            limits=x
        )

    } else if (max(abs(x)) >= log2(100000)) {
        message(
            "Something may be wrong with your DESeq model to have fold ",
            "changes > 10^5..."
        )
    }
}
