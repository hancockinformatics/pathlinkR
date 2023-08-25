#' INTERNAL Break long strings at spaces
#'
#' @param x Character to be truncated
#' @param l Desired maximum length for the output character
#'
#' @return Character vector
#'
#' @import dplyr
#'
#' @importFrom purrr map_chr
#' @importFrom stringr str_length str_sub str_replace
#'
#' @description Trims a character string to the desired length, without breaking
#'   in the middle of a word (i.e. chops at the nearest space). Appends an
#'   ellipsis at the end to indicate some text has been removed.
#'
#' @references None.
#'
#' @seealso <https://github.com/hancockinformatics/pathlinkR>
#'
.truncNeatly <- function(x, l=60) {
    map_chr(
        x,
        ~if (is.na(.x)) {
            return(NA_character_)
        } else if (str_length(.x) <= l) {
            return(.x)
        } else {
            shortened <- .x %>%
                as.character() %>%
                str_sub(., start=1, end=l) %>%
                str_replace(., pattern="\\s([^\\s]*)$", replacement="...")
            return(shortened)
        }
    )
}
