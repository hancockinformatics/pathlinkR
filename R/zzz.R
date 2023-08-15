.onAttach <- function(...) {

    packageStartupMessage(
        paste0(
            "Thanks for using pathlinkR v", utils::packageVersion("pathlinkR"),
            "! If you encounter any bugs or problems, please submit an issue ",
            "at the Github page: https://github.com/hancockinformatics/",
            "pathlinkR/issues"
        ) %>% stringr::str_wrap(width=getOption("width"))
    )
}
