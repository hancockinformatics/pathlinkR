.onAttach <- function(...) {

  packageStartupMessage(
    paste0(
      "Thanks for using pathnet v", utils::packageVersion("pathnet"), "! ",
      "If you encounter any bugs or problems, please submit an issue at the ",
      "Github page: https://github.com/hancockinformatics/pathnet/issues"
    ) %>% stringr::str_wrap(width = getOption("width"))
  )
}
