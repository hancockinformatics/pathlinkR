test_that("NA input returns NA_character_", {
  expect_identical(.truncNeatly(NA), NA_character_)
  expect_identical(.truncNeatly(NA_character_, l=10), NA_character_)
})

test_that("strings shorter than or equal to l are returned unchanged", {
  expect_identical(.truncNeatly("short", l=10), "short")
  expect_identical(.truncNeatly("", l=10), "")
  # exactly l characters long
  expect_identical(.truncNeatly("1234567890", l=10), "1234567890")
})

test_that("strings longer than l are truncated at the nearest space with ellipsis", {
  expect_identical(.truncNeatly("hello world", l=10), "hello...")
  expect_identical(
    .truncNeatly("The quick brown fox jumps", l=15),
    "The quick..."
  )
})

test_that("default l=60 is used when not specified", {
  short_string <- strrep("a", 60)
  long_string  <- paste(strrep("a", 55), "more words here")

  expect_identical(.truncNeatly(short_string), short_string)

  result <- .truncNeatly(long_string)
  expect_true(str_length(result) <= 60 + 3) # allow for "..." appended
  expect_true(grepl(x=result, pattern= "\\.\\.\\.$"))
})

test_that("no ellipsis is added if there's no space within the first l characters", {
  x <- "Supercalifragilisticexpialidocious extra words here"
  result <- .truncNeatly(x, l=10)

  # substring has no space in the first 10 characters, so regex won't match
  expect_identical(result, str_sub(x, 1, 10))
  expect_false(grepl(x=result, pattern= "\\.\\.\\.$"))
})

test_that("function is vectorized over a character vector", {
  x <- c(
    NA,
    "short",
    "hello world",
    "The quick brown fox jumps"
  )

  result <- .truncNeatly(x, l=10)

  expect_length(result, 4)
  expect_identical(result[1], NA_character_)
  expect_identical(result[2], "short")
  expect_identical(result[3], "hello...")
  expect_identical(result[4], "The quick...")
})

test_that("empty character vector returns empty character vector", {
  expect_identical(.truncNeatly(character(0)), character(0))
})

test_that("l=0 truncates everything without spaces to itself (no match possible)", {
  # str_sub with end=0 returns "", which has length <= 0? Actually str_length("") <= 0 is TRUE
  expect_identical(.truncNeatly("anything", l=0), "")
})

test_that("non-character input coercible to character is handled inside truncation branch", {
  # as.character() is applied inside the truncation branch, so numeric input longer
  # than l (as a string) should still be processed correctly.
  x <- 123456789012345
  result <- .truncNeatly(x, l=5)
  expect_identical(result, str_sub(as.character(x), 1, 5))
})
