test_that("string truncation works", {
  expect_equal(.truncNeatly("This is a test string", 18), "This is a test...")
})

test_that("NA in means NA out", {
    expect_equal(.truncNeatly(x=NA), NA_character_)
})
