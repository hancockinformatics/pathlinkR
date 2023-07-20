test_that("breaks are of the correct class", {
  expect_type(.eruptionBreaks(c(4, 5)), "environment")
})
