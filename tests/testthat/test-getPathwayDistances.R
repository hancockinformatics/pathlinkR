test_that("we have basic functionality", {
    data("sigoraDatabase")

    exGetPathwayDistances <- getPathwayDistances(
        pathwayData = dplyr::slice_head(
            dplyr::arrange(sigoraDatabase, pathwayId),
            prop = 0.05
        ),
        distMethod = "jaccard"
    )

    expect_equal(nrow(exGetPathwayDistances), ncol(exGetPathwayDistances))
})

test_that("an error is returned with wrong input type", {
    expect_error(getPathwayDistances(pathwayData = "wrong"))
})

test_that("a lack of requisite columns produces an error", {
    testtable <- tibble("a" = c(1, 2, 3), "b" = c("x", "y", "z"))
    expect_error(getPathwayDistances(pathwayData = testtable))
})
