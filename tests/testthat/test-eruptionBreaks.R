test_that("breaks are of the correct class", {
    expect_type(.eruptionBreaks(c(-5, 5)), "environment")
})

test_that("the correct message is returned", {
    expect_message(
        .eruptionBreaks(c(-20, 20)),
        regexp=paste0(
            "Something may be wrong with your DESeq model"
        )
    )
})
