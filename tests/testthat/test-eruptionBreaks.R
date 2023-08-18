test_that("breaks are of the correct class", {
    expect_type(.eruptionBreaks(c(-5, 5)), "environment")
    expect_type(.eruptionBreaks(c(-6.5, 6.5)), "environment")
    expect_type(.eruptionBreaks(c(-9, 9)), "environment")
    expect_type(.eruptionBreaks(c(-13, 13)), "environment")
    expect_type(.eruptionBreaks(c(-16, 16)), "environment")
})

test_that("the correct message is returned", {
    expect_message(
        .eruptionBreaks(c(-20, 20)),
        regexp=paste0(
            "Something may be wrong with your DESeq model"
        )
    )
})
