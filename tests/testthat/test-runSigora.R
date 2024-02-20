test_that("runSigora works on a basic level", {
    data("exampleDESeqResults")

    testGenes <- rownames(exampleDESeqResults[[1]])[seq(500)]
    expect_no_error(.runSigora(testGenes, gpsRepo="default"))
})
