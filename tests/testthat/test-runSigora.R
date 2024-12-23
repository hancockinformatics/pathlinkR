test_that("runSigora works on a basic level", {
    data("exampleDESeqResults")

    testGenes <- rownames(exampleDESeqResults[[1]])[seq(500)]
    testResult <- .runSigora(testGenes, gpsRepo="default")

    expect_equal(dim(testResult), c(331, 8))
})
