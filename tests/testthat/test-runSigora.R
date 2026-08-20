test_that("runSigora works on a basic level", {
    data("exampleDESeqResultsHS")

    testGenes <- rownames(exampleDESeqResultsHS[[1]])[seq(500)]
    testResult <- .runSigora(testGenes, gpsRepo="reaH", gpsLevel=4)

    expect_equal(dim(testResult), c(331, 8))
})

test_that("runSigora works on a basic level", {
    data("exampleDESeqResultsMM")

    testGenes <- rownames(exampleDESeqResultsMM[[1]])[seq(500)]
    testResult <- .runSigora(testGenes, gpsRepo="reaM", gpsLevel=4)

    expect_equal(dim(testResult), c(280, 8))
})

