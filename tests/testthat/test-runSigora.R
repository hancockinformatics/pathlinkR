test_that("runSigora works on a basic level", {
    testGenes <- rownames(deseqExampleList[[1]])[seq(500)]

    testGenesSmall <- testGenes[seq(10)]


    expect_no_error(.runSigora(testGenes, gpsRepo="default"))
})
