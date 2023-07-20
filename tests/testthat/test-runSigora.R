testGenes <- rownames(deseqExampleList[[1]])[seq(500)]

testGenesSmall <- testGenes[seq(10)]

test_that("runSigora works on a basic level", {
  expect_no_error(.runSigora(testGenes))
})
