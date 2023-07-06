testResultReactomepa <- enrichPathway(
    inputList = deseqExampleList,
    analysis = "reactomepa"
)

testResultHallmark <- enrichPathway(
    inputList = deseqExampleList,
    analysis = "hallmark",
    split = FALSE
)

test_result_sigora <- sigoraExamples

expectedColnames <- c(
    "comparison",
    "direction",
    "pathwayId",
    "pathwayName",
    "pValue",
    "pValueAdjusted",
    "genes",
    "numCandidateGenes",
    "numBgGenes",
    "geneRatio",
    "totalGenes",
    "topPathways"
)


test_that("we get the right number of dimensions from each method", {
    expect_equal(dim(testResultSigora), c(66, 12))
    expect_equal(dim(testResultReactomepa), c(121, 12))
    expect_equal(dim(testResultHallmark), c(12, 12))
})

test_that("we have the right columns for each method", {
    expect_setequal(
        colnames(testResultSigora),
        expectedColnames
    )
    expect_setequal(
        colnames(testResultReactomepa),
        expectedColnames
    )
    expect_setequal(
        colnames(testResultHallmark),
        expectedColnames
    )
})
