test_that("ReactomePA enrichment works as expected", {
    testResultReactomepa <- enrichPathway(
        inputList = deseqExampleList,
        analysis = "reactomepa"
    )

    expect_equal(dim(testResultReactomepa), c(121, 12))

    expect_setequal(
        colnames(testResultReactomepa),
        c(
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
    )
})

test_that("we have the right columns for each method", {
    testResultHallmark <- enrichPathway(
        inputList = deseqExampleList,
        analysis = "hallmark",
        split = FALSE
    )

    expect_equal(dim(testResultHallmark), c(12, 12))

    expect_setequal(
        colnames(testResultHallmark),
        c(
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
    )
})
