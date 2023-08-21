test_that("Sigora enrichment works as expected", {
    testResultReactomepa <- pathwayEnrichment(
        inputList=deseqExampleList,
        analysis="sigora"
    )

    expect_equal(dim(testResultReactomepa), c(66, 12))

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
            "topLevelPathway"
        )
    )
})

test_that("ReactomePA enrichment works as expected", {
    testResultReactomepa <- pathwayEnrichment(
        inputList=deseqExampleList,
        analysis="reactomepa"
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
            "topLevelPathway"
        )
    )
})

test_that("Hallmark enrichment works as expected", {
    testResultHallmark <- pathwayEnrichment(
        inputList=deseqExampleList,
        analysis="hallmark",
        split=FALSE
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
            "topLevelPathway"
        )
    )
})
