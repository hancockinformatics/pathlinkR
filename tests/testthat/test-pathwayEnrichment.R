test_that("Sigora enrichment works as expected", {
    data("exampleDESeqResults")

    suppressMessages(
        testResultSigora <- pathwayEnrichment(
            inputList=exampleDESeqResults[1],
            analysis="sigora"
        )
    )

    expect_equal(dim(testResultSigora), c(30, 12))

    expect_setequal(
        colnames(testResultSigora),
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

test_that("Sigora enrichment works with KEGG", {
    data("exampleDESeqResults")

    suppressMessages(
        testResultSigoraKEGG <- pathwayEnrichment(
            inputList=exampleDESeqResults[1],
            analysis="sigora",
            gpsRepo="kegH"
        )
    )

    expect_equal(dim(testResultSigoraKEGG), c(30, 12))

    expect_setequal(
        colnames(testResultSigoraKEGG),
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
    data("exampleDESeqResults")

    suppressMessages(
        testResultReactomepa <- pathwayEnrichment(
            inputList=exampleDESeqResults,
            analysis="reactomepa"
        )
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
    data("exampleDESeqResults")

    suppressMessages(
        testResultHallmark <- pathwayEnrichment(
            inputList=exampleDESeqResults,
            analysis="hallmark",
            split=FALSE
        )
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
