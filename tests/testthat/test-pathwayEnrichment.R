test_that("Sigora enrichment works as expected", {
    data("exampleDESeqResults")

    suppressMessages(
        testResultSigora <- pathwayEnrichment(
            inputList=exampleDESeqResults[1],
            species="human",
            analysis="sigora",
            gpsRepo="reaH"
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
            species="human",
            analysis="sigora",
            gpsRepo="kegH",
            gpsLevel=2
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

    expect_no_error(
        suppressMessages(
            testResultReactomepa <- pathwayEnrichment(
                inputList=exampleDESeqResults,
                species="human",
                analysis="reactomepa"
            )
        )
    )
})

test_that("Hallmark enrichment works as expected", {
    data("exampleDESeqResults")

    suppressMessages(
        testResultHallmark <- pathwayEnrichment(
            inputList=exampleDESeqResults,
            species="human",
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

test_that("a no-gene scenario works properly", {
    data("exampleDESeqResults")

    exampleDESeqResults[[1]]$padj <- exampleDESeqResults[[1]]$padj + 1

    testResultNoGenes <- pathwayEnrichment(
        inputList=exampleDESeqResults,
        species="human",
        analysis="reactomepa",
        verbose=TRUE
    )

    expect_length(unique(testResultNoGenes$comparison), 1)
})
