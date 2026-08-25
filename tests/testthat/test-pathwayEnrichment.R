test_that("Sigora enrichment works as expected", {
    data("exampleDESeqResultsHS")
    data("reaH", package="sigora")

    suppressMessages(
        testResultSigora <- pathwayEnrichment(
            inputList=exampleDESeqResultsHS[1],
            species="human",
            analysis="sigora",
            gpsRepo=reaH,
            verbose=FALSE
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
    data("exampleDESeqResultsHS")
    data("kegH", package="sigora")

    suppressMessages(
        testResultSigoraKEGG <- pathwayEnrichment(
            inputList=exampleDESeqResultsHS[1],
            species="human",
            analysis="sigora",
            gpsRepo=kegH,
            gpsLevel=2,
            verbose=FALSE
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
    data("exampleDESeqResultsHS")

    expect_no_error(
        suppressMessages(
            testResultReactomepa <- pathwayEnrichment(
                inputList=exampleDESeqResultsHS,
                species="human",
                analysis="reactomepa",
                verbose=FALSE
            )
        )
    )
})

test_that("Hallmark enrichment works as expected", {
    data("exampleDESeqResultsHS")

    suppressMessages(
        testResultHallmark <- pathwayEnrichment(
            inputList=exampleDESeqResultsHS,
            species="human",
            analysis="hallmark",
            split=FALSE,
            verbose=FALSE
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
    data("exampleDESeqResultsHS")

    exampleDESeqResultsHS[[1]]$padj <- exampleDESeqResultsHS[[1]]$padj + 1

    testResultNoGenes <- pathwayEnrichment(
        inputList=exampleDESeqResultsHS,
        species="human",
        analysis="reactomepa",
        verbose=FALSE
    )

    expect_length(unique(testResultNoGenes$comparison), 1)
})

test_that("Sigora enrichment works as expected", {
    data("exampleDESeqResultsMM")
    data("reaM", package="sigora")

    suppressMessages(
        testResultSigora <- pathwayEnrichment(
            inputList=exampleDESeqResultsMM[1],
            species="mouse",
            analysis="sigora",
            gpsRepo=reaM,
            verbose=FALSE
        )
    )

    expect_equal(dim(testResultSigora), c(12, 12))

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
    data("exampleDESeqResultsMM")
    data("kegM", package="sigora")

    suppressMessages(
        testResultSigoraKEGG <- pathwayEnrichment(
            inputList=exampleDESeqResultsMM[1],
            species="mouse",
            analysis="sigora",
            gpsRepo=kegM,
            gpsLevel=2,
            verbose=FALSE
        )
    )

    expect_equal(dim(testResultSigoraKEGG), c(10, 12))

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
    data("exampleDESeqResultsMM")

    expect_no_error(
        suppressMessages(
            testResultReactomepa <- pathwayEnrichment(
                inputList=exampleDESeqResultsMM,
                species="mouse",
                analysis="reactomepa",
                verbose=FALSE
            )
        )
    )
})

test_that("a no-gene scenario works properly", {
    data("exampleDESeqResultsMM")

    exampleDESeqResultsMM[[1]]$padj <- exampleDESeqResultsMM[[1]]$padj + 1

    testResultNoGenes <- pathwayEnrichment(
        inputList=exampleDESeqResultsMM,
        species="mouse",
        analysis="reactomepa",
        verbose=FALSE
    )

    expect_length(unique(testResultNoGenes$comparison), 1)
})
