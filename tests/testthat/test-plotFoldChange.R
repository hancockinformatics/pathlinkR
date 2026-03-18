test_that("fold change heatmaps are working with some customization", {
    data("exampleDESeqResults")

    expect_no_error(
        plotFoldChange(
            inputList=exampleDESeqResults,
            pathName="Interferon alpha/beta signaling",
            hideNonsigFC=FALSE,
            invert=TRUE,
            clusterColumns=TRUE
        )
    )
})

test_that("argument 'returnData' is working.", {
    data("exampleDESeqResults")

    expect_type(
        plotFoldChange(
            inputList = exampleDESeqResults,
            pathName = "Interferon alpha/beta signaling",
            hideNonsigFC = FALSE,
            invert = TRUE,
            clusterColumns = TRUE,
            returnData = TRUE
        ),
        "list"
    )
})

test_that("fold change heatmaps are working with custom font sizes", {
    data("exampleDESeqResults")

    expect_no_error(
        plotFoldChange(
            inputList=exampleDESeqResults,
            pathName="PD-1 signaling",
            hideNonsigFC=FALSE,
            clusterColumns=TRUE
        )
    )
})
