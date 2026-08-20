test_that("fold change heatmaps are working with some customization", {
    data("exampleDESeqResultsHS")

    expect_no_error(
        plotFoldChange(
            inputList=exampleDESeqResultsHS,
            species="human",
            pathName="Interferon alpha/beta signaling",
            hideNonsigFC=FALSE,
            invert=TRUE,
            clusterColumns=TRUE,
            returnData=TRUE
        )
    )
})

test_that("argument 'returnData' is working.", {
    data("exampleDESeqResultsHS")

    expect_type(
        plotFoldChange(
            inputList=exampleDESeqResultsHS,
            species="human",
            pathName="Interferon alpha/beta signaling",
            hideNonsigFC=FALSE,
            invert=TRUE,
            clusterColumns=TRUE,
            returnData=TRUE
        ),
        "list"
    )
})

test_that("fold change heatmaps are working with custom font sizes", {
    data("exampleDESeqResultsHS")

    expect_no_error(
        plotFoldChange(
            inputList=exampleDESeqResultsHS,
            species="human",
            pathName="PD-1 signaling",
            hideNonsigFC=FALSE,
            clusterColumns=TRUE,
            returnData=TRUE
        )
    )
})

test_that("fold change heatmaps are working with some customization", {
    data("exampleDESeqResultsMM")

    expect_no_error(
        plotFoldChange(
            inputList=exampleDESeqResultsMM,
            species="mouse",
            pathName="Interferon alpha/beta signaling",
            hideNonsigFC=FALSE,
            invert=TRUE,
            clusterColumns=TRUE,
            returnData=TRUE
        )
    )
})

test_that("argument 'returnData' is working.", {
    data("exampleDESeqResultsMM")

    expect_type(
        plotFoldChange(
            inputList=exampleDESeqResultsMM,
            species="mouse",
            pathName="Interferon alpha/beta signaling",
            hideNonsigFC=FALSE,
            invert=TRUE,
            clusterColumns=TRUE,
            returnData=TRUE
        ),
        "list"
    )
})

test_that("fold change heatmaps are working with custom font sizes", {
    data("exampleDESeqResultsMM")

    expect_no_error(
        plotFoldChange(
            inputList=exampleDESeqResultsMM,
            species="mouse",
            pathName="PD-1 signaling",
            hideNonsigFC=FALSE,
            clusterColumns=TRUE,
            returnData=TRUE
        )
    )
})
