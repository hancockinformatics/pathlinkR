test_that("a volcano plot works with some of customizations", {
    data("exampleDESeqResults", "sigoraDatabase")

    interferonGenes <- sigoraDatabase %>%
        filter(pathwayName == "Interferon Signaling") %>%
        pull(ensemblGeneId)

    expect_no_error(
        eruption(
            rnaseqResult=exampleDESeqResults[[1]],
            xaxis=c(-4, 4),
            yaxis=c(0, 8),
            highlightGenes=interferonGenes
        )
    )
})

test_that("eruption returns a ggplot object", {

  data("exampleDESeqResults")

  p <- eruption(
    exampleDESeqResults[[1]],
    columnFC = "logFC",
    columnP = "padj"
  )

  expect_s3_class(p, "ggplot")
})

test_that("xaxis and yaxis must have length two", {

  data("exampleDESeqResults")

  expect_error(
    eruption(
      exampleDESeqResults[[1]],
      columnFC = "logFC",
      columnP = "padj",
      xaxis = c(-2, 2, 3)
    ),
    "length-two"
  )

  expect_error(
    eruption(
      exampleDESeqResults[[1]],
      columnFC = "logFC",
      columnP = "padj",
      yaxis = 1
    ),
    "length-two"
  )
})

test_that("plot title is applied", {

  data("exampleDESeqResults")

  p <- eruption(
    exampleDESeqResults[[1]],
    columnFC = "logFC",
    columnP = "padj",
    title = "My Volcano"
  )

  expect_equal(
    p$labels$title,
    "My Volcano"
  )
})

test_that("cutoff lines are placed correctly", {

  data("exampleDESeqResults")

  p <- eruption(
    exampleDESeqResults[[1]],
    columnFC = "logFC",
    columnP = "padj",
    fcCutoff = 2,
    pCutoff = 0.05
  )

  built <- ggplot_build(p)

  hline <- built$data[[4]]
  vline <- built$data[[5]]

  expect_equal(
    unique(hline$yintercept),
    -log10(0.05)
  )

  expect_equal(
    sort(unique(vline$xintercept)),
    c(-1, 1)
  )
})

test_that("manual axis limits are respected", {

  data("exampleDESeqResults")

  p <- eruption(
    exampleDESeqResults[[1]],
    columnFC = "logFC",
    columnP = "padj",
    xaxis = c(-1, 1),
    yaxis = c(0, 3)
  )

  built <- ggplot_build(p)

  expect_equal(
    built$layout$panel_params[[1]]$x.range,
    c(-1.1, 1.1)
  )

  expect_equal(
    built$layout$panel_params[[1]]$y.range,
    c(-0.15, 3.15)
  )
})
