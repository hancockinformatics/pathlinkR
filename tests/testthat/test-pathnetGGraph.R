startingPathways <- createFoundation(
    mat = pathwayDistancesJaccard,
    maxDistance = 0.8
)

exPathnet <- createPathnet(
    sigoraResult = filter(
        sigoraExamples,
        comparison == "COVID Pos Over Time"
    ),
    foundation = startingPathways,
    trim = TRUE,
    trimOrder = 1
)

set.seed(1)

exPathnetGGraph <- pathnetGGraph(
    exPathnet,
    labelProp = 0.1,
    nodeLabelSize = 4,
    nodeLabelOverlaps = 8,
    segColour = "red"
)

test_that("we get the right plot", {
    vdiffr::expect_doppelganger("pathnetGGraphExample", exPathnetGGraph)
})
