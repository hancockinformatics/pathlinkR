test_starting_pathways <- create_foundation(
    mat = pathway_distances_jaccard,
    max_distance = 0.8
)

test_ex_pathway_network_input <- sigora_examples %>%
    filter(comparison == "COVID Pos Over Time")

test_my_pathway_network <- create_pathnet(
    sigora_result = sigora_examples,
    foundation = test_starting_pathways,
    trim = TRUE,
    trim_order = 1
)

test_that("there are the right number of nodes", {
    expect_length(test_my_pathway_network, 168)
})
