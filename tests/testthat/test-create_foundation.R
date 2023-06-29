test_starting_pathways <- create_foundation(
    mat = pathway_distances_jaccard,
    max_distance = 0.8
)

test_that("The foundation has the right dimensions", {
  expect_equal(dim(test_starting_pathways), c(6852, 5))
})

test_that("The foundation has the right columns", {
    expect_setequal(
        colnames(test_starting_pathways),
        c(
            "pathway_name_1",
            "pathway_name_2",
            "distance",
            "pathway_1",
            "pathway_2"
        )
    )
})
