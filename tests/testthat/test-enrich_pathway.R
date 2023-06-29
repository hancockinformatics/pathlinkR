test_result_reactomepa <- enrich_pathway(
    input_list = deseq_example_list,
    analysis = "reactomepa"
)

test_result_hallmark <- enrich_pathway(
    input_list = deseq_example_list,
    analysis = "hallmark",
    split = FALSE
)

test_result_sigora <- sigora_examples

expected_colnames <- c(
    "comparison",
    "direction",
    "pathway_id",
    "pathway_description",
    "p_value",
    "p_value_adjusted",
    "genes",
    "num_candidate_genes",
    "num_bg_genes",
    "gene_ratio",
    "total_genes",
    "top_pathways"
)


test_that("We get the right number of dimensions from each method", {
    expect_equal(dim(test_result_sigora), c(66, 12))
    expect_equal(dim(test_result_reactomepa), c(121, 12))
    expect_equal(dim(test_result_hallmark), c(12, 12))
})

test_that("We have the right columns for each method", {
    expect_setequal(
        colnames(test_result_sigora),
        expected_colnames
    )
    expect_setequal(
        colnames(test_result_reactomepa),
        expected_colnames
    )
    expect_setequal(
        colnames(test_result_hallmark),
        expected_colnames
    )
})
