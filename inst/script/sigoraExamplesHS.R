# Load packages -----------------------------------------------------------

devtools::load_all()


# Run enrichment ----------------------------------------------------------

sigoraExamplesHS <- pathwayEnrichment(
    inputList=exampleDESeqResultsHS,
    species="human",
    filterInput=TRUE,
    split=TRUE
)


# Save the data -----------------------------------------------------------

usethis::use_data(sigoraExamplesHS, overwrite=TRUE)
