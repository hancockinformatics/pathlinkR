# Load packages -----------------------------------------------------------

devtools::load_all()


# Run enrichment ----------------------------------------------------------

sigoraExamples <- pathwayEnrichment(
    inputList=exampleDESeqResults,
    species="human",
    filterInput=TRUE,
    split=TRUE
)

glimpse(sigoraExamples)


# Save the data -----------------------------------------------------------

usethis::use_data(sigoraExamples, overwrite=TRUE)
