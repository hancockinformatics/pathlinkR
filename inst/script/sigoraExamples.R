# Load packages -----------------------------------------------------------

devtools::load_all()


# Run enrichment ----------------------------------------------------------

sigoraExamples <- pathwayEnrichment(
    inputList=exampleDESeqResults,
    filterInput=TRUE,
    split=TRUE
)


# Save the data -----------------------------------------------------------

usethis::use_data(sigoraExamples, overwrite=TRUE)
