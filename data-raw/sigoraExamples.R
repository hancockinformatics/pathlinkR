# Load packages -----------------------------------------------------------

devtools::load_all()


# Run enrichment ----------------------------------------------------------

sigoraExamples <- pathwayEnrichment(
    inputList=deseqExampleList,
    filterInput=TRUE,
    split=TRUE
)


# Save the data -----------------------------------------------------------

usethis::use_data(sigoraExamples, overwrite=TRUE)
