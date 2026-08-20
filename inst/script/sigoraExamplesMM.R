# Load packages -----------------------------------------------------------

devtools::load_all()


# Run enrichment ----------------------------------------------------------

data(reaM, package="sigora")

sigoraExamplesMM <- pathwayEnrichment(
    inputList=exampleDESeqResultsMM,
    species="mouse",
    filterInput=TRUE,
    split=TRUE,
    analysis="sigora",
    gpsRepo="reaM",
    gpsLevel=4,
    verbose=TRUE
)


# Save the data -----------------------------------------------------------

usethis::use_data(sigoraExamplesMM, overwrite=TRUE)
