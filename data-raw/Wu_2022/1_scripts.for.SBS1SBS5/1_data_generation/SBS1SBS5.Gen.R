# SBS1SBS5.Gen.R
# Generate 20 SBS1-SBS5 spectra datasets with varying correlation and mutation count ratio.
require(SynSigGen)

# Set working directory to the folder which contains results of
# computational approaches on SBS1-SBS5-correlated data sets
# before running this script.
#
# PATH <- paste0("<path_to_results_on_SBS1-SBS5-correlated_datasets>")
#
# setwd(PATH)
SynSigGen::CreateSBS1SBS5CorrelatedSyntheticData(
  top.level.dir = "./0.Input_datasets",
  regress.dir = NULL,
  overwrite = FALSE,
  add.info = TRUE,
  unlink = FALSE)
