# SBS1SBS5.Gen.R
# Generate 20 SBS1-SBS5 spectra datasets with varying correlation and mutation count ratio.
require(SynSigGen)

# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)
SynSigGen::CreateSBS1SBS5CorrelatedSyntheticData(
  top.level.dir = "../",
  regress.dir = NULL,
  overwrite = FALSE,
  add.info = TRUE,
  unlink = FALSE)
