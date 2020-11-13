
# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)

## Load required packages
library(ICAMS)
library(SynSigEval)

## Specify slopes and Rsqs for the datasets
slopes <- c(0.1,0.5,1,2,10)
Rsqs <- c(0.1,0.2,0.3,0.6)
datasetNames <- c()
for(slope in slopes)
  for(Rsq in Rsqs)
    datasetNames <- c(datasetNames, paste0("S.",slope,".Rsq.",Rsq))


for(datasetName in datasetNames){
  for(nrun in 1:20){
    out.dir <- paste0(datasetName,"/sp.sp/ExtrAttrExact/EMu.results/run.",nrun)
    dir.create(out.dir, recursive = TRUE)
    SynSigRun::CreateEMuOutput(
      catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
      out.dir = out.dir,
      overwrite = F)
  }
}
