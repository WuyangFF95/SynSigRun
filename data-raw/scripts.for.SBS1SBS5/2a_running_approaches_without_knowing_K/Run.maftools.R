
# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)

## Load required packages
library(ICAMS)
library(SynSigRun)
library(NMF)



## Specify slopes and Rsqs for the datasets
slopes <- c(0.1,0.5,1,2,10)
Rsqs <- c(0.1,0.2,0.3,0.6)
datasetNames <- c()
for(slope in slopes)
  for(Rsq in Rsqs)
    datasetNames <- c(datasetNames, paste0("S.",slope,".Rsq.",Rsq))

## Specify 1 seed used in software running
## Note: maftools has been fixed to use seed 123456!
seedsInUse <- c(123456)



## Run Extraction and attribution packages
## sigproextractor (Python package) and MultiModalMuSig (Julia package)
## needs to be run with external script.
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){
    Runmaftools(input.catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
              out.dir = paste0(datasetName,"/sp.sp/ExtrAttr/maftools.results/seed.",seedInUse),
              seedNumber = seedInUse,
              K.range = c(1,10),
              overwrite = T)
  }
}
