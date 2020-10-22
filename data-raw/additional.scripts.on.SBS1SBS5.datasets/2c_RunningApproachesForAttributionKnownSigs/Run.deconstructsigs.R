
# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)

## Load required packages
library(ICAMS)
library(SynSigRun)



## Specify slopes and Rsqs for the datasets
slopes <- c(0.1,0.5,1,2,10)
Rsqs <- c(0.1,0.2,0.3,0.6)
datasetNames <- c()
for(slope in slopes)
  for(Rsq in Rsqs)
    datasetNames <- c(datasetNames, paste0("S.",slope,".Rsq.",Rsq))

## Specify 20 seeds used in software running
seedsInUse <- c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753)



## Run Extraction and attribution packages
## sigproextractor (Python package) and MultiModalMuSig (Julia package)
## needs to be run with external script.
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){
    RundeconstructSigsAttributeOnly(
      input.catalog = paste0(datasetName, "/sp.sp/ground.truth.syn.catalog.csv"),
      gt.sigs.file = paste0(datasetName, "/sp.sp/ground.truth.syn.sigs.csv"),
      out.dir = paste0(datasetName, "/sp.sp/Attr/deconstructSigs.results/seed.", seedInUse),
      seedNumber = seedInUse,
      overwrite = T)
  }
}

