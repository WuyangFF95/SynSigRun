
# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)


topLevelFolder4Data <- "../research_data/0.Input_datasets"
topLevelFolder4Run <- "../research_data/2a.Full_output_K_unspecified"


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

## Specify 20 seeds used in software running
seedsInUse <- c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753)

## Create helmsman-formatted input-catalog under
## <dataset.top.level.dir>/sp.sp/ExtrAttr/helmsman.results
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    out.dir <- paste0(topLevelFolder4Run,"/helmsman.NMF.results/",datasetName,"/seed.",seedInUse)
    dir.create(out.dir, recursive = TRUE)
    SynSigEval::CreatehelmsmanOutput(
      catalog = paste0(topLevelFolder4Data,"/",datasetName,"/ground.truth.syn.catalog.csv"),
      out.dir = out.dir,
      overwrite = T)
  }
}

