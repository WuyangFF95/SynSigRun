# Set working directory to the folder which contains results of
# computational approaches on SBS1-SBS5-correlated data sets
# before running this script.
#
# PATH <- paste0("<path_to_results_on_SBS1-SBS5-correlated_datasets>")
#
# setwd(PATH)
topLevelFolder4Data <- "./0.Input_datasets"
topLevelFolder4Run <- "./2b.Full_output_K_as_2"


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

## Specify 20 seeds used in software running
seedsInUse <- c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753)



## Run SomaticSignatures.NMF
## NMF will have an error if K.range includes 1!
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){

    out.dir <- paste0(topLevelFolder4Run,"/SomaticSignatures.NMF.results/",datasetName,"/seed.",seedInUse)
    if(file.exists(paste0(out.dir,"/inferred.exposures.csv"))) next

    cat("\n===========================================\n")
    cat(paste0("Running SomaticSignatures.NMF on data set ",datasetName," using seed ",seedInUse,"...\n"))
    cat("\n===========================================\n")

    RunSomaticSignatures(
      input.catalog = paste0(topLevelFolder4Data,"/",datasetName,"/ground.truth.syn.catalog.csv"),
      out.dir = out.dir,
      CPU.cores = 10,
      seedNumber = seedInUse,
      K.exact = 2,
      nrun.extract = 1,
      overwrite = TRUE)
  }
}
