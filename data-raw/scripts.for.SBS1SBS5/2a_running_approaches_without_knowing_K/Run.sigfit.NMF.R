
# Set working directory to the folder which contains results of
# computational approaches on SBS1-SBS5-correlated data sets
# before running this script.
#
# PATH <- paste0("<path_to_results_on_SBS1-SBS5-correlated_datasets>")
#
# setwd(PATH)
topLevelFolder4Data <- "./0.Input_datasets"
topLevelFolder4Run <- "./2a.Full_output_K_unspecified"


## Load required packages
library(ICAMS)
library(SynSigRun)
library(sigfit)
library(rstan)
library(rstantools)

## For faster computation, enable parallel computing in rstan.
##
## Allow precompiled rstan program to be written temporarily,
## so that parallel call may be faster
rstan::rstan_options(auto_write = TRUE)



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



## Run approaches sigfit.NMF
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){

    out.dir = paste0(topLevelFolder4Run,"/sigfit.NMF.results/",datasetName,"/seed.",seedInUse)
    if(file.exists(paste0(out.dir,"/inferred.exposures.csv"))) next

    cat("\n===========================================\n")
    cat(paste0("Running sigfit on data set ",datasetName," using seed ",seedInUse,"...\n"))
    cat("\n===========================================\n")

    Runsigfit(
      input.catalog = paste0(topLevelFolder4Data,"/",datasetName,"/ground.truth.syn.catalog.csv"),
      out.dir = out.dir,
      model = "nmf",
      CPU.cores = 10,
      seedNumber = seedInUse,
      K.range = c(2, 10),
      overwrite = TRUE)
  }
}

