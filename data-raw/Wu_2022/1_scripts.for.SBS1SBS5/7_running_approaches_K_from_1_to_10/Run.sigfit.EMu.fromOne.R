
# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)

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



## Run approach sigfit.EMu
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){

    cat("\n===========================================\n")
    cat(paste0("Running sigfit on data set ",datasetName," using seed ",seedInUse,"...\n"))
    cat("\n===========================================\n")

    Runsigfit(
      input.catalog = paste0(datasetName, "/sp.sp/ground.truth.syn.catalog.csv"),
      out.dir = paste0(datasetName,
                       "/sp.sp/ExtrAttrFromOne/sigfit.EMu.results/seed.",
                       seedInUse),
      model = "emu",
      CPU.cores = 10,
      seedNumber = seedInUse,
      K.range = c(1, 10),
      overwrite = TRUE)
  }
}

