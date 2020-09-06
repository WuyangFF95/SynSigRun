
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

## Specify 20 seeds used in software running
seedsInUse <- c(123456)

## Enable memory sharing in NMF
install.extras("NMF")


## Run MutationalPatterns which requires package NMF.
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){
    RunMutationalPatterns(input.catalog = paste0(datasetName, "/sp.sp/ground.truth.syn.catalog.csv"),
                          out.dir = paste0(datasetName, "/sp.sp/ExtrAttrExact/MutationalPatterns.results/seed.", seedInUse),
                          CPU.cores = 10,
                          K.exact = 2,
                          overwrite = TRUE)
  }
}

