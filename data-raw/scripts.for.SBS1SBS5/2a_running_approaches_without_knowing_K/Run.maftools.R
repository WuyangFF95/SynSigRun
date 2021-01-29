
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
## Note: The seed for maftools is hard-coded as 123456.
seedsInUse <- c(123456)

## Enable memory sharing in NMF
## Applicable for R 3.6
## Not applicable for R 4.0+
#install.extras("NMF")


## NMF will have an error if K.exact equals to 1
## Run approach maftools
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){

    out.dir <- paste0(topLevelFolder4Run,"/maftools.results/",datasetName,"/seed.",seedInUse)
    if(file.exists(paste0(out.dir,"/inferred.exposures.csv"))) next

    cat("\n===========================================\n")
    cat(paste0("Running maftools on data set ",datasetName," using seed ",seedInUse,"...\n"))
    cat("\n===========================================\n")


    Runmaftools(input.catalog = paste0(topLevelFolder4Data,"/",datasetName,"/ground.truth.syn.catalog.csv"),
              out.dir = out.dir,
              CPU.cores = 10,
              K.range = c(2,10),
              overwrite = T)
  }
}
