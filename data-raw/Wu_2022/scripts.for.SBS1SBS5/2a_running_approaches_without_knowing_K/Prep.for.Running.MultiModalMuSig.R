
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

## Create MultimodalMuSig-formatted input-catalog under
## <topLevelFolder4Run>/MultimodalMuSig.LDA.results/<datasetName> and
## <topLevelFolder4Run>/MultimodalMuSig.MMCTM.results/<datasetName>
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    CreateMultiModalMuSigOutput(
      catalog = paste0(topLevelFolder4Data,"/",datasetName,"/ground.truth.syn.catalog.csv"),
      out.dir = paste0(topLevelFolder4Run,"/MultiModalMuSig.LDA.results/",
	            datasetName,"/seed.",seedInUse),
      overwrite = T)
    CreateMultiModalMuSigOutput(
      catalog = paste0(topLevelFolder4Data,"/",datasetName,"/ground.truth.syn.catalog.csv"),
      out.dir = paste0(topLevelFolder4Run,"/MultiModalMuSig.MMCTM.results/",
	            datasetName,"/seed.",seedInUse),
      overwrite = T)
  }
}
