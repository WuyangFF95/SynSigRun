
## Run 2b_running_approaches_given_K/Run.MultiModalMuSig.exact.jl
## Before running this script.
##
## Run this script before running Summarize.extr.attr.K.given.R



# Set working directory to the folder which contains results of
# computational approaches on SBS1-SBS5-correlated data sets
# before running this script.
#
# PATH <- paste0("<path_to_results_on_SBS1-SBS5-correlated_datasets>")
# setwd(PATH)
topLevelFolder4Data <- "./0.Input_datasets"
topLevelFolder4Run <- "./2b.Full_output_K_as_2"


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


## After running, convert MultiModalMuSig-formatted tsv files
## to SynSigEval/ICAMS-formatted csv files.
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    for(approach in c("LDA","CTM")){
      ## Convert signatures
      signatures <- SynSigEval::MMCatalog2ICAMS(
        paste0(resultDir,"/extracted.signatures.tsv"),
        region = "unknown",
        catalog.type = "counts.signature")
      ## extracted signatures need to be normalized.
      for(sigName in colnames(signatures)){
        signatures[,sigName] <- signatures[,sigName] / sum(signatures[,sigName])
      }
      ICAMS::WriteCatalog(
        signatures,
        paste0(resultDir,"/extracted.signatures.csv"))

      ## Convert exposures
      rawExposure <- SynSigEval::ReadExposureMM(
        paste0(resultDir, "/inferred.exposures.tsv"))
      ## The sum of exposure of each spectrum needs to
      ## be normalized to the total number of mutations
      ## in each spectrum.
      spectra <- ICAMS::ReadCatalog(
        file = paste0(topLevelFolder4Data, "/", datasetName,
                      "/ground.truth.syn.catalog.csv"),
        catalog.type = "counts",
        strict = FALSE)
      exposureCounts <- rawExposure
      for(sample in colnames(exposureCounts)){
        exposureCounts[,sample] <- rawExposure[,sample] / sum(rawExposure[,sample]) * sum(spectra[,sample])
      }

      SynSigGen::WriteExposure(
        exposureCounts,
        paste0(resultDir, "/inferred.exposures.csv"))
    }
  }
}
