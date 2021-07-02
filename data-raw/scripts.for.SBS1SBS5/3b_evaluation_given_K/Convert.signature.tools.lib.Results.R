
## Run 2b_running_approaches_given_K/Run.signature.tools.lib.exact.R
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


## After running EMu, convert EMu-formatted tsv files
## to SynSigEval/ICAMS-formatted csv files.
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    ## Grep the names of EMu-formatted files.
    ## When the K is different, the output file would be different.
    ## Extracted signature file name:
    ## _{K}_ml_spectra.txt
    ## Attributed exposure file name:
    ## _{K}_assigned.txt
    resultDir <- paste0(topLevelFolder4Run,"/signature.tools.lib.results/",
      datasetName,"/seed.",seedInUse)


    files <- list.files(resultDir)

    signatureFile <- files[grep(pattern = "Sigs_plot_test_ns2_nboots20_MC.tsv",x = files)]


    ## Convert signatures
    signatures <- utils::read.table(
      paste0(resultDir,"/",signatureFile),
      sep = "\t")
    signatures <- ICAMS::as.catalog(
      signatures,
      region = "unknown",
      catalog.type = "counts.signature")
    ## extracted signatures need to be normalized.
    for(sigName in colnames(signatures)){
      signatures[,sigName] <- signatures[,sigName] / sum(signatures[,sigName])
    }
    ICAMS::WriteCatalog(
      signatures,
      paste0(resultDir,"/extracted.signatures.csv"))

  }
}
