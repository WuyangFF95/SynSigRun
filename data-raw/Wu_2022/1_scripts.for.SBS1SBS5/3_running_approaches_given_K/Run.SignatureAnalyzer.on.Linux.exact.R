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



## Run SignatureAnalyzer extraction and attribution code.
## Note: This script is only tested on Linux
## and might not work well on Windows.
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){
    ## Run extraction and attribution
    ## SignatureAnalyzer needs to run 20 parallel runs,
    ## and pick the best run as the final result.
    out.dir <- paste0(topLevelFolder4Run,"/SignatureAnalyzer.results/",datasetName,"/seed.",seedInUse)
    if(file.exists(paste0(out.dir,"/best.run/sa.output.exp.csv"))) next

    message("\n\n########################################################\n\n")
    message(paste0("Begin running SignatureAnalyzer with maxK = 2",datasetName," using seed ",seedInUse,"...\n"))
    message("\n\n########################################################\n\n")

    SynSigRun:::SAMultiRunOneCatalog(
      num.runs = 20,
      signatureanalyzer.code.dir = paste0(usethis::proj_path(),"/data-raw/SignatureAnalyzer.052418"),
      input.catalog = paste0(topLevelFolder4Data,"/",datasetName,"/ground.truth.syn.catalog.csv"),
      out.dir = out.dir,
      maxK = 2,
      tol = 1e-7,
      test.only = FALSE,
      delete.tmp.files = TRUE,
      overwrite = FALSE,
      mc.cores = 20,
      verbose = FALSE,
      seed = seedInUse)

    SynSigRun:::CopyBestSignatureAnalyzerResult(
      sa.results.dir = paste0(datasetName,"/sp.sp/ExtrAttrExact/SignatureAnalyzer.results/seed.",seedInUse),
      verbose = TRUE,
      overwrite = FALSE)
  }
}
