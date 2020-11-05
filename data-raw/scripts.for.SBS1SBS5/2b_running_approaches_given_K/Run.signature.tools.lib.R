# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)


## Load required packages
library(ICAMS)
library(signature.tools.lib)
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



## Run Extraction and attribution packages
## sigproextractor (Python package) and MultiModalMuSig (Julia package)
## needs to be run with external script.

## NMF will have an error if K.range includes 1!
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){

    cat("\n===========================================\n")
    cat(paste0("Running signature.tools.lib\ on data set ",datasetName," using seed ",seedInUse,"...\n"))
    cat("\n===========================================\n")

    input.catalog = paste0(datasetName, "/sp.sp/ground.truth.syn.catalog.csv")
    out.dir = paste0(datasetName, "/sp.sp/ExtrAttrExact/signature.tools.lib.results/seed.", seedInUse,"/")
    CPU.cores = 10
    K.exact = 2

    spectra <- ICAMS::ReadCatalog(
	  input.catalog,
	  strict = FALSE)
    ## convSpectra: convert the ICAMS-formatted spectra catalog
    ## into a matrix which signature.tools.lib accepts:
    ## 1. Remove the catalog related attributes in convSpectra
    ## 2. Transpose the catalog
    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL


    retval <- signature.tools.lib::SignatureExtraction(
      cat = convSpectra,
      outFilePath = out.dir,
      nrepeats = 200,
      nboots = 20,
      nparallel = CPU.cores,
      nsig = K.exact,
      mut_thr = 0,
      type_of_extraction = "generic",
      project = "test",
      parallel = TRUE,
      nmfmethod = "brunet")

  save(retval,paste0(out.dir,"/retval.Rdata"))
  }
}
