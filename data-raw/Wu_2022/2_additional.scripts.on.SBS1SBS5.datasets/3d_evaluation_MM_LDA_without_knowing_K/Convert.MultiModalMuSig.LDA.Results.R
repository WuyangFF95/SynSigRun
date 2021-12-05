
## Run 2a_running_approaches_without_knowing_K/Run.MultiModalMuSig.jl
## and 2b_running_approaches_given_K/Run.MultiModalMuSig.jl
## Before running this script.
##
## Run this script before running Summarize.all.results.R



# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)



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
                
## Specify etas
etas <- c("0.1", "0.2", "0.3", "0.4",
          "0.5", "0.6", "0.75", "1.0",
          "2.5", "5.0", "7.5", "10.0",
          "15.0", "20.0", "30.0", "50.0", "75.0")


## After running, convert MultiModalMuSig-formatted tsv files
## to SynSigEval/ICAMS-formatted csv files.
for(datasetName in datasetNames){
  for(eta in etas){
    for(seedInUse in seedsInUse){
        ## Convert signatures
        signatures <- SynSigEval::MMCatalog2ICAMS(
          paste0(datasetName,
                 "/sp.sp/ExtrAttrMMLDATest/MMLDA.alpha.0.1.eta.",eta,".results/",
                 "seed.",seedInUse,
                 "/extracted.signatures.tsv"),
          region = "unknown",
          catalog.type = "counts.signature")
        ICAMS::WriteCatalog(
          signatures,
          paste0(datasetName,
                 "/sp.sp/ExtrAttrMMLDATest/MMLDA.alpha.0.1.eta.",eta,".results/",
                 "seed.",seedInUse,
                 "/extracted.signatures.csv"))

        ## Convert exposures
        exposure <- SynSigEval::ReadExposureMM(
          paste0(datasetName,
                 "/sp.sp/ExtrAttrMMLDATest/MMLDA.alpha.0.1.eta.",eta,".results/",
                 "seed.",seedInUse,
                 "/inferred.exposures.tsv"))
        SynSigGen::WriteExposure(
          exposure,
          paste0(datasetName,
                 "/sp.sp/ExtrAttrMMLDATest/MMLDA.alpha.0.1.eta.",eta,".results/",
                 "seed.",seedInUse,
                 "/inferred.exposures.csv"))
    }
  }
}

