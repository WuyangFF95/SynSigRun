
# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)

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
                
## List all eta values we want to test for MM.LDA
alpha <- 0.1
etas <- c(2.5,"5.0",7.5,"10.0","15.0","20.0","30.0","50.0","75.0")

## List all alpha values we want to test for MM.CTM
alphas <- c(2.5,"5.0",7.5,"10.0","15.0","20.0","30.0","50.0","75.0")

## Create MultimodalMuSig-formatted input-catalog under
## <dataset.top.level.dir>/sp.sp/ExtrAttrMMLDATest/
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    for(eta in etas) {
      CreateMultiModalMuSigOutput(
        catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
        out.dir = paste0(datasetName,"/sp.sp/ExtrAttrMMLDATest/MMLDA.alpha.",
                         alpha,".eta.",eta,".results",
                         "/seed.",seedInUse),
        overwrite = T)
    
      CreateMultiModalMuSigOutput(
        catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
        out.dir = paste0(datasetName,"/sp.sp/ExtrAttrExactMMLDATest/MMLDA.alpha.",
                         alpha,".eta.",eta,".results",
                         "/seed.",seedInUse),
        overwrite = T)
    }

    for(alpha in alphas) {
      CreateMultiModalMuSigOutput(
        catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
        out.dir = paste0(datasetName,"/sp.sp/ExtrAttrMMCTMTest/MMCTM.alpha.",
                         alpha,".results",
                         "/seed.",seedInUse),
        overwrite = T)
      CreateMultiModalMuSigOutput(
        catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
        out.dir = paste0(datasetName,"/sp.sp/ExtrAttrExactMMCTMTest/MMCTM.alpha.",
                                 alpha,".results",
                         "/seed.",seedInUse),
        overwrite = T)
    }
  }
}