
## Load required packages
library(ICAMS)
library(SynSigRun)
library(mutSignatures)



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
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){
    RunmutSignatures(
      input.catalog = paste0(datasetName, "/sp.sp/ground.truth.syn.catalog.csv"),
      out.dir = paste0(datasetName, "/sp.sp/Attr/mutSignatures.results/seed.", seedInUse),
      CPU.cores = 30,
      iterations = 1000,
      seedNumber = seedInUse,
      K.range = c(1,10),
      overwrite = T)


    RunmutSignaturesAttributeOnly(
      input.catalog = paste0(datasetName, "/sp.sp/ground.truth.syn.catalog.csv"),
      gt.sigs.file = paste0(datasetName, "/sp.sp/ground.truth.syn.sigs.csv"),
      out.dir = paste0(datasetName, "/sp.sp/Attr/mutSignatures.results/seed.", seedInUse),
      seedNumber = seedInUse,
      overwrite = T
    )
  }
}

