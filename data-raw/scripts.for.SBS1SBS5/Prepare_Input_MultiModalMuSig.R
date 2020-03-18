require(SynSigRun)


## Specify dataset names
slopes <- c(0.1,0.5,1,2,10)
Rsqs <- c(0.1,0.2,0.3,0.6)

## Specify 20 seeds used in software running
seedsInUse <- c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753)


## Create an input file for each dataset in each folder
## containing runs of different seeds.
for(slope in slopes)
  for(Rsq in Rsqs){

    datasetName <- paste0("S.",slope,".Rsq.",Rsq)

    gt.syn.catalog <- ICAMS:::ReadCatalog(
      paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
      region = "genome",
      catalog.type = "counts",
      strict = FALSE)

    for(seedInUse in seedsInUse){
    CreateMultiModalMuSigOutput(catalog = gt.syn.catalog,
                                out.dir = paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.LDA.results/seed.",seedInUse),
                                overwrite = T)
    CreateMultiModalMuSigOutput(catalog = gt.syn.catalog,
                                out.dir = paste0(datasetName,"/sp.sp/ExtrAttr/MultiModalMuSig.CTM.results/seed.",seedInUse),
                                overwrite = T)
    }
}

