require(SynSig)


## Specify dataset names
slopes <- c(0.1,0.5,1)
Rsqs <- c(0.1,0.2,0.3,0.6)
for(slope in slopes)
  for(Rsq in Rsqs){

    datasetName <- paste0("S.",slope,".Rsq.",Rsq)

    gt.syn.catalog <- ICAMS:::ReadCatSNS96(path = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
                                           strict = FALSE)

    CreateMultiModalMuSigOutput(catalog = gt.syn.catalog,
                                out.dir = paste0(datasetName,"/sp.sp/MultiModalMuSig.results"),
                                overwrite = T)
  }

##

