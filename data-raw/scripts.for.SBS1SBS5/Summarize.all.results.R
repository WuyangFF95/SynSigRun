require(ICAMS)
require(SynSigEval)

## Specify default options
options(stringsAsFactors = F)

## Specify dataset names
slopes <- c(0.1,0.5,1,2,10)
Rsqs <- c(0.1,0.2,0.3,0.6)
datasetNames <- character(0)

for(slope in slopes)
  for(Rsq in Rsqs){
    datasetNames <- c(datasetNames,
                      paste0("S.",slope,".Rsq.",Rsq))
}

## Specify tool Names
# Tools which can do both extraction and attribution,
# excluding SP and SA
if(0){
  extrAttrToolNames <- c("hdp","sigfit.nmf","signeR",
                         "MutationalPatterns","SomaticSignatures")
  # Tools can do attribution,
  attrToolNames <- c("decompTumor2Sig","deconstructSigs","mSigAct",
                     "MutationalPatterns","mutSignatures",
                     "SignatureEstimation.QP","SignatureEstimation.SA",
                     "YAPSA")
} else{ ## Temporary workout
  extrAttrToolNames <- c("signeR","hdp","sigfit.nmf")
  # Tools can do attribution,
  attrToolNames <- c("decompTumor2Sig","deconstructSigs","mSigAct",
                     "MutationalPatterns","mutSignatures",
                     "SignatureEstimation.QP","SignatureEstimation.SA",
                     "YAPSA","sigfit.nmf")
}

## Specify seeds used in analysis.
## Specify 20 seeds used in software running
seedsInUse <- c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753)


## Part I: Run Summarize functions in SynSigEval
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    ## Summarize R-based Extraction and attribution tools.
    for(extrAttrToolName in extrAttrToolNames){
      SynSigEval::SummarizeSigOneExtrAttr96Subdir(
        run.dir = paste0(datasetName,"/sp.sp/ExtrAttr/",extrAttrToolName,
                            ".results/seed.",seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        overwrite = T)
    }
    ## Summarize R-based attribution-only tools.
    for(attrToolName in attrToolNames){
      SynSigEval::SummarizeSigOneAttr96Subdir(
        run.dir = paste0(datasetName,"/sp.sp/Attr/",attrToolName,
                            ".results/seed.",seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        overwrite = T)
    }
    ## Summarize sigproextractor
    SynSigEval::SummarizeSigOneSPSubdir(
      run.dir = paste0(datasetName,
                       "/sp.sp/ExtrAttr/sigproextractor.results/seed.",
                       seedInUse,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
    ## Summarize helmsman
    SynSigEval::SummarizeSigOnehelmsmanSubdir(
      run.dir = paste0(datasetName,
                       "/sp.sp/ExtrAttr/helmsman.results/seed.",
                       seedInUse,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
    ## Summarize SignatureAnalyzer
    SynSigEval:::SummarizeSigOneSASubdir(
      run.dir = paste0(datasetName,
                       "/sp.sp/ExtrAttr/SignatureAnalyzer.results/seed.",
                       seedInUse,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      which.run = "",overwrite = T)
  }
  ## Summarize maftools
  SynSigEval::SummarizeSigOneExtrAttr96Subdir(
    run.dir = paste0(datasetName,
                     "/sp.sp/ExtrAttr/maftools.results/seed.123456","/"),
    ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
    overwrite = T)
}

## Part II: Write summary table for 20 seeds for each tool with each dataset
for(datasetName in datasetNames){
  for(extrAttrToolName in c(extrAttrToolNames,"sigproextractor","helmsman","SignatureAnalyzer")){
    SynSigEval::SummarizeMultiRuns(
      tool.dir = paste0(datasetName,"/sp.sp/ExtrAttr/",extrAttrToolName,
                          ".results/"),
      run.names = paste0("seed.",seedsInUse))
  }
  SynSigEval::SummarizeMultiRuns(
    tool.dir = paste0(datasetName,"/sp.sp/ExtrAttr/","maftools",
                      ".results/"),
    run.names = paste0("seed.","123456"))
  for(attrToolName in attrToolNames){
    SynSigEval::SummarizeMultiRuns(
      tool.dir = paste0(datasetName,"/sp.sp/Attr/",attrToolName,
                        ".results/"),
      run.names = paste0("seed.",seedsInUse))
  }
}

## Part III: Write Summary table of multiple tools for each dataset
for(datasetName in datasetNames){
    SynSigEval::SummarizeMultiToolsOneDataset(
      third.level.dir = paste0(datasetName,"/sp.sp/ExtrAttr/"),
      tool.dirnames = paste0(c("sigproextractor","SignatureAnalyzer","helmsman","maftools",
                               extrAttrToolNames),".results"))
    SynSigEval::SummarizeMultiToolsOneDataset(
      third.level.dir = paste0(datasetName,"/sp.sp/Attr/"),
      tool.dirnames = paste0(attrToolNames,".results"))
}

## Part IV: Generate a combined summary table for results
## of multiple datasets, from multiple tools
SummarizeMultiToolsMultiDatasets(dataset.dirs = datasetNames,
                       second.third.level.dirname = "sp.sp/ExtrAttr",
                       out.dir = "./FinalExtrAttrSummary", overwrite = T)
SummarizeMultiToolsMultiDatasets(dataset.dirs = datasetNames,
                       second.third.level.dirname = "sp.sp/Attr",
                       out.dir = "./FinalAttrSummary", overwrite = T)


