
# Run Convert.EMu.Results.R and
# Convert.MultiModalMuSig.Results.R
# before running this script.


# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
# SynSigRun home can be retrieved by usethis::proj_path
#
# PATH <- paste0(usethis::proj_path,"/data-raw/scripts.for.SBS1SBS5")
# setwd(PATH)


require(ICAMS)
require(SynSigEval)

## Specify default options
options(stringsAsFactors = F)

## Specify dataset names
slopes <- c("0.1","0.5","1","2","10")
Rsqs <- c("0.1","0.2","0.3","0.6")
datasetNames <- character(0)

for(slope in slopes){
  for(Rsq in Rsqs){
    datasetNames <- c(datasetNames,
                      paste0("S.",slope,".Rsq.",Rsq))
  }
}
## Specify names of computational approaches
# R-based tools which can do both extraction and attribution,
# excluding SignatureAnalyzer (due to special folder structure)
# and maftools (its seed is hard-coded)
RBasedExtrAttrToolNames <- c("hdp",
                             "mutSpec.NMF",
                             "sigfit.EMu","sigfit.NMF",
                             "signature.tools.lib","signeR",
                             "SomaticSignatures.NMF",
                             "TCSM")
# Python or other language based tools.
# excluding maftools and MutationalPatterns (seed is hard-coded) and EMu (cannot designate seed)
otherExtrAttrToolNames <- c("MultiModalMuSig.CTM","MultiModalMuSig.LDA")

# List computational approaches with hard-coded seed or do not accept seeds.
toolNameWOSeed <- "EMu"
toolNamesWFixedSeed <- c("maftools","MutationalPatterns")


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
    for(extrAttrToolName in RBasedExtrAttrToolNames){
      SynSigEval::SummarizeSigOneExtrAttrSubdir(
        run.dir = paste0(datasetName,"/sp.sp/ExtrAttrExact/",extrAttrToolName,
                            ".results/seed.",seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        overwrite = T)
    }
    ## Summarize non-R Extraction and attribution tools.
    for(extrAttrToolName in otherExtrAttrToolNames){
      SynSigEval::SummarizeSigOneExtrAttrSubdir(
        run.dir = paste0(datasetName,"/sp.sp/ExtrAttrExact/",extrAttrToolName,
                         ".results/seed.",seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        overwrite = T)
    }
    ## Summarize SigProExtractor
    SynSigEval::SummarizeSigOneSigProExtractorSubdir(
      run.dir = paste0(datasetName,
                       "/sp.sp/ExtrAttrExact/SigProExtractor.results/seed.",
                       seedInUse,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
    ## Summarize helmsman.NMF
    SynSigEval::SummarizeSigOnehelmsmanSubdir(
      run.dir = paste0(datasetName,
                       "/sp.sp/ExtrAttrExact/helmsman.NMF.results/seed.",
                       seedInUse,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
    ## Summarize SignatureAnalyzer
    {
      SynSigEval::CopyBestSignatureAnalyzerResult(
        sa.results.dir = paste0(datasetName,
                                "/sp.sp/ExtrAttrExact/SignatureAnalyzer.results/seed.",
                                seedInUse,"/"),
        overwrite = T)

      SynSigEval:::SummarizeSigOneSASubdir(
        run.dir = paste0(datasetName,
                         "/sp.sp/ExtrAttrExact/SignatureAnalyzer.results/seed.",
                         seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        which.run = "/best.run/",
        overwrite = T)
    }
  }
  ## Summarize maftools and MutationalPatterns which have a hard-coded seed.
  for(toolName in toolNamesWFixedSeed) {
    SynSigEval::SummarizeSigOneExtrAttrSubdir(
      run.dir = paste0(datasetName,
                       "/sp.sp/ExtrAttrExact/",toolName,".results/seed.123456","/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
  }
  ## Summarize EMu
  for(nrun in 1:20){
    SynSigEval::SummarizeSigOneExtrAttrSubdir(
      run.dir = paste0(datasetName,"/sp.sp/ExtrAttrExact/","EMu",
                       ".results/run.",nrun,"/"),
      ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
      overwrite = T)
  }
}

## Part II: Write summary table for 20 seeds for each tool with each dataset
otherExtrAttrToolNames <- c("helmsman.NMF","MultiModalMuSig.CTM","MultiModalMuSig.LDA","SigProExtractor","SignatureAnalyzer")
for(datasetName in datasetNames){
  ## For each dataset, summarize 20 runs
  ## using different seeds by EMu
  for(extrAttrToolName in c(
    RBasedExtrAttrToolNames,otherExtrAttrToolNames)){
    SynSigEval::SummarizeMultiRuns(
      datasetName = datasetName,
      toolName = extrAttrToolName,
      tool.dir = paste0(datasetName,"/sp.sp/ExtrAttrExact/",extrAttrToolName,
                          ".results/"),
      run.names = paste0("seed.",seedsInUse))
  }
  ## For each dataset, summarize 20 runs
  ## (without seeds) by EMu
  SynSigEval::SummarizeMultiRuns(
    datasetName = datasetName,
    toolName = "EMu",
    tool.dir = paste0(datasetName,"/sp.sp/ExtrAttrExact/","EMu",
                      ".results/"),
    run.names = paste0("run.",1:20))
  ## For each dataset, summarize 1 run by maftools and MutationalPatterns
  for(toolName in toolNamesWFixedSeed){
    SynSigEval::SummarizeMultiRuns(
      datasetName = datasetName,
      toolName = toolName,
      tool.dir = paste0(datasetName,"/sp.sp/ExtrAttrExact/",toolName,
                        ".results/"),
      run.names = paste0("seed.","123456"))
  }
}

## Part III: Write Summary table of multiple tools for each dataset
for(slope in slopes){
  for(Rsq in Rsqs){
    datasetName <- paste0("S.",slope,".Rsq.",Rsq)
    SynSigEval::SummarizeMultiToolsOneDataset(
      third.level.dir = paste0(datasetName,"/sp.sp/ExtrAttrExact/"),
      toolName = c(RBasedExtrAttrToolNames,otherExtrAttrToolNames,
                   "EMu",toolNamesWFixedSeed),
      tool.dirnames = paste0(c(RBasedExtrAttrToolNames,otherExtrAttrToolNames,
                               "EMu",toolNamesWFixedSeed),".results"),
      datasetGroup = Rsq,
      datasetGroupName = "Pearson's R^2",
      datasetSubGroup = slope,
      datasetSubGroupName = "SBS1:SBS5 mutation count ratio")
  }
}
## Part IV: Generate a summary table and boxplot for results
## of multiple datasets, from each separate tool.
datasetGroup <- rep(c(0.1,0.2,0.3,0.6),5)
names(datasetGroup) <- datasetNames
datasetSubGroup <- rep(c(0.1,0.5,1,2,10),each = 4)
names(datasetSubGroup) <- datasetNames


for(toolName in c(
  RBasedExtrAttrToolNames,otherExtrAttrToolNames,
  toolNamesWFixedSeed,"EMu")){
  SummarizeOneToolMultiDatasets(
    dataset.dirs = datasetNames,
    datasetGroup = datasetGroup,
    datasetGroupName = "Pearson's R^2",
    datasetSubGroup = datasetSubGroup,
    datasetSubGroupName = "SBS1:SBS5 mutation count ratio",
    toolName = toolName,
    tool.dirname = paste0("sp.sp/ExtrAttrExact/",toolName,".results/"),
    out.dir = paste0("FinalToolWiseSummary/ExtrAttrExact/",toolName,"/"),
    overwrite = T)
}



## Part V: Generate a combined summary table for results
## of multiple datasets, from multiple tools
FinalExtrAttrExact <- SummarizeMultiToolsMultiDatasets(
  dataset.dirs = datasetNames,
  second.third.level.dirname = "sp.sp/ExtrAttrExact",
  out.dir = "./FinalExtrAttrExactSummary", overwrite = T)


