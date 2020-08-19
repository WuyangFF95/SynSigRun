
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
## Specify tool Names
# Here, toolNames refer to runs on different eta 
# for MultiModalMuSig.LDA, which are "MMLDA.alpha.0.1.eta.<eta>".
etas <- c("0.1", "0.2", "0.3", "0.4",
          "0.5", "0.6", "0.75", "1.0",
          "2.5", "5.0", "7.5", "10.0",
          "15.0", "20.0", "30.0", "50.0", "75.0")

otherExtrAttrToolNames <- paste0("MMLDA.alpha.0.1.eta.",etas)


## Specify seeds used in analysis.
## Specify 20 seeds used in software running
seedsInUse <- c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753)



## Part I: Run Summarize functions in SynSigEval
for(datasetName in datasetNames){
  for(seedInUse in seedsInUse){
    ## Summarize MultiModalMuSig.LDA results
    for(extrAttrToolName in otherExtrAttrToolNames){
      SynSigEval::SummarizeSigOneExtrAttrSubdir(
        run.dir = paste0(datasetName,"/sp.sp/ExtrAttrExactMMLDATest/",extrAttrToolName,
                         ".results/seed.",seedInUse,"/"),
        ground.truth.exposure.dir = paste0(datasetName,"/sp.sp/"),
        overwrite = T)
    }
  }
}

## Part II: Write summary table for 20 seeds for each tool with each dataset
for(datasetName in datasetNames){
  ## For each dataset, summarize 20 runs
  ## using different seeds by EMu
  for(extrAttrToolName in otherExtrAttrToolNames){
    SynSigEval::SummarizeMultiRuns(
      datasetName = datasetName,
      toolName = extrAttrToolName,
      tool.dir = paste0(datasetName,"/sp.sp/ExtrAttrExactMMLDATest/",extrAttrToolName,
                          ".results/"),
      run.names = paste0("seed.",seedsInUse))
  }
}

## Part III: Write Summary table of multiple tools for each dataset
for(slope in slopes){
  for(Rsq in Rsqs){
    datasetName <- paste0("S.",slope,".Rsq.",Rsq)
    SynSigEval::SummarizeMultiToolsOneDataset(
      third.level.dir = paste0(datasetName,"/sp.sp/ExtrAttrExactMMLDATest/"),
      toolName = otherExtrAttrToolNames,
      tool.dirnames = paste0(c(otherExtrAttrToolNames),".results"),
      datasetGroup = Rsq,
      datasetGroupName = "Pearson's R^2",
      datasetSubGroup = slope,
      datasetSubGroupName = "SBS1:SBS5 mutation count ratio")
  }
}
## Part IV: Generate a summary table and boxplot for results
## of multiple datasets, from each separate tool.
datasetGroup <- rep(c("0.1","0.2","0.3","0.6"),5)
names(datasetGroup) <- datasetNames
datasetSubGroup <- rep(c("0.1","0.5","1","2","10"),each = 4)
names(datasetSubGroup) <- datasetNames


for(toolName in otherExtrAttrToolNames){
  SummarizeOneToolMultiDatasets(
    dataset.dirs = datasetNames,
    datasetGroup = datasetGroup,
    datasetGroupName = "Pearson's R^2",
    datasetSubGroup = datasetSubGroup,
    datasetSubGroupName = "SBS1:SBS5 mutation count ratio",
    toolName = toolName,
    tool.dirname = paste0("sp.sp/ExtrAttrExactMMLDATest/",toolName,".results/"),
    out.dir = paste0("FinalToolWiseSummary/ExtrAttrExactMMLDATest/",toolName,"/"),
    overwrite = T)
}



## Part V: Generate a combined summary table for results
## of multiple datasets, from multiple tools
FinalExtrAttr <- SummarizeMultiToolsMultiDatasets(
  dataset.dirs = datasetNames,
  second.third.level.dirname = "sp.sp/ExtrAttrExactMMLDATest/",
  out.dir = "./FinalExtrAttrExactMMLDATestSummary", overwrite = T)


