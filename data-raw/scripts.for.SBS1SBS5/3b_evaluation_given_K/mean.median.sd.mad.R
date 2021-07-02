require(data.table)
require(gtools)


# Set working directory to the folder which contains top-level
# Summaries of results on SBS1-SBS5-correlated data sets with K specified as 2.
#
# PATH <- paste0("<path_to_results_on_SBS1-SBS5-correlated_datasets>/1a.Top_level_summary_for_K_unspecified")
# setwd(PATH)
oldWd <- getwd()
newWd <- "./1a.Top_level_summary_for_K_as_2/"
setwd(newWd)
load(paste0("FinalSummary.RDa"))


summaryList <- list()

for(measure in c("TPR","PPV","compositeMeasure")){

  dtCurrent <- FinalSummary$FinalExtr[[measure]]
  summaryList[[measure]] <- data.frame()

  for(tool in mixedsort(unique(dtCurrent$toolName))){
    nrows <- which(dtCurrent$toolName == tool)
    currValues <- dtCurrent$value[nrows]
    summaryList[[measure]] <- rbind(
      summaryList[[measure]],
      data.frame(
        "Name.of.Computational.Approach" = tool,
        mean = mean(currValues),
        median = stats::median(currValues),
        "standard.deviation" = stats::sd(currValues),
        "MAD" = stats::mad(currValues)
      )
    )
  }
}

## Total signatures extracted equals to
## TP + FP
{

  dtTP <- FinalSummary$FinalExtr$truePos
  dtFP <- FinalSummary$FinalExtr$falsePos
  summaryList[["totalSigs"]] <- data.frame()

  for(tool in mixedsort(unique(dtTP$toolName))){
    nrows <- which(dtTP$toolName == tool)
    currValues <- dtTP$value[nrows] + dtFP$value[nrows]
    summaryList[["totalSigs"]] <- rbind(
      summaryList[["totalSigs"]],
      data.frame(
        "Name.of.Computational.Approach" = tool,
        mean = mean(currValues),
        median = stats::median(currValues),
        "standard.deviation" = stats::sd(currValues),
        "MAD" = stats::mad(currValues)
      )
    )
  }
}



for(measure in c("cosSim","NumSigsSimilar")){
  summaryList[[measure]] <- list()
  for(gtSigName in c("SBS1","SBS5")){
    dtCurrent <- FinalSummary$FinalExtr[[measure]][[gtSigName]]
    summaryList[[measure]][[gtSigName]] <- data.frame()

    for(tool in mixedsort(unique(dtCurrent$toolName))){
      nrows <- which(dtCurrent$toolName == tool)
      currValues <- dtCurrent$value[nrows]
      summaryList[[measure]][[gtSigName]] <- rbind(
        summaryList[[measure]][[gtSigName]],
        data.frame(
          "Name.of.Computational.Approach" = tool,
          mean = mean(currValues),
          median = stats::median(currValues),
          "standard.deviation" = stats::sd(currValues),
          "MAD" = stats::mad(currValues)
        )
      )
    }
  }
}

dir.create("summary.tables", recursive = T)

for(measure in c("TPR","PPV","compositeMeasure")){
  data.table::fwrite(
    summaryList[[measure]],
    paste0("summary.tables/summary.of.",measure,".csv"))
}


for(gtSigName in c("SBS1","SBS5")){
  data.table::fwrite(
    summaryList$cosSim[[gtSigName]],
    paste0("summary.tables/summary.of.cosine.similarity.to.",gtSigName,".csv"))
}

for(gtSigName in c("SBS1","SBS5")){
  data.table::fwrite(
    summaryList$NumSigsSimilar[[gtSigName]],
    paste0("summary.tables/summary.of.number.extracted.sigs.resembling.",gtSigName,".csv"))
}

data.table::fwrite(
  summaryList$totalSigs,
  paste0("summary.tables/summary.of.total.number.of.sigs.extracted.csv"))

setwd(oldWd)
