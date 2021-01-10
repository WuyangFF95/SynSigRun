require(data.table)

oldWd <- getwd()
setwd("../../practice/3_Signature_Challenge/3.2_SBS1-SBS5-Correlation_Project/new_results/FinalExtrAttrSummary")
load("FinalSummary.RDa")

dtList <- list()
summaryList <- list()

for(measure in c("TPR","PPV","compositeMeasure")){

  dtList[[measure]] <- FinalSummary$FinalExtr[[measure]]
  summaryList[[measure]] <- data.frame()

  for(tool in unique(dtList[[measure]]$toolName)){
    nrows <- which(dtList[[measure]]$toolName == tool)
    currValues <- dtList[[measure]]$value[nrows]
    summaryList[[measure]] <- rbind(
      summaryList[[measure]],
      data.frame(
        "Name.of.Computational.Approach" = tool,
        mean = mean(currValues),
        median = stats::median(currValues),
        "standard.deviation" = stats::sd(currValues),
        "median.absolute.deviation" = stats::mad(currValues, constant = 1),
        "scaled.MAD" = stats::mad(currValues)
      )
    )
  }
}



for(gtSigName in c("SBS1","SBS5")){

  dtList[[gtSigName]] <- FinalSummary$FinalExtr$cosSim[[gtSigName]]
  summaryList[[gtSigName]] <- data.frame()

  for(tool in unique(dtList[[gtSigName]]$toolName)){
    nrows <- which(dtList[[gtSigName]]$toolName == tool)
    currValues <- dtList[[gtSigName]]$value[nrows]
    summaryList[[gtSigName]] <- rbind(
      summaryList[[gtSigName]],
      data.frame(
        "Name.of.Computational.Approach" = tool,
        mean = mean(currValues),
        median = stats::median(currValues),
        "standard.deviation" = stats::sd(currValues),
        "median.absolute.deviation" = stats::mad(currValues, constant = 1),
        "scaled.MAD" = stats::mad(currValues)
      )
    )
  }
}


for(measure in c("TPR","PPV","compositeMeasure")){
  data.table::fwrite(summaryList[[measure]],
                     paste0("summary.of.",measure,".csv") )
}

for(gtSigName in c("SBS1","SBS5")){
  data.table::fwrite(summaryList[[measure]],
                     paste0("summary.of.cosine.similarity.to.",gtSigName,".csv") )
}

setwd(oldWd)
