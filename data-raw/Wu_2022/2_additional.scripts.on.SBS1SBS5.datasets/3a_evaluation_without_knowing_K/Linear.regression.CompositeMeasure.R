
## Require tidyverse for using "%>%"
require(tidyverse)
require(lattice)

setwd("../../practice/3_Signature_Challenge/3.2_SBS1-SBS5-Correlation_Project/new_results/FinalToolWiseSummary/ExtrAttr/")


## Specify slopes as datasetSubGroup values
## Specify Rsqs as datasetGroup values
## In order to regress correctly, slopes and Rsqs
## must be set as numeric
slopes <- c(0.1,0.5,1,2,10)
Rsqs <- c(0.1,0.2,0.3,0.6)

## Specify seeds used in analysis.
## Specify 20 seeds used in software running
## To correctly regress, must be set as character.
seedsInUse <- as.character(c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476, 20897,
                27847, 34637, 49081, 75679, 103333,
                145879, 200437, 310111, 528401, 1076753))



extrAttrToolNames <-
  c("hdp","MutationalPatterns","sigfit.EMu",
    "sigfit.NMF","signeR","TCSM",
    "helmsman.NMF","MultiModalMuSig.CTM",
    "MultiModalMuSig.LDA","SigProExtractor","SignatureAnalyzer")
toolNameWOSeed <- "EMu"
toolNameWFixedSeed <- "maftools"
composites4Regression <- list()

## Linear regression for the MEDIAN performance
## of each computational approach.
##
## Extraction performance: COMPOSITE measure
##
## Attribution performance: Manhattan distance
##
for(toolName in extrAttrToolNames){

  ## Load one tool summary
  load(paste0(toolName,"/OneToolSummary.RDa"))
  sigNames <- setdiff(names(OneToolSummary$cosSim),"combined")

  ## Calculate COMPOSITE measure for each Rsq and slope and seed
  composite4Regression <- data.frame()

  for(slope in slopes){
    for(Rsq in Rsqs){

      for (seedInUse in seedsInUse){

        values <- numeric(length = 4)
        names(values) <- c("SBS1","SBS5","TPR","1-FDR")


        TPR <-  OneToolSummary$TPR %>%
          dplyr::filter(seed == paste0("seed.",seedInUse),
                        datasetGroup == Rsq,
                        datasetSubGroup == slope)
        values["TPR"] <- TPR$value


        FDR <- OneToolSummary$FDR %>%
          dplyr::filter(seed == paste0("seed.",seedInUse),
                        datasetGroup == Rsq,
                        datasetSubGroup == slope)
        values["1-FDR"] <- 1 - FDR$value

        cosSimSBS1 <-  OneToolSummary$cosSim$SBS1 %>%
          dplyr::filter(seed == paste0("seed.",seedInUse),
                        datasetGroup == Rsq,
                        datasetSubGroup == slope)
        values["SBS1"] <- cosSimSBS1$value


        cosSimSBS5 <- OneToolSummary$cosSim$SBS5 %>%
          dplyr::filter(seed == paste0("seed.",seedInUse),
                        datasetGroup == Rsq,
                        datasetSubGroup == slope)
        values["SBS5"] <- cosSimSBS5$value

        current <- data.frame(
          seed = seedInUse,
          datasetGroup = Rsq,
          datasetSubGroup = slope,
          value = sum(values)
        )

        composite4Regression <-
          rbind(composite4Regression, current)
      }
    }
  }

  utils::write.csv(
    composite4Regression,
    file = paste0(toolName,"/COMPOSITE4Regression.csv"),
    row.names = F)

  composites4Regression[[toolName]] <- composite4Regression
}


## Do linear regression and plotting for Composite measure against R^2.
models2Rsq <- list()

for(toolName in extrAttrToolNames){

  models2Rsq[[toolName]] <- lm(value ~ datasetGroup, data = composites4Regression[[toolName]])
  summary(models2Rsq[[toolName]])

  current <- composites4Regression[[toolName]]

  pdf(paste0(toolName,"/CompositeMeasure2Rsq.pdf"))

  #### Page 1: Plot general trend diagnostics
  ## Plot data points and fitted line
  plot(y = current$value, x  = current$datasetGroup, pch = 20, col = "blue",
       main = paste0(toolName,": composite measure AGAINST R^2"),
       xlab = "Pearson's R^2", ylab = "Composite measure for extraction")
  abline(models2Rsq[[toolName]], lty = "dotted")

  ## Plot line chart for median of each Pearson's R^2.
  medians <- numeric(0)
  means <- numeric(0)

  for(Rsq in Rsqs){
      currValues <- current %>%
        dplyr::filter(datasetGroup == Rsq) %>%
        dplyr::select(value)

      currMedian <- median(currValues[,1])
      currMean <- mean(currValues[,1])
      names(currMedian) <- Rsq
      names(currMean) <- Rsq
      medians <- c(medians, currMedian)
      means <- c(means, currMean)
  }
  lines(x = names(medians),y = medians,
        type = "o",pch = 4,col = "red")
  lines(x = names(means),y = means,
        type = "o",pch = 4,col = "blue")

  #### Page 2: Plot trend diagnostics of different Rsq in lattice plot
  ## Plot data points and fitted line
  latticeObj <- lattice::xyplot(
    value ~ datasetGroup | as.factor(datasetSubGroup),
    data = current,
    type=c("p", "r"),
    pch = 20, col = "blue",
    lty = "dotted",
    main = paste0("composite measure AGAINST R^2 given slope"),
    xlab = "Pearson's R^2", ylab = "Composite measure for extraction")

  plot(latticeObj)
  grDevices::dev.off()

}


## Do linear regression and plotting for Composite measure
## against SBS1:SBS5 count ratio.
models2slope <- list()

for(toolName in extrAttrToolNames){

  models2slope[[toolName]] <- lm(value ~ datasetSubGroup, data = composites4Regression[[toolName]])
  summary(models2slope[[toolName]])

  current <- composites4Regression[[toolName]]

  pdf(paste0(toolName,"/CompositeMeasure2SBS1SBS5ratio.pdf"))

  #### Page 1: plot overall trend diagnostics
  ## Plot data points and fitted line
  plot(y = current$value, x  = current$datasetSubGroup, pch = 20, col = "blue",
       main = paste0(toolName,": composite measure AGAINST SBS1:SBS5 ratio"),
       xlab = "SBS1:SBS5 mutation count ratio", ylab = "Composite measure for extraction")
  abline(models2slope[[toolName]], lty = "dotted")
  ## Plot line chart for median of each Pearson's R^2.
  medians <- numeric(0)
  means <- numeric(0)
  for(slope in slopes){
    currValues <- current %>%
      dplyr::filter(datasetSubGroup == slope) %>%
      dplyr::select(value)
    currMedian <- median(currValues[,1])
    currMean <- mean(currValues[,1])
    names(currMedian) <- slope
    names(currMean) <- slope
    medians <- c(medians, currMedian)
    means <- c(means, currMean)
  }
  lines(x = names(medians),y = medians,
        type = "o",pch = 4,col = "red")
  lines(x = names(means),y = means,
        type = "o",pch = 4,col = "blue")

  #### Page 2: Plot trend diagnostics of different slopes
  ## in lattice bivariate xyplot with multi-facets
  ## Plot data points and fitted line
  latticeObj <- lattice::xyplot(
    value ~ datasetSubGroup | as.factor(datasetGroup),
    data = current,
    type=c("p", "r"), ## Plot points and regression line.
    pch = 20, col = "blue",
    lty = "dotted",
    main = paste0("composite measure AGAINST SBS1:SBS5 ratio given R^2"),
    xlab = "SBS1:SBS5 count ratio", ylab = "Composite measure for extraction")

  plot(latticeObj)
  grDevices::dev.off()

}

save.image("Linear.regression.Composite.Rdata")

