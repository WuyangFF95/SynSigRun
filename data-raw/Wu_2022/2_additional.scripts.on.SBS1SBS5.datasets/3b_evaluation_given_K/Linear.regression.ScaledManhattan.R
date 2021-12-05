
## Require tidyverse for using "%>%"
require(tidyverse)
require(lattice)

setwd("../../practice/3_Signature_Challenge/3.2_SBS1-SBS5-Correlation_Project/new_results/FinalToolWiseSummary/ExtrAttrExact/")


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
  c("hdp","sigfit.EMu",
    "sigfit.NMF","signeR","TCSM",
    "helmsman.NMF","MultiModalMuSig.CTM",
    "MultiModalMuSig.LDA","SigProExtractor","SignatureAnalyzer")
toolNameWOSeed <- "EMu"
toolNameWFixedSeed <- c("maftools","MutationalPatterns")

for(toolName in c(extrAttrToolNames,toolNameWFixedSeed,toolNameWOSeed)){
  dirToBeCreated <- paste0("../../TrendDiag/ExtrAttrExact/",toolName,"/")
  if(!dir.exists(dirToBeCreated))
    dir.create(dirToBeCreated,recursive = T)
}


## Linear regression for the MEDIAN performance
## of each computational approach.
##
## Extraction performance: COMPOSITE measure
##
## Attribution performance: Manhattan distance
##
ForRegression <- list()
for(toolName in c(extrAttrToolNames,toolNameWFixedSeed,toolNameWOSeed)){

  ## Load one tool summary
  load(paste0(toolName,"/OneToolSummary.RDa"))
  sigNames <- setdiff(names(OneToolSummary$cosSim),"combined")


  ForRegression[[toolName]] <- list()

  ForRegression[[toolName]]$ManhattSBS1 <-  OneToolSummary$ManhattanDist$SBS1 %>%
    dplyr::select(seed, value, datasetGroup, datasetSubGroup)

  ForRegression[[toolName]]$ManhattSBS5 <- OneToolSummary$ManhattanDist$SBS5 %>%
    dplyr::select(seed, value, datasetGroup, datasetSubGroup)

  indexes <- c("ManhattSBS1","ManhattSBS5")

  for(index in indexes){

    ForRegression[[toolName]][[index]]$seed <-
      gsub("seed\\.","",ForRegression[[toolName]][[index]]$seed)

    ForRegression[[toolName]][[index]]$datasetGroup <-
      as.numeric(paste(ForRegression[[toolName]][[index]]$datasetGroup))

    ForRegression[[toolName]][[index]]$datasetSubGroup <-
      as.numeric(paste(ForRegression[[toolName]][[index]]$datasetSubGroup))

    utils::write.csv(
      ForRegression[[toolName]][[index]],
      file = paste0(toolName,"/",index,"ForRegression.csv"),
      quote = F,
      row.names = F)
  }

}


## Do linear regression and plotting for Composite measure against R^2.
models2Rsq <- list()

for(toolName in c(extrAttrToolNames,toolNameWFixedSeed,toolNameWOSeed)){

  for(index in indexes)  {

    models2Rsq[[toolName]] <-
      lm(value ~ datasetGroup,
         data = ForRegression[[toolName]][[index]])
    summary(models2Rsq[[toolName]])

    current <- ForRegression[[toolName]][[index]]

    pdf(paste0("../../TrendDiag/ExtrAttrExact/",toolName,"/",index,"AgainstRsq.pdf"))

    #### Page 1: Plot general trend diagnostics
    ## Plot data points and fitted line
    plot(y = current$value, x  = current$datasetGroup, pch = 20, col = "blue",
         main = paste0(toolName,":",index," AGAINST R^2"),
         xlab = "Pearson's R^2", ylab = index)
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
      panel = function(x, y, ...){
        panel.xyplot(x, y,
                     type=c("p", "r"), ## Plot points and regression line.
                     pch = 20, col = "black",
                     lty = "dotted", ...)
        panel.linejoin(x, y,
                       fun = mean,
                       horizontal = FALSE,
                       pch = 4, col = "blue",
                       lty = "dashed", ...)
        panel.linejoin(x, y,
                       fun = median,
                       horizontal = FALSE,
                       pch = 4, col = "red",
                       lty = "dashed", ...)

      },
      main = paste0(index," AGAINST R^2 given slope"),
      xlab = "Pearson's R^2", ylab = index)

    plot(latticeObj)
    grDevices::dev.off()

  }

}


## Do linear regression and plotting for Composite measure
## against SBS1:SBS5 count ratio.
models2slope <- list()

for(toolName in c(extrAttrToolNames,toolNameWFixedSeed,toolNameWOSeed)){
  for(index in indexes){

    models2slope[[toolName]] <-
      lm(value ~ datasetSubGroup,
         data = ForRegression[[toolName]][[index]])
    summary(models2slope[[toolName]])

    current <- ForRegression[[toolName]][[index]]

    pdf(paste0("../../TrendDiag/ExtrAttrExact/",toolName,"/",index,"AgainstSBS1SBS5ratio.pdf"))

    #### Page 1: plot overall trend diagnostics
    ## Plot data points and fitted line
    plot(y = current$value, x  = current$datasetSubGroup, pch = 20, col = "blue",
         main = paste0(toolName,":",index," AGAINST SBS1:SBS5 ratio"),
         xlab = "SBS1:SBS5 mutation count ratio", ylab = index)
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
      panel = function(x, y, ...){
        panel.xyplot(x, y,
                     type=c("p", "r"), ## Plot points and regression line.
                     pch = 20, col = "black",
                     lty = "dotted", ...)
        panel.linejoin(x, y,
                       fun = mean,
                       horizontal = FALSE,
                       pch = 4, col = "blue",
                       lty = "dashed", ...)
        panel.linejoin(x, y,
                       fun = median,
                       horizontal = FALSE,
                       pch = 4, col = "red",
                       lty = "dashed", ...)

      },
      main = paste0(index," AGAINST SBS1:SBS5 ratio given R^2"),
      xlab = "SBS1:SBS5 count ratio", ylab = index)

    plot(latticeObj)
    grDevices::dev.off()

  }
}
save.image("../../TrendDiag/ExtrAttrExact/Linear.regression.Manhattan.Distances.Rdata")

