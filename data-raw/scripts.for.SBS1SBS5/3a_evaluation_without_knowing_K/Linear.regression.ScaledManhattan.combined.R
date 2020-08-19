
## Require tidyverse for using "%>%"
require(tidyverse)
require(lattice)
require(dplyr)
require(ggplot2)
require(ggpubr)

setwd("../../practice/3_Signature_Challenge/3.2_SBS1-SBS5-Correlation_Project/new_results/FinalExtrAttrSummary")


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

if(!dir.exists(paste0("../../TrendDiag/ExtrAttr/")))
  dir.create(paste0("../../TrendDiag/ExtrAttr/"),recursive = T)


## Load final summary for all approaches.
load("FinalSummary.RDa")

## Calculate means for separate attribution performances
## (scaled Manhattan distances)
## conditioning on toolNames, Rsqs (datasetGroups)
## and SBS1:SBS5 ratios (datasetSubGroups)
summaries <- list()


for(index in c("SBS1","SBS5")){

  summaries[[index]] <-
    FinalSummary$FinalAttr[[index]]

  summaries[[index]] <- summaries[[index]] %>%
    dplyr::select(value,toolName,datasetGroup,datasetSubGroup) %>%
    dplyr::mutate(datasetGroup = as.numeric(datasetGroup),
                  datasetSubGroup = as.numeric(datasetSubGroup)) %>%
    dplyr::group_by(toolName,datasetGroup,datasetSubGroup) %>%
    dplyr::summarise(Mean = mean(value),.groups = "keep")

}

## Plot trend of means for each computational approach.
ggVsRsq <- list()
ggVsRatio <- list()
for(index in c("TPR","FDR","SBS1","SBS5")){
  ggVsRsq[[index]] <- ggplot2::ggplot(
    summaries[[index]],
    aes(x = datasetGroup)) +
    ## by specifying group aesthetics = toolName,
    ## we draw a line for each toolName.
    ## This is done implicitly when specifying colour = toolName,
    ## yet we still specify it explicitly, as it will become
    ## a disaster if we specify group = datasetSubGroup.
    ggplot2::geom_line(aes(y = Mean, colour = toolName)) +
    ggplot2::labs(x = "Pearson's R squared") +
    ggplot2::facet_wrap(facets = vars(datasetSubGroup))

  ggVsRatio[[index]] <- ggplot2::ggplot(
    summaries[[index]],
    aes(x = datasetSubGroup)) +
    ## by specifying group aesthetics = toolName,
    ## we draw a line for each toolName.
    ## This is done implicitly when specifying colour = toolName,
    ## yet we still specify it explicitly, as it will become
    ## a disaster if we specify group = datasetGroup.
    ggplot2::geom_line(aes(y = Mean, colour = toolName, group = toolName)) +
    ggplot2::labs(x = "SBS1:SBS5 exposure ratio") +
    ggplot2::facet_wrap(facets = vars(datasetGroup))
}

## Plot these summary plots separately, with customized x labels and y labels.
{
pdf("../TrendDiag/ExtrAttr/ScaledManhattanSBS1AgainstRsq.pdf")
  ggVsRsq$SBS1 <- ggVsRsq$SBS1 + ggplot2::labs(y = "Scaled Manhattan distance of SBS1")
plot(ggVsRsq$SBS1)
dev.off()

pdf("../TrendDiag/ExtrAttr/ScaledManhattanSBS5AgainstRsq.pdf")
ggVsRsq$SBS5 <- ggVsRsq$SBS5 + ggplot2::labs(y = "Scaled Manhattan distance of SBS5")
plot(ggVsRsq$SBS5)
dev.off()
}

{
pdf("../TrendDiag/ExtrAttr/ScaledManhattanSBS1AgainstSBS1SBS5Ratio.pdf")
  ggVsRatio$SBS1 <- ggVsRatio$SBS1 + ggplot2::labs(y = "Scaled Manhattan distance of SBS1")
plot(ggVsRatio$SBS1)
dev.off()

pdf("../TrendDiag/ExtrAttr/ScaledManhattanSBS5AgainstSBS1SBS5Ratio.pdf")
ggVsRatio$SBS5 <- ggVsRatio$SBS5 + ggplot2::labs(y = "Scaled Manhattan distance of SBS5")
plot(ggVsRatio$SBS5)
dev.off()
}

## Plot multiple summary plots into one page.
## These PDF files with multiple panels are ready for publication.
{
  pdf("../TrendDiag/ExtrAttr/ScaledManhattanCombined.pdf", width = 10, height = 10)
  ggObj <- ggpubr::ggarrange(
    ggVsRsq$SBS1 + rremove("legend"),
    ggVsRsq$SBS5 + rremove("legend"),
    ggVsRatio$SBS1 + rremove("legend"),
    ggVsRatio$SBS5,
    labels = c("A","B","C","D"),
    font.label = list(size = 14, color = "black", face = "bold", family = "sans"),
    nrow = 2, ncol = 2,
    common.legend = T)
  plot(ggObj)
  dev.off()
}

save.image("../TrendDiag/ExtrAttr/Manhattan.Combined.Rdata")
