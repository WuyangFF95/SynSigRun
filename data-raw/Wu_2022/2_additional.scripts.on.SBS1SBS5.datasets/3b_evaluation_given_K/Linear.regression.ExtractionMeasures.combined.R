
## Require tidyverse for using "%>%"
require(tidyverse)
require(lattice)
require(dplyr)
## Require plotting package
require(ggplot2)
## Require package combining multiple grid plots
require(ggpubr)
## Require package altering axis scales
require(scales)
## Require package altering color
require(RColorBrewer)

setwd("../../practice/3_Signature_Challenge/3.2_SBS1-SBS5-Correlation_Project/new_results/FinalExtrAttrExactSummary")


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

if(!dir.exists(paste0("../TrendDiag/ExtrAttrExact/")))
  dir.create(paste0("../TrendDiag/ExtrAttrExact/"),recursive = T)


## Load final summary for all approaches.
load("FinalSummary.RDa")

## Calculate means for separate extraction performances
## conditioning on toolNames, Rsqs (datasetGroups)
## and SBS1:SBS5 ratios (datasetSubGroups)
summaries <- list()

for(index in c("TPR","PPV")){

  summaries[[index]] <-
    FinalSummary$FinalExtr[[index]]

  summaries[[index]] <- summaries[[index]] %>%
    dplyr::select(value,toolName,datasetGroup,datasetSubGroup) %>%
    dplyr::mutate(datasetGroup = as.numeric(datasetGroup),
                  datasetSubGroup = as.numeric(datasetSubGroup)) %>%
    dplyr::group_by(toolName,datasetGroup,datasetSubGroup) %>%
    dplyr::summarise(Mean = mean(value),.groups = "keep")
}

for(index in c("SBS1","SBS5")){

  summaries[[index]] <-
    FinalSummary$FinalExtr$cosSim[[index]]

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
for(index in c("TPR","PPV","SBS1","SBS5")){

  numTools <- length(unique(summaries[[index]]$toolName))

  ggVsRsq[[index]] <- ggplot2::ggplot(
    summaries[[index]],
    aes(x = datasetGroup)) +
    ## by specifying group aesthetics = toolName,
    ## we draw a line for each toolName.
    ## This is done implicitly when specifying colour = toolName,
    ggplot2::geom_line(aes(y = Mean, colour = toolName, linetype = toolName)) +
    ggplot2::labs(
      x = "SBS1-SBS5 correlation (in Pearson's R squared)",
      ## Change title of legend
      colour = "Computational approaches",
      linetype = "Computational approaches"
    ) +
    ## Change colour scale for lines
    scale_color_manual(values = rainbow(numTools)) +
    ## Make line segment in legend longer.
    ggplot2::theme(
      legend.position = "right",
      legend.key.width = grid::unit(0.5,"inches")
    ) +
    ## Change x axis to log2 scaled.
    ## Show tick marks exactly at possible Pearson's R^2 values.
    ## (0.1, 0.2, 0.3, 0.6)
    ggplot2::scale_x_continuous(
      trans = log2_trans(),
      breaks = c(0.1, 0.2, 0.3, 0.6)
    ) +
    ## Fix limit of y axis to be 0 to 1
    ## for all extraction measures
    ggplot2::scale_y_continuous(limits = c(0,1)) +
    ggplot2::facet_wrap(facets = vars(datasetSubGroup)) +
    ## Add facet title on the top,
    ## implemented by adding a sub-title
    ggplot2::labs(subtitle = "SBS1:SBS5 exposure ratio") +
    ggplot2::theme(
      plot.title = element_blank(),
      plot.subtitle = element_text(hjust = 0.5))

  ggVsRatio[[index]] <- ggplot2::ggplot(
    summaries[[index]],
    aes(x = datasetSubGroup)) +
    ## by specifying group aesthetics = toolName,
    ## we draw a line for each toolName.
    ## This is done implicitly when specifying colour = toolName,
    ggplot2::geom_line(aes(y = Mean, colour = toolName, linetype = toolName)) +
    ggplot2::labs(
      x = "SBS1:SBS5 exposure ratio",
      ## Change title of legend
      colour = "Computational approaches",
      linetype = "Computational approaches"
    ) +
    ## Change colour scale for lines
    scale_color_manual(values = rainbow(numTools)) +
    ## Make line segment in legend longer.
    ggplot2::theme(
      legend.position = "right",
      legend.key.width = grid::unit(0.5,"inches")
    ) +
    ## Change x axis to log10 scaled.
    ## Show tick marks exactly at possible SBS1:SBS5 exposure ratio values.
    ## (0.1, 0.5, 1, 2, 10)
    ggplot2::scale_x_continuous(
      trans = log10_trans(),
      breaks = c(0.1, 0.5, 1, 2, 10)
      ) +
    ## Fix limit of y axis to be 0 to 1
    ## for all extraction measures
    ggplot2::scale_y_continuous(limits = c(0,1)) +
    ggplot2::facet_wrap(facets = vars(datasetGroup)) +
    ## Add facet title on the top,
    ## implemented by adding a sub-title
    ggplot2::labs(subtitle = "SBS1-SBS5 correlation (in Pearson's R squared)") +
    ggplot2::theme(
      plot.title = element_blank(),
      plot.subtitle = element_text(hjust = 0.5))
}

## Plot these summary plots separately, with customized x labels and y labels.
{
pdf("../TrendDiag/ExtrAttrExact/cosSimSBS1AgainstRsq.pdf",10,8)
  ggVsRsq$SBS1 <- ggVsRsq$SBS1 + ggplot2::labs(y = "Cosine similarity to SBS1")
plot(ggVsRsq$SBS1)
dev.off()

pdf("../TrendDiag/ExtrAttrExact/cosSimSBS5AgainstRsq.pdf",10,8)
ggVsRsq$SBS5 <- ggVsRsq$SBS5 + ggplot2::labs(y = "Cosine similarity to SBS5")
plot(ggVsRsq$SBS5)
dev.off()

pdf("../TrendDiag/ExtrAttrExact/TPRAgainstRsq.pdf",10,8)
ggVsRsq$TPR <- ggVsRsq$TPR + ggplot2::labs(y = "True Positive Rate")
plot(ggVsRsq$TPR)
dev.off()

pdf("../TrendDiag/ExtrAttrExact/PPVAgainstRsq.pdf",10,8)
ggVsRsq$PPV <- ggVsRsq$PPV + ggplot2::labs(y = "Positive Predictive Value")
plot(ggVsRsq$PPV)
dev.off()
}

{
pdf("../TrendDiag/ExtrAttrExact/cosSimSBS1AgainstSBS1SBS5Ratio.pdf",10,8)
  ggVsRatio$SBS1 <- ggVsRatio$SBS1 + ggplot2::labs(y = "Cosine similarity to SBS1")
plot(ggVsRatio$SBS1)
dev.off()

pdf("../TrendDiag/ExtrAttrExact/cosSimSBS5AgainstSBS1SBS5Ratio.pdf",10,8)
ggVsRatio$SBS5 <- ggVsRatio$SBS5 + ggplot2::labs(y = "Cosine similarity to SBS5")
plot(ggVsRatio$SBS5)
dev.off()

pdf("../TrendDiag/ExtrAttrExact/TPRAgainstSBS1SBS5Ratio.pdf",10,8)
ggVsRatio$TPR <- ggVsRatio$TPR + ggplot2::labs(y = "True Positive Rate")
plot(ggVsRatio$TPR)
dev.off()

pdf("../TrendDiag/ExtrAttrExact/PPVAgainstSBS1SBS5Ratio.pdf",10,8)
ggVsRatio$PPV <- ggVsRatio$PPV + ggplot2::labs(y = "Positive Predictive Value")
plot(ggVsRatio$PPV)
dev.off()
}

## Plot multiple summary plots into one page.
## These PDF files with multiple panels are ready for publication.
{
  pdf("../TrendDiag/ExtrAttrExact/cosSimCombined.SBS1.pdf", width = 10, height = 15)
  ggObj <- ggpubr::ggarrange(
    ggVsRsq$SBS1 + rremove("legend"),
    ggVsRatio$SBS1,
    labels = c("A","B"),
    font.label = list(size = 14, color = "black", face = "bold", family = "sans"),
    nrow = 2, ncol = 1,
    legend = "right",
    common.legend = T)
  plot(ggObj)
  dev.off()
}

{
  pdf("../TrendDiag/ExtrAttrExact/cosSimCombined.SBS5.pdf", width = 10, height = 15)
  ggObj <- ggpubr::ggarrange(
    ggVsRsq$SBS5 + rremove("legend"),
    ggVsRatio$SBS5,
    labels = c("A","B"),
    font.label = list(size = 14, color = "black", face = "bold", family = "sans"),
    nrow = 2, ncol = 1,
    legend = "right",
    common.legend = T)
  plot(ggObj)
  dev.off()
}

save.image("../TrendDiag/ExtrAttrExact/Extraction.Measures.combined.Rdata")
