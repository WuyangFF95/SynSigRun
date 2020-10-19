

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
## Specify other global variables
dataset.dirs = datasetNames
second.third.level.dirname = "sp.sp/ExtrAttrExact"
out.dir = "./FinalExtrAttrExactSummary"
overwrite = T

## Summarizing attribution Scaled Manhattan distance results.
{
  FinalAttr <- list()
  ## Combine attribution assessment onto multiple sheets.
  ## Each sheet shows Scaled Manhattan distance for one mutational signature.
  for(datasetDir in dataset.dirs){
    thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
    ## Add multiTools <- NULL to please R check
    multiTools <- NULL
    load(paste0(thirdLevelDir,"/multiTools.RDa"))

    gtSigNames <- multiTools$gtSigNames
    sigNums <- length(gtSigNames)

    if(length(FinalAttr) == 0){
      for(gtSigName in gtSigNames) {
        FinalAttr[[gtSigName]] <- data.frame()
      }
    }

    ## Combine Scaled Manhattan distance
    for(gtSigName in gtSigNames){
      FinalAttr[[gtSigName]] <- rbind(
        FinalAttr[[gtSigName]],
        multiTools$ManhattanDist[[gtSigName]])
    }
  }
}


## Plot general png and pdf for attribution Scaled Manhattan distance summary
## Plot a general violin + beeswarm plot for multiple signatures
## in all runs and in all datasets.
{
  plotDFList <- list()

  ## For ground-truth signature,
  ## Create a data.frame integrating results of
  ## all runs and for all datasets
  gtSigNames <- multiTools$gtSigNames
  sigNums <- length(gtSigNames)

  for(gtSigName in gtSigNames){
    plotDFList[[gtSigName]] <- data.frame()
  }

  ## For each dataset, combine the gtSigName values into plotDFList[[gtSigName]]
  for(datasetDir in dataset.dirs){
    thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
    ## Add multiTools <- NULL to please R check
    multiTools <- NULL
    load(paste0(thirdLevelDir,"/multiTools.RDa"))

    for(gtSigName in gtSigNames){
      plotDFList[[gtSigName]] <- rbind(plotDFList[[gtSigName]],multiTools$ManhattanDist[[gtSigName]])
    }
  }

  ## Combine all plotDFList[[gtSigName]] into plotDFList$Combined
  plotDFList$combined <- data.frame()
  for(gtSigName in gtSigNames){
    plotDFList$combined <- rbind(plotDFList$combined,plotDFList[[gtSigName]])
  }

  ## Convert plotDFList$combined$datasetGroup and
  ## Let their levels follow gtools::mixedsort() fashion
  ## So that the order of the facet labels will be more reasonable for readers.
  plotDFList$combined$datasetGroup <- factor(
    plotDFList$combined$datasetGroup,
    levels = gtools::mixedsort(unique(plotDFList$combined$datasetGroup)))

  if(!is.null(multiTools$datasetSubGroup) & !is.null(multiTools$datasetSubGroupName)) {
    plotDFList$combined$datasetSubGroup <- factor(
      plotDFList$combined$datasetSubGroup,
      levels = gtools::mixedsort(unique(plotDFList$combined$datasetSubGroup)))
  }

  ggplotList <- list()
  ## Plot a multi-facet ggplot for all gtSigNames and all runs.
  {
    ## Generate a ggplot object based on plotDFList$combined
    ggplotList$general <- ggplot2::ggplot(
      plotDFList$combined,
      ggplot2::aes(x = .data$toolName, y = .data$value)) +
      ## Draw geom_violin and geom_quasirandom
      ggplot2::geom_violin(
        ## Change filling color to white
        fill = "#FFFFFF",
        #ggplot2::aes(fill = gtSigName),
        ## Maximize the violin plot width
        scale = "width",
        ## Make bandwidth larger
        #position = "dodge",
        #width = 1.2
        ## Hide outliers
        #outlier.shape = NA
      ) +
      #ggbeeswarm::geom_quasirandom(
      #  groupOnX = TRUE, size = 0.3
      #  ,ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
      #) +
      ## Show median of the Scaled Manhattan distance distribution
      ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
      ## Show mean of the extraction meaasure distribution, as a blue diamond.
      ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
      ## Add title for general violin + beeswarm plot
      ggplot2::ggtitle(label = "Scaled Manhattan distance between inferred and grond-truth exposures",
                       subtitle = "for all computational approaches, ratios and correlation values.") +
      ## Change axis titles
      ggplot2::labs(x = "Computational approach",
                    y = "Scaled Manhattan distance") +
      ## Rotate the names of tools,
      ## move axis.text.x right below the tick marks
      ## and remove legends
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        ## Rotate the axis.text.x (names of tools),
        angle = 90,
        ## move axis.text.x right below the tick marks
        hjust = 1, vjust = 0.5),
        ## remove legends.
        legend.position = "none") +
      ## Split the plot into multiple facets,
      ## according to different gtSigNames
      ggplot2::facet_wrap(
        ggplot2::vars(gtSigName),
        ## Force facet_wrap to have 2 columns
        ncol = 2,
        scales = "free",
        ## Let facets be plotted vertically
        dir = "v"
      ) +
      ## Restrict the decimal numbers of values of measures to be 2
      ggplot2::scale_y_continuous(
        ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
        #limits = c(0, max(plotDFList$combined$value)),
        limits = c(0,5),
        labels =function(x) sprintf("%d", x))
  }
  ## Plot a multi-facet ggplot,
  ## facets are separated by gtSigNames and datasetGroup
  ## (in example, it refers to slope.)
  if(!is.null(multiTools$datasetSubGroup) & !is.null(multiTools$datasetSubGroupName)) {
    bys <- c("datasetGroup","datasetSubGroup")
  } else {
    bys <- c("datasetGroup")
  }

  for(by in bys)  {

    ## The value of "datasetGroupName" or "datasetSubGroupName"
    ## which is the caption of "datasetGroup"
    byCaption <- eval(parse(
      text = paste0("multiTools$",by,"Name")))


    ## Generate a ggplot object based on plotDFList$combined
    ggplotList[[by]] <- ggplot2::ggplot(
      plotDFList$combined,
      ggplot2::aes(x = .data$toolName, y = .data$value))
    ## Draw geom_violin and geom_quasirandom
    ggplotList[[by]] <- ggplotList[[by]] +
      ggplot2::geom_violin(
        ## Change filling color to white
        fill = "#FFFFFF",
        #ggplot2::aes(fill = gtSigName),
        ## Maximize the violin plot width
        scale = "width"
        #,
        ## Make bandwidth larger
        #position = "dodge",
        #width = 1.2
        ## Hide outliers
        #outlier.shape = NA
      ) +
      #ggbeeswarm::geom_quasirandom(
      #  groupOnX = TRUE, size = 0.3
      #  ,ggplot2::aes(color = grDevices::hcl(h = 300,c = 35,l = 60)) ## A purple color, albeit deeper than default hcl colors.
      #) +
      ## Show median of the Scaled Manhattan distance distribution
      ggplot2::stat_summary(fun.y="median", geom="point", shape = 21, fill = "red") +
      ## Show mean of the extraction meaasure distribution, as a blue diamond.
      ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="blue") +
      ## Add title for general violin + beeswarm plot
      ggplot2::ggtitle(
        label = paste0("Scaled Manhattan distance summary plot as a function of "),
        subtitle = paste0("ground-truth signature names and ",byCaption,".")) +
      ## Change axis titles
      ggplot2::labs(x = "Computational approach",
                    y = "Scaled Manhattan distance") +
      ## Rotate the axis.text.x (names of tools),
      ## move axis.text.x right below the tick marks
      ## and remove legends
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        ## Rotate the axis.text.x (names of tools),
        angle = 90,
        ## move axis.text.x right below the tick marks
        hjust = 1, vjust = 0.5),
        ## remove legends.
        legend.position = "none") +
      ## Split the plot into multiple facets,
      ## according to different gtSigNames
      ggplot2::facet_grid(rows =  ggplot2::vars(gtSigName),
                          cols = eval(parse(text = paste0("ggplot2::vars(",by,")"))),
                          scales = "free") +
      ## Restrict the decimal numbers of values of measures to be 2
      ggplot2::scale_y_continuous(
        ## For scaled Manhattan distance, set ylim from 0 to the maximum of Manhattan distance value
        #limits = c(0, max(plotDFList$combined$value)),
        limits = c(0,5),
        labels =function(x) sprintf("%d", x))
  }

  ## Plot violin + beeswarm plots in pdf format
  grDevices::pdf(paste0(out.dir,"/Manhattan.Dist.violins.ylim.restricted.pdf"), pointsize = 1)
  for(by in names(ggplotList)){
    print(ggplotList[[by]])
  }
  grDevices::dev.off()
}
