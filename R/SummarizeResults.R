CopyWithChecks <- function(from, to.dir, overwrite = FALSE) {
  if (!file.exists(from)) {
    warning("Cannot find", from, "\n\nSkipping\n\n")
  } else {
    copy.res <- file.copy(
      from = from, to = to.dir, overwrite = overwrite)
    if (!copy.res)
      cat("Copy from", from, "to directory", to.dir, "failed\n\n")
  }
}


#' Assess/evaluate one result subfolder from a software package
#'
#' Note: For summarizing sigproextractor or SignatureAnalyzer,
#' users should use sigproextractor(SigProfiler-Python) v0.0.5.43+
#' and SignatureAnalyzer 2018-Apr-18.
#'
#' @param run.dir Lowest level path to result of a run. That is,
#' \code{top.dir}/sp.sp/ExtrAttr/SignatureAnalyzer.results/seed.1/. or
#' \code{top.dir}/sa.sa.96/Attr/YAPSA.results/seed.691/
#' Here, \code{top.dir} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory within the \code{run.dir}
#' which stores the software output.
#'
#' @param ground.truth.exposure.dir Folder which stores ground-truth exposures.
#' Usually, it refers to \code{sub.dir}, i.e. \code{run.dir}/../../../
#'
#' @param extracted.sigs.path Path to extracted sigs file, e.g.
#' \code{<run.dir>/SBS96/Selected_Solution/De_Novo_Solution/signatures.PCAWG.format.csv}.
#'
#' @param attributed.exp.path Path to attributed exposures file.
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @param summary.folder.name The name of the folder containing summary results.
#' Usually, it equals to "summary".
#'
#' @keywords internal
#'
#' @importFrom utils capture.output sessionInfo

SummarizeSigOneSubdir <-
  function(run.dir,
           ground.truth.exposure.dir,
           extracted.sigs.path,
           attributed.exp.path = NULL,
           # TODO(Steve): copy this to the summary and do analysis on how much
           # extracted signature contributes to exposures.
           # read.extracted.sigs.fn = NULL,
           # read.ground.truth.sigs.fn = NULL,
           # write.cat.fn = NULL,
           # plot.pdf.fn = NULL,
           overwrite = FALSE,
           summary.folder.name = "summary") {
    ## Specify default catalog treatment functions
    # if(is.null(read.extracted.sigs.fn)) read.extracted.sigs.fn = ICAMS::ReadCatalog
    # if(is.null(read.ground.truth.sigs.fn)) read.ground.truth.sigs.fn = ICAMS::ReadCatalog
    # if(is.null(write.cat.fn)) write.cat.fn = ICAMS::WriteCatalog
    # if(is.null(plot.pdf.fn)) plot.pdf.fn = ICAMS::PlotCatalogToPdf


    ## Output path - path to dump the ReadAndAnalyzeSigs() results
    outputPath <- paste0(run.dir, "/", summary.folder.name)

    ## Analyze signature extraction similarity
    sigAnalysis <-
      ReadAndAnalyzeSigs(
        extracted.sigs = extracted.sigs.path,
        ground.truth.sigs =
          paste0(ground.truth.exposure.dir,"/ground.truth.syn.sigs.csv"),
        ground.truth.exposures =
          paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv") # ,
        # read.extracted.sigs.fn = read.ground.truth.sigs.fn,
        # read.ground.truth.sigs.fn = read.ground.truth.sigs.fn
        )

    if (dir.exists(outputPath)) {
      if (!overwrite) stop(outputPath, " already exists")
    }
    suppressWarnings(dir.create(outputPath))

    # Copies ground.truth exposures from second.level.dir
    # to outputPath == run.dir/<summary.folder.name>.
    CopyWithChecks(
      from = paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv"),
      to.dir = outputPath,
      overwrite = TRUE)

    # Writes bi-directional matching and cos.sim calculation
    write.csv(sigAnalysis$match1,
              file = paste(outputPath,"match1.csv",sep = "/"))
    write.csv(sigAnalysis$match2,
              file = paste(outputPath,"match2.csv",sep = "/"))

    # Writes ground truth and extracted signatures
    # write.cat.fn(
    ICAMS::WriteCatalog(
      sigAnalysis$gt.sigs,
      paste(outputPath,"ground.truth.sigs.csv",sep = "/"),
      )
    # write.cat.fn
    ICAMS::WriteCatalog(
      sigAnalysis$ex.sigs,
      paste(outputPath,"extracted.sigs.csv",sep = "/"))

    # Dumps other outputs into "other.results.txt"
    capture.output(
      cat("Average cosine similarity\n"),
      sigAnalysis$avg,
      cat("\nNumber of ground-truth signatures\n"),
      ncol(sigAnalysis$gt.sigs),
      cat("\nNumber of extracted signatures\n"),
      ncol(sigAnalysis$ex.sigs),
      cat("\nsigAnalysis$extracted.with.no.best.match\n"),
      sigAnalysis$extracted.with.no.best.match,
      cat("\nsigAnalysis$ground.truth.with.no.best.match\n"),
      sigAnalysis$ground.truth.with.no.best.match,
      file = paste0(outputPath,"/other.results.txt"))

    ## Obsolete
    if(FALSE){
      #if (class(plot.pdf.fn) == "function") {
      # Output ground-truth sigs to a PDF file
      plot.pdf.fn(sigAnalysis$gt.sigs,
                  paste0(outputPath,"/ground.truth.sigs.pdf"))

      # Output extracted sigs to a PDF file
      plot.pdf.fn(sigAnalysis$ex.sigs,
                  paste0(outputPath,"/extracted.sigs.pdf"))
    }

    ## Plot signatures as "counts.signatures" typed catalog
    # Output ground-truth sigs to a PDF file
    ICAMS::PlotCatalogToPdf(sigAnalysis$gt.sigs,
                paste0(outputPath,"/ground.truth.sigs.pdf"))
    # Output extracted sigs to a PDF file
    ICAMS::PlotCatalogToPdf(sigAnalysis$ex.sigs,
                paste0(outputPath,"/extracted.sigs.pdf"))


    ## Analyze exposure attribution
    # To be compatible with PCAWG project which only studies
    # signature extraction not exposure attribution,
    # errors will not be thrown if !is.null(attributed.exp.path) == F.
    # Here we shouldn't use "exists("attritbuted.exp.path")" because
    # attributed.exp.path is defaulted to be NULL, but is always defined
    # therefore exists.
    if(!is.null(attributed.exp.path)) {

      if(file.exists(attributed.exp.path)) {
        exposureDiff <- ReadAndAnalyzeExposures(
          extracted.sigs = extracted.sigs.path,
          ground.truth.sigs =
            paste0(ground.truth.exposure.dir,"/ground.truth.syn.sigs.csv"),
          attributed.exp.path = attributed.exp.path,
          ground.truth.exposures =
            paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv"))

        # Write results of exposure attribution analysis
        write.csv(exposureDiff,
                  file = paste0(outputPath,"/exposureDifference.csv"),
                  quote = T)

        # Copy attributed exposures to summary folder.
        CopyWithChecks(attributed.exp.path,
                       paste0(outputPath,"/attributed.exposures.csv"),
                       overwrite = overwrite)
      }
      else {
        warning("Cannot find", attributed.exp.path, "\n\nSkipping\n\n")
      }
    }

    ## Log of system time and session info
    capture.output(Sys.time(), sessionInfo(),
                   file = paste0(outputPath,"/log.txt"))

    ## Save Signature extraction summary into RDa file,
    ## for reuse in SummarizeMultiRuns().
    save(sigAnalysis,
         file = paste0(outputPath,"/sigAnalysis.RDa"))
    ## Save exposure attribution summary into RDa file,
    ## for reuse in SummarizeMultiRuns().
    save(exposureDiff,
         file = paste0(outputPath,"/exposureDiff.RDa"))

    invisible(sigAnalysis) # So we have something to check in tests
  }


#' Assess/evaluate multiple summarized runs for one dataset from one software package.
#'
#' Summarize results from each software tool in \code{tool.dir}/\code{run.names}
#' (generated by running a software tool),
#' combine them into \code{tool.dir}.
#'
#' @param datasetName Name of the dataset. (e.g. "S.0.1.Rsq.0.1").
#' Usually, it is has the same name as \code{basename(top.dir)}.
#'
#' @param toolName Name of software package. (e.g. "sigproextractor")
#'
#' @param tool.dir Fourth level path from the \code{top.dir}. Expected to have
#' multiple runs with different names (e.g. "seed.1")
#' That is,
#' \code{top.dir}/sp.sp/ExtrAttr/sa.results/. or
#' \code{top.dir}/sa.sa.96/Attr/deconstructSigs.results/
#'
#' Here, \code{top.dir} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory within the \code{tool.dir}
#' which stores the software output.
#'
#' @param run.names A character vector records the list of \code{run.dir},
#' or fifth level directories from the dataset top-level folder.
#' E.g., c("seed.1","seed.691")
#'
#' @return A list contain c(\code{mean},\code{sd}) of multiple runs:
#' Cosine similarity
#' True Positives(TP): Ground-truth signatures which are active in
#' the spectra, and extracted.
#' False Negatives(FN): Ground-truth signatures not extracted.
#' False Positives(FP): Signatures wrongly extracted, not resembling
#' any ground-truth signatures.
#' True Positive Rate (TPR, Sensitivity): TP / (TP + FN)
#' False Discovery Rate (FDR): FP / (FP + TP)
#'
#' @details Also writes multiple files into folder \code{tool.dir}:
#'
#'
#' @importFrom utils capture.output sessionInfo
#'
#' @export
SummarizeMultiRuns <-
  function(
    datasetName,
    toolName,
    tool.dir,
    run.names){

    ## Indexes for signature extraction in multiple runs
    cosSim <- numeric(0)
    truePos <- numeric(0)
    falsePos <- numeric(0)
    falseNeg <- numeric(0)
    TPR <- numeric(0)
    FDR <- numeric(0)


    for(runName in run.names){
      ## Load directories
      runDir <- paste0(tool.dir,"/",runName)
      summaryDir <- paste0(runDir,"/summary")
      sigAnalysisFile <- paste0(summaryDir,"/sigAnalysis.RDa")
      load(file = sigAnalysisFile)

      cosSim <- c(cosSim,sigAnalysis$avg)

      gtSigsNames <- rownames(sigAnalysis$match2)
      falseNegNames <- sigAnalysis$ground.truth.with.no.best.match
      falsePosNames <- sigAnalysis$extracted.with.no.best.match
      truePosNames <- setdiff(gtSigsNames,falseNegNames)

      falseNeg <- c(falseNeg,length(falseNegNames))
      falsePos <- c(falsePos,length(falsePosNames))
      truePos <- c(truePos, length(truePosNames))

      currentTPR <- length(truePosNames) / length(gtSigsNames)
      currentFDR <- length(falsePosNames) / (length(truePosNames) + length(falsePosNames))
      TPR <- c(TPR, currentTPR)
      FDR <- c(FDR, currentFDR)
    }

  names(cosSim) <- run.names
  names(falseNeg) <- run.names
  names(falsePos) <- run.names
  names(truePos) <- run.names
  names(TPR) <- run.names
  names(FDR) <- run.names

  multiRun <- list()
  ## Save name of the software package and the dataset.
  multiRun$datasetName <- datasetName
  multiRun$toolName <- toolName
  ## Save extraction indexes on multiple runs
  multiRun$cosSim <- cosSim
  multiRun$falseNeg <- falseNeg
  multiRun$falsePos <- falsePos
  multiRun$truePos <- truePos
  multiRun$TPR <- TPR
  multiRun$FDR <- FDR



  ## Calculate mean and SD for indexes of signature extraction
  multiRun$meanSD <- matrix(nrow = 6, ncol = 2)
  indexes <- c("cosSim","falseNeg","falsePos",
                   "truePos","TPR","FDR")
  rownames(multiRun$meanSD) <- indexes
  colnames(multiRun$meanSD) <- c("mean","stdev")
  for(index in indexes){
    currentMean <- mean(multiRun[[index]])
    currentStdev <- stats::sd(multiRun[[index]])
    multiRun$meanSD[index,] <- c(currentMean, currentStdev)
  }

  ## Calculate fivenums for signature extraction
  multiRun$fivenum <- matrix(nrow = 6, ncol = 5)
  indexes <- c("cosSim","falseNeg","falsePos",
                   "truePos","TPR","FDR")
  rownames(multiRun$fivenum) <- indexes
  colnames(multiRun$fivenum) <- c("min","lower-hinge","median","upperhinge","max")
  for(index in indexes){
    currentFiveNum <- stats::fivenum(multiRun[[index]])
    multiRun$fivenum[index,] <- currentFiveNum
  }

  ## Plot boxplot for signature extraction
  grDevices::pdf(paste0(tool.dir,"/boxplot.extraction.indexes.pdf"))
  indexes <- c("cosSim","falseNeg","falsePos",
                   "truePos","TPR","FDR")
  titles <- c("Average cosine similarity",
              "False negatives",
              "False positives",
              "True positives",
              "True Positive Rate (sensitivity)",
              "False Discovery Rate (FDR)")
  subtitles <- c("","Number of ground-truth signatures not extracted",
                 "Number of signatures extracted, but different from ground-truth signatures",
                 "Number of ground-truth signatures extracted",
                 "True Positives / (True Positives + False Negatives)",
                 "False Positives / (True Positives + False Positives)")
  for(indexNum in seq(1,length(indexes))){
    graphics::boxplot(
      multiRun[[ indexes[indexNum] ]],
      main = titles[indexNum],
      sub = subtitles[indexNum])
  }
  grDevices::dev.off()


  ## Indexes for exposure attribution in multiple runs
  ManhattanDist <- matrix(nrow = length(gtSigsNames), ncol = length(run.names))
  rownames(ManhattanDist) <- gtSigsNames
  colnames(ManhattanDist) <- run.names
  for(runName in run.names){
    runDir <- paste0(tool.dir,"/",runName)
    summaryDir <- paste0(runDir,"/summary")
    exposureDiffFile <- paste0(summaryDir,"/exposureDiff.RDa")
    load(file = exposureDiffFile)
    ManhattanDist[gtSigsNames,runName] <- exposureDiff[gtSigsNames,"Manhattan.distance"]
  }
  multiRun$ManhattanDist <- ManhattanDist

  ## Calculate mean and SD for indexes of exposure attribution
  meanSDMD <- matrix(nrow = length(gtSigsNames), ncol = 2)
  rownames(meanSDMD) <- gtSigsNames
  colnames(meanSDMD) <- c("mean","stdev")
  for(sig in gtSigsNames){
    meanSDMD[sig,"mean"] <- mean(ManhattanDist[sig,])
    meanSDMD[sig,"stdev"] <- stats::sd(ManhattanDist[sig,])
  }
  rownames(meanSDMD) <- paste0(rownames(meanSDMD),".Manhattan.Dist")
  multiRun$meanSDMD <- meanSDMD

  ## Calculate fivenums for exposure attribution
  multiRun$fivenumMD <- matrix(nrow = length(gtSigsNames), ncol = 5)
  rownames(multiRun$fivenumMD) <- gtSigsNames
  colnames(multiRun$fivenumMD) <- c("min","lower-hinge","median","upperhinge","max")
  for(sig in gtSigsNames){
    multiRun$fivenumMD[sig,] <- stats::fivenum(ManhattanDist[sig,])
  }

  ## Plot boxplot for exposure attribution
  grDevices::pdf(paste0(tool.dir,"/boxplot.attribution.indexes.pdf"))
  for(sig in gtSigsNames){


      ggplotObj <- ggplot2::ggplot(
        data.frame(value = ManhattanDist[sig,],
                   sigName = sig),
        ggplot2::aes(x = sigName, y = value))
      ggplotObj <- ggplotObj +
        ggplot2::ggtitle(paste0("L1-difference of exposure of signature ",sig))
      ggplotObj <- ggplotObj +
        ggplot2::geom_boxplot(
          notch = FALSE,
          position = ggplot2::position_dodge(0.9))
  }
  grDevices::dev.off()

  ## Save data and results
  save(multiRun,file = paste0(tool.dir,"/multiRun.RDa"))
  write.csv(x = multiRun$ManhattanDist,
            file = paste0(tool.dir,"/ManhattanDist.csv"))
  write.csv(x = multiRun$meanSD,
            file = paste0(tool.dir,"/meanSD.csv"))
  write.csv(x = multiRun$meanSDMD,
            file = paste0(tool.dir,"/meanSD.Manhattan.dist.csv"))
  write.csv(x = multiRun$fivenum,
            file = paste0(tool.dir,"/fivenum.csv"))
  write.csv(x = multiRun$fivenumMD,
            file = paste0(tool.dir,"/fivenum.Manhattan.dist.csv"))
  invisible(multiRun)
}



#' Combine results for a single dataset, from different software tools.
#'
#' Summarize results from each software tool in \code{third.level.dir}/\code{tool.dirnames}
#' (generated by \code{\link{SummarizeMultiRuns}}),
#' combine them into \code{third.level.dir}.
#'
#'
#' @param third.level.dir Third level path distinguishing de-novo extraction
#' + attribution packages from attribution-only packages.
#' Examples:
#' \code{top.dir}/sp.sp/ExtrAttr/
#' \code{top.dir}/sa.sa/Attr/
#'
#' @param toolNames Names of software package. (e.g. "sigproextractor")
#'
#' @param tool.dirnames Third level path from the \code{top.dir}. Expected to have
#' summarized results generated by \code{\link{SummarizeMultiRuns}}.
#' (multiRun.RDa, ManhattanDist.csv, meanSD.csv, meanSD.Manhattan.dist.csv)
#' Examples:
#' \code{"signeR.results"} (Under \code{third.level.dir} "ExtrAttr")
#' \code{"deconstructSigs.results"} (Under \code{third.level.dir} "Attr")
#'
#' Here, \code{top.dir} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory within the \code{tool.names}
#' which stores the software output.
#'
#' @param datasetGroups Numeric or character vector specifying the group
#' each dataset belong to.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider slope as the group:
#' c("slope=0.1","slope=0.5","slope=1","slope=2","slope=5","slope=10")
#' Default: "Default"
#'
#' @param datasetSubGroups Numeric or character vector differentiating
#' datasets within each group.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider Pearson's R^2
#' as the subgroup:
#' c("Rsq=0.1","Rsq=0.2","Rsq=0.3","Rsq=0.6")
#' Default: Names of datasets, which are \code{basename(dataset.dirs)}
#'
#' @return A list contain c(\code{mean},\code{sd}) of multiple runs:
#' Cosine similarity
#' True Positives(TP): Ground-truth signatures which are active in
#' the spectra, and extracted.
#' False Negatives(FN): Ground-truth signatures not extracted.
#' False Positives(FP): Signatures wrongly extracted, not resembling
#' any ground-truth signatures.
#' True Positive Rate (TPR, Sensitivity): TP / (TP + FN)
#' False Discovery Rate (FDR): FP / (FP + TP)
#'
#' @export
#'
#' @importFrom utils capture.output sessionInfo
SummarizeMultiToolsOneDataset <- function(
  third.level.dir,
  toolNames,
  tool.dirnames,
  datasetGroups,
  datasetSubGroups){

  multiTools <- list()
  combMeanSD <- NULL
  combMeanSDMD <- NULL

  for(toolNumber in 1:length(toolNames)){

    toolName <- toolNames[toolNumber]
    toolDirName <- tool.dirnames[toolNumber]
    toolPath <- paste0(third.level.dir,"/",toolDirName)
    load(paste0(toolPath,"/multiRun.RDa"))

    ## Combine multi-runs and multi-tools for each index
    indexes <- c("cosSim","falseNeg","falsePos",
                 "truePos","TPR","FDR")
    for(index in indexes){
      tmp <- data.frame(seed = names(multiRun[[index]]),
                        index = index,
                        value = multiRun[[index]],
                        toolName = toolName,
                        datasetName = datasetName,
                        datasetGroups = datasetGroups,
                        datasetSubGroups = datasetSubGroups,
                        stringsAsFactors = FALSE)

      rownames(tmp) <- NULL

      ## Create a data.frame for each index,
      ## and summarize multi-Run, multiDataset values
      ## for each index.
      if(is.null(multiTools[[index]])){
        multiTools[[index]] <- data.frame()
      }
      multiTools[[index]] <- rbind(multiTools[[index]],tmp)
    }


    ## meanSDMD contains mean and standard deviation
    ## for each extraction index.
    meanSD <- multiRun$meanSD
    colnames(meanSD) <- paste0(toolDirName,".", colnames(meanSD))
    if(is.null(meanSD)){
      combMeanSD <- meanSD
    } else{
      combMeanSD <- cbind(combMeanSD,meanSD)
    }

    ## meanSDMD contains mean and standard deviation
    ## for Manhattan Distance between ground-truth exposures
    ## and attributed exposures for each ground-truth signature
    meanSDMD <- multiRun$meanSDMD
    colnames(meanSDMD) <- paste0(toolDirName,".", colnames(meanSDMD))
    if(is.null(meanSDMD)){
      combMeanSDMD <- meanSDMD
    } else{
      combMeanSDMD <- cbind(combMeanSDMD,meanSDMD)
    }
  }

  multiTools$combMeanSD <- combMeanSD
  multiTools$combMeanSDMD <- combMeanSDMD

  save(multiTools,file = paste0(third.level.dir,"/multiTools.RDa"))
  write.csv(x = multiTools$combMeanSD,
            file = paste0(third.level.dir,"/combined.meanSD.csv"))
  write.csv(x = multiTools$combMeanSDMD,
            file = paste0(third.level.dir,"/combined.meanSD.Manhattan.dist.csv"))
  invisible(multiTools)
}

#' Combine results for multiple datasets, from different software tools.
#'
#' Summarize results from each software tool in \code{third.level.dir}/\code{tool.dirnames}
#' (generated by \code{\link{SummarizeMultiRuns}}),
#' combine them into \code{third.level.dir}.
#'
#' @param dataset.dirs Paths of top-level dataset directories trees you want
#' to investigate.
#' E.g. "./S.0.1.Rsq.0.1"
#'
#' @param second.third.level.dirname Name of the second.level.dir (e.g. "sp.sp")
#' and the third.level.dir (e.g. "ExtrAttr") to be investigated.
#'
#' Examples are: "sp.sp/ExtrAttr", "sa.sa.96/Attr"
#'
#' Note: \code{multiTools.RDa} are expected to be exist under
#' \code{dataset.dirs}/\code{second.third.level.dirname}.
#'
#' @param out.dir Path of the output directory.
#'
#' @param overwrite Whether to overwrite the contents in out.dir if
#' it already exists. (Default: FALSE)
#'
#' @export
#'
SummarizeMultiToolsMultiDatasets <-
  function(dataset.dirs,
           second.third.level.dirname,
           out.dir,
           overwrite = FALSE){

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Summarizing extraction results.
    FinalExtr <- list()
    ## Combine extraction assessment onto 7 sheets:
    for(datasetDir in dataset.dirs){
      thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
      load(paste0(thirdLevelDir,"/multiTools.RDa"))

      indexNum <- nrow(multiTools$combMeanSD)
      if(length(FinalExtr) == 0){
        for(index in seq(1,indexNum)) {
          FinalExtr[[index]] <- data.frame()
        }
        names(FinalExtr) <- rownames(multiTools$combMeanSD)
      }
      current <- list()
      for(index in seq(1,indexNum)){
        current[[index]] <- multiTools$combMeanSD[index,,drop = F]
        rownames(current[[index]]) <- datasetDir
        FinalExtr[[index]] <- rbind(FinalExtr[[index]],current[[index]])
      }
    }
    for(summaryFileName in names(FinalExtr)){
      write.csv(FinalExtr[[summaryFileName]],
                file = paste0(out.dir,"/",summaryFileName,".csv"))
    }

    ## Plot general png and pdf for extraction summary
    ## Plot a general boxplot for multiple indexes
    if(FALSE){
      ggplotObj <- ggplot2::ggplot(
        OneToolSummary[["extraction"]],
        ggplot2::aes(x = toolName, y = value))
      ggplotObj <- ggplotObj +
        ggplot2::geom_boxplot(
          notch = FALSE,
          position = ggplot2::position_dodge(0.9),
          ggplot2::aes(fill = index)) +
        ggplot2::scale_fill_manual(
          values = grDevices::topo.colors(length(indexes)))
      ## Add title for general boxplot
      ggplotObj <- ggplotObj +
        ggplot2::ggtitle(label = "Boxplot for multiple indexes and multiple datasets.")
      print(ggplotObj)
    }


    ## Summarizing attribution results.
    FinalAttr <- list()
    ## Combine attribution assessment onto multiple sheets.
    ## Each sheet shows Manhattan distance for one mutational signature.
    for(datasetDir in dataset.dirs){
      thirdLevelDir <- paste0(datasetDir,"/",second.third.level.dirname)
      load(paste0(thirdLevelDir,"/multiTools.RDa"))

      indexNum <- nrow(multiTools$combMeanSDMD)
      if(length(FinalAttr) == 0){
        for(index in seq(1,indexNum)) {
          FinalAttr[[index]] <- data.frame()
        }
        names(FinalAttr) <- rownames(multiTools$combMeanSDMD)
      }
      current <- list()
      for(index in seq(1,indexNum)){
        current[[index]] <- multiTools$combMeanSDMD[index,,drop = F]
        rownames(current[[index]]) <- datasetDir
        FinalAttr[[index]] <- rbind(FinalAttr[[index]],current[[index]])
      }
    }
    for(summaryFileName in names(FinalAttr)){
      write.csv(FinalAttr[[summaryFileName]],
                file = paste0(out.dir,"/",summaryFileName,".csv"))
    }



    invisible(list(FinalExtr = FinalExtr,
                FinalAttr = FinalAttr))
  }


#' Combine results for multiple datasets, from one software tools.
#'
#' Summarize results from each software tool in \code{third.level.dir}/\code{tool.dirnames}
#' (generated by \code{\link{SummarizeMultiRuns}}),
#' combine them into \code{third.level.dir}.
#'
#' @param dataset.dirs Paths of top-level dataset directories trees you want
#' to investigate.
#' E.g. "./S.0.1.Rsq.0.1"
#'
#' @param datasetGroups Numeric or character vector specifying the group
#' each dataset belong to.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider slope as the group:
#' c("slope=0.1","slope=0.5","slope=1","slope=2","slope=5","slope=10")
#' Default: "Default"
#'
#' @param datasetSubGroups Numeric or character vector differentiating
#' datasets within each group.
#' E.g. For SBS1-SBS5 correlated datasets, we can consider Pearson's R^2
#' as the subgroup:
#' c("Rsq=0.1","Rsq=0.2","Rsq=0.3","Rsq=0.6")
#' Default: Names of datasets, which are \code{basename(dataset.dirs)}
#'
#' @param tool.dirname Name of the second.level.dir (e.g. "sp.sp"),
#' third.level.dir (e.g. "ExtrAttr") and tool.dir
#' (e.g. "sigproextractor.results") to be investigated.
#'
#' One example: "sp.sp/ExtrAttr/sigproextractor.results"
#'
#' Note: this function expects the summary generated by
#' \code{SummarizeSigOneSubdir} under \code{dataset.dirs}/\code{tool.dirname}
#'
#' @param out.dir Path of the output directory.
#'
#' @param overwrite Whether to overwrite the contents in out.dir if
#' it already exists. (Default: FALSE)
#'
#' @export
#'
SummarizeOneToolMultiDatasets <-
  function(dataset.dirs,
           datasetGroups = NULL,
           datasetSubGroups = NULL,
           tool.dirname,
           out.dir,
           overwrite = FALSE){

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exists")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Wrap all datasets into one group, if datasetGroups is NULL.
    ## Re-order the dataset.group for better visualization of
    ## ggplot facets.
    datasetNames <- basename(dataset.dirs)

    if(is.null(datasetGroups))
      datasetGroups <- rep("Default",length(dataset.dirs))
    datasetGroups <- factor(
      datasetGroups,
      levels = gtools::mixedsort(unique(datasetGroups)))
    names(datasetGroups) <- datasetNames

    if(is.null(datasetSubGroups))
      datasetSubGroups <- datasetNames
    datasetSubGroups <- factor(
      datasetSubGroups,
      levels = gtools::mixedsort(unique(datasetSubGroups)))
    names(datasetSubGroups) <- datasetNames


    ## Summarizing extraction results for one software package.
    OneToolSummary <- list()

    ## Combine extraction assessment for multiple datasets
    ## in multiple runs onto 1 sheet:
    OneToolSummary[["extraction"]] <- data.frame()

    for(datasetDir in dataset.dirs){
      thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
      toolName <- strsplit(basename(tool.dirname),".results")[[1]]
      load(paste0(thirdLevelDir,"/multiRun.RDa"))
      indexNum <- nrow(multiRun$meanSD)
      indexes <- rownames(multiRun$meanSD)
      datasetName <- basename(datasetDir)

      for(index in indexes){
        tmp <- data.frame(seed = names(multiRun[[index]]),
                          index = index,
                          value = multiRun[[index]],
                          toolName = toolName,
                          datasetName = datasetName,
                          datasetGroups = datasetGroups[datasetName],
                          datasetSubGroups = datasetSubGroups[datasetName],
                          stringsAsFactors = FALSE)

        rownames(tmp) <- NULL

        ## Create a data.frame for each index,
        ## and summarize multi-Run, multiDataset values
        ## for each index.
        if(is.null(OneToolSummary[[index]])){
          OneToolSummary[[index]] <- data.frame()
        }
        OneToolSummary[[index]] <- rbind(OneToolSummary[[index]],tmp)
      }
    }

    for(index in indexes){
      tmp <- data.frame(OneToolSummary[[index]],
                        stringsAsFactors = FALSE)
      rownames(tmp) <- NULL

      if(nrow(OneToolSummary[["extraction"]]) == 0 |
         ncol(OneToolSummary[["extraction"]]) == 0 |
         is.null(dim(OneToolSummary[["extraction"]])) ) {
        OneToolSummary[["extraction"]] <- tmp
      } else {
        OneToolSummary[["extraction"]] <-
          rbind(OneToolSummary[["extraction"]],tmp)
      }
    }


    ## Draw boxplot for extraction indexes
    ## Designate titles and subtitles for each page
    titles <- c("Average cosine similarity",
                "False negatives",
                "False positives",
                "True positives",
                "True Positive Rate (sensitivity)",
                "False Discovery Rate (FDR)")
    names(titles) <- indexes
    subtitles <- c("","Number of ground-truth signatures not extracted",
                   "Number of signatures extracted, but different from ground-truth signatures",
                   "Number of ground-truth signatures extracted",
                   "True Positives / (True Positives + False Negatives)",
                   "False Positives / (True Positives + False Positives)")
    names(subtitles) <- indexes


    ## Create a list to store ggplot2 boxplot objects
    ggplotObj <- list()
    ## Plot a general boxplot for multiple indexes
    if(FALSE){
      ggplotObj[["general"]] <- ggplot2::ggplot(
        OneToolSummary[["extraction"]],
        ggplot2::aes(x = toolName, y = value))
      ggplotObj[["general"]] <- ggplotObj[["general"]] +
        ggplot2::geom_boxplot(
          notch = FALSE,
          position = ggplot2::position_dodge(0.9),
          ggplot2::aes(fill = index)) +
        ggplot2::scale_fill_manual(
          values = grDevices::topo.colors(length(indexes)))
      ## Add title for general boxplot
      ggplotObj[["general"]] <- ggplotObj[["general"]] +
        ggplot2::ggtitle(label = "Boxplot for multiple indexes and multiple datasets.")
    }
    ## Plot a value~datasetSubGroups boxplot for each index.
    for(index in indexes){
      indexNum <- which(indexes == index)
      ggplotObj[[index]] <- ggplot2::ggplot(
        OneToolSummary[[index]],
        ggplot2::aes(x = datasetSubGroups, y = value))
      ggplotObj[[index]] <- ggplotObj[[index]] +
        ggplot2::geom_boxplot(
          notch = FALSE,
          position = ggplot2::position_dodge(0.9),
          ggplot2::aes(fill = index)) +
        ggplot2::scale_fill_manual(
          values = grDevices::topo.colors(length(indexes))[indexNum]) +
        ggplot2::theme(
          ## Rotate the text for better visualization
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
          ## Remove legend
          legend.position = "none") +
        ggplot2::facet_wrap(~datasetGroups)
      ## Add title for value~datasetSubGroups boxplot
      ggplotObj[[index]] <- ggplotObj[[index]] +
        ggplot2::ggtitle(label = titles[index],subtitle = subtitles[index])

      print(ggplotObj[[index]])
    }


    ## Output average cosine similarity in a png file
    grDevices::png(paste0(out.dir,"/boxplot.onetool.cosine.similarity.png"))
    print(ggplotObj[["cosSim"]])
    grDevices::dev.off()

    ## Output multiple extractionindexes in a pdf file
    grDevices::pdf(paste0(out.dir,"/boxplot.onetool.extraction.indexes.pdf"))
    if(FALSE) print(ggplotObj[["general"]])
    for(index in indexes) print(ggplotObj[[index]])
    grDevices::dev.off()


    ## Write Summary tables
    for(summaryFileName in names(OneToolSummary)){
      write.csv(OneToolSummary[[summaryFileName]],
                file = paste0(out.dir,"/",summaryFileName,".csv"),
                quote = F, row.names = F)
    }

  invisible(NULL)
}
