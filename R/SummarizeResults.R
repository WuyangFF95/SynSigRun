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
#' @param read.extracted.sigs.fn Function to read the extracted sigs file.
#' e.g. [ICAMS]\code{ReadCatalog}
#'
#' @param read.ground.truth.sigs.fn Function to read the ground-truth sigs file.
#' e.g. [ICAMS]\code{ReadCatalog}
#'
#' @param plot.pdf.fn If a function, use it to plot PDFs of the ground truth and
#' extracted signatures.
#' e.g. [ICAMS]\code{PlotCatalogToPdf}
#'
#' @param write.cat.fn Function to write a catalog to disk, for example
#' [ICAMS]\code{WriteCatalog}
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
           read.extracted.sigs.fn = NULL,
           read.ground.truth.sigs.fn = NULL,
           write.cat.fn = NULL,
           plot.pdf.fn = NULL,
           overwrite = FALSE,
           summary.folder.name = "summary") {
    ## Specify default catalog treatment functions
    if(is.null(read.extracted.sigs.fn)) read.extracted.sigs.fn = ICAMS::ReadCatalog
    if(is.null(read.ground.truth.sigs.fn)) read.ground.truth.sigs.fn = ICAMS::ReadCatalog
    if(is.null(write.cat.fn)) write.cat.fn = ICAMS::WriteCatalog
    if(is.null(plot.pdf.fn)) plot.pdf.fn = ICAMS::PlotCatalogToPdf


    ## Output path - path to dump the ReadAndAnalyzeSigs() results
    outputPath <- paste0(run.dir, "/", summary.folder.name)

    ## Analyze signature extraction similarity
    sigAnalysis <-
      ReadAndAnalyzeSigs(
        extracted.sigs = extracted.sigs.path,
        ground.truth.sigs =
          paste0(ground.truth.exposure.dir,"/ground.truth.syn.sigs.csv"),
        ground.truth.exposures =
          paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv"),
        read.extracted.sigs.fn = read.ground.truth.sigs.fn,
        read.ground.truth.sigs.fn = read.ground.truth.sigs.fn)

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
    write.cat.fn(
      sigAnalysis$gt.sigs,
      paste(outputPath,"ground.truth.sigs.csv",sep = "/"))
    write.cat.fn(
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

    if (class(plot.pdf.fn) == "function") {
      # Output ground-truth sigs to a PDF file
      plot.pdf.fn(sigAnalysis$gt.sigs,
                  paste0(outputPath,"/ground.truth.sigs.pdf"))

      # Output extracted sigs to a PDF file
      plot.pdf.fn(sigAnalysis$ex.sigs,
                  paste0(outputPath,"/extracted.sigs.pdf"))
    }

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
            paste0(ground.truth.exposure.dir,"/ground.truth.syn.exposures.csv"),
          read.extracted.sigs.fn = read.ground.truth.sigs.fn,
          read.ground.truth.sigs.fn = read.ground.truth.sigs.fn)

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


#' Assess/evaluate multiple summarized runs from a software package.
#'
#' Summarize results from each software tool in \code{tool.dir}/\code{run.names}
#' (generated by running a software tool),
#' combine them into \code{tool.dir}.
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
  function(tool.dir,
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
      falsePos <- c(truePos,length(falsePosNames))
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

  multiRun$cosSim <- cosSim
  multiRun$falseNeg <- falseNeg
  multiRun$falsePos <- falsePos
  multiRun$truePos <- truePos
  multiRun$TPR <- TPR
  multiRun$FDR <- FDR



  ## Calculate mean and SD for indexes of signature extraction
  multiRun$meanSD <- matrix(nrow = 6, ncol = 2)
  toCalculate <- c("cosSim","falseNeg","falsePos",
                   "truePos","TPR","FDR")
  rownames(multiRun$meanSD) <- toCalculate
  colnames(multiRun$meanSD) <- c("mean","stdev")
  for(current in toCalculate){
    currentMean <- mean(multiRun[[current]])
    currentStdev <- stats::sd(multiRun[[current]])
    multiRun$meanSD[current,] <- c(currentMean, currentStdev)
  }

  ## Calculate fivenums for signature extraction
  multiRun$fivenum <- matrix(nrow = 6, ncol = 5)
  toCalculate <- c("cosSim","falseNeg","falsePos",
                   "truePos","TPR","FDR")
  rownames(multiRun$fivenum) <- toCalculate
  colnames(multiRun$fivenum) <- c("min","lower-hinge","median","upperhinge","max")
  for(current in toCalculate){
    currentFiveNum <- fivenum(multiRun[[current]])
    multiRun$fivenum[current,] <- currentFiveNum
  }

  ## Plot boxplot for signature extraction
  pdf(paste0(tool.dir,"/boxplot.extraction.indexes.pdf"))
  toCalculate <- c("cosSim","falseNeg","falsePos",
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
  for(ii in seq(1,length(toCalculate))){
    boxplot(multiRun[[ toCalculate[ii] ]],
            main = titles[ii],
            sub = subtitles[ii])
  }
  dev.off()


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
    multiRun$fivenumMD[sig,] <- fivenum(ManhattanDist[sig,])
  }

  ## Plot boxplot for exposure attribution
  pdf(paste0(tool.dir,"/boxplot.attribution.indexes.pdf"))
  for(sig in gtSigsNames){
    boxplot(ManhattanDist[sig,],
            main = paste0("L1-difference of exposure of signature ",sig))
  }
  dev.off()

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
SummarizeMultiToolsOneDataset <- function(third.level.dir,
                                          tool.dirnames){

  multiTools <- list()
  combMeanSD <- NULL
  combMeanSDMD <- NULL

  for(toolDirName in tool.dirnames){

    toolPath <- paste0(third.level.dir,"/",toolDirName)
    load(paste0(toolPath,"/multiRun.RDa"))

    meanSD <- multiRun$meanSD
    colnames(meanSD) <- paste0(toolDirName,".", colnames(meanSD))
    if(is.null(meanSD)){
      combMeanSD <- meanSD
    } else{
      combMeanSD <- cbind(combMeanSD,meanSD)
    }

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



    return(list(FinalExtr,FinalAttr))
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
           tool.dirname,
           out.dir,
           overwrite = FALSE){

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Summarizing extraction results for one software package.
    OneToolSummary <- list()


    ## Combine extraction assessment for multiple datasets
    ## in multiple runs onto 1 sheet:
    OneToolSummary[["extraction"]] <- data.frame()
    for(index in indexes){
      OneToolSummary[[index]] <- data.frame()
    }

    for(datasetDir in dataset.dirs){
      thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
      load(paste0(thirdLevelDir,"/multiRun.RDa"))
      indexNum <- nrow(multiRun$meanSD)
      indexes <- rownames(multiRun$meanSD)

      for(index in indexes){
        tmp <- multiRun$meanSD[index,,drop = F]
        rownames(tmp) <- datasetDir
        OneToolSummary[[index]] <- rbind(OneToolSummary[[index]],tmp)
      }
    }

    for(index in indexes){
      colnames(OneToolSummary[[index]]) <- paste0(index,colnames(OneToolSummary[[index]]))
      if( dim(OneToolSummary[["extraction"]]) == c(0,0) ) {
        OneToolSummary[["extraction"]] <- OneToolSummary[[index]]
      } else{
        OneToolSummary[["extraction"]] <- cbind(OneToolSummary[["extraction"]], OneToolSummary[[index]])
      }
    }

    ## Draw boxplot for extraction indexes
    pdf(paste0(out.dir,"/boxplot.onetool.extraction.indexes.pdf"))

    toCalculate <- c("cosSim","falseNeg","falsePos",
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


    for(datasetDir in dataset.dirs){

      ## Load multiRun for each dataset.
      thirdLevelDir <- paste0(datasetDir,"/",tool.dirname)
      load(paste0(thirdLevelDir,"/multiRun.RDa"))
      indexNum <- nrow(multiRun$meanSD)
      indexes <- rownames(multiRun$meanSD)

      ## Build a dataframe contain two columns:
      ## Column 1: the value of index (e.g. 1)
      ## Column 2: the name of index (e.g. cosSim, TPR)
      boxplotDF <- data.frame()
      for(index in indexes){
        tmpDF <- data.frame(value = unname(multiRun[[index]]),type = index)
        boxplotDF <- rbind(boxplotDF,tmpDF)
      }
      boxplot(value~type,
              data = boxplotDF)
    }
    dev.off()



    ## Write Summary tables
    for(summaryFileName in names(OneToolSummary)){
      write.csv(OneToolSummary[[summaryFileName]],
                file = paste0(out.dir,"/",summaryFileName,".csv"))
    }


    return(NULL)
}
