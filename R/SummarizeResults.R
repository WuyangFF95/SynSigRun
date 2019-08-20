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
#' It defaults to be \code{sub.dir}, i.e. \code{run.dir}/../../
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
           read.extracted.sigs.fn,
           read.ground.truth.sigs.fn,
           write.cat.fn,
           plot.pdf.fn,
           overwrite = FALSE,
           summary.folder.name = "summary") {

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
#' @param overwrite If TRUE overwrite existing directories and files.
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
           run.names,
           overwrite = T){

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


  multiRun <- list()
  multiRun$cosSim <- cosSim
  multiRun$falseNeg <- falseNeg
  multiRun$falsePos <- falsePos
  multiRun$truePos <- truePos
  multiRun$TPR <- TPR
  multiRun$FDR <- FDR

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

  meanSDMD <- matrix(nrow = length(gtSigsNames), ncol = 2)
  rownames(meanSDMD) <- gtSigsNames
  colnames(meanSDMD) <- c("mean","stdev")
  for(sig in gtSigsNames){
    meanSDMD[sig,"mean"] <- mean(ManhattanDist[sig,])
    meanSDMD[sig,"stdev"] <- stats::sd(ManhattanDist[sig,])
  }
  rownames(meanSDMD) <- paste0(rownames(meanSDMD),".Manhattan.Dist")
  multiRun$meanSDMD <- meanSDMD

  save(multiRun,file = paste0(tool.dir,"/multiRun.RDa"))
  write.csv(x = multiRun$ManhattanDist,
            file = paste0(tool.dir,"/ManhattanDist.csv"))
  write.csv(x = multiRun$meanSD,
            file = paste0(tool.dir,"/meanSD.csv"))
  write.csv(x = multiRun$meanSDMD,
            file = paste0(tool.dir,"/meanSD.Manhattan.dist.csv"))
  invisible(multiRun)
}



#' Combine \code{\link{SummarizeMultiRuns}} folders from a software package.
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
#' @param overwrite If TRUE overwrite existing directories and files.
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
SummarizeMultiTools <- function(third.level.dir,
                                tool.dirnames,
                                overwrite = T){

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
