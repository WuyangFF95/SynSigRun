#' Install SomaticSignatures from Bioconductor,
#' also installing its dependent package, NMF.
#'
#' @keywords internal
InstallSomaticSignatures <- function(){
  message("Installing SomaticSignatures from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    utils::install.packages("BiocManager")
  BiocManager::install("SomaticSignatures")
}



#' Run SomaticSignatures.NMF extraction and attribution on a spectra catalog file
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#'
#' @param CPU.cores Number of CPUs to use in running
#' SomaticSignatures.NMF. For a server, 30 cores would be a good
#' choice; while for a PC, you may only choose 2-4 cores.
#' By default (CPU.cores = NULL), the CPU.cores would be equal
#' to \code{(parallel::detectCores())/2}, total number of CPUs
#' divided by 2.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run SomaticSignatures. Setting seed can make the
#' attribution of SomaticSignatures repeatable.
#' Default: 1.
#'
#' @param K.exact,K.range \code{K.exact} is the exact value for
#' the number of signatures active in spectra (K).
#' Specify \code{K.exact} if you know exactly how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell SomaticSignatures.NMF to search the best
#' signature number active in spectra, K, in this range of Ks.
#' Specify \code{K.range} if you don't know how many signatures
#' are active in the \code{input.catalog}.
#'
#' WARNING: You must specify only one of \code{K.exact} or \code{K.range}!
#'
#' Default: NULL
#'
#' @param nrun.est.K Number of NMF runs for each possible number of signature.
#' This is used in the step to estimate the most plausible number
#' of signatures in input spectra catalog.
#'
#' @param nrun.extract number of NMF runs for extracting signatures and inferring
#' exposures.
#'
#' @param pConstant A small positive value (a.k.a. pseudocount)
#' to add to every entry in the \code{input.catalog}.
#' Specify a value ONLY if an "non-conformable arrays error"
#' is raised.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return A list contains: \itemize{
#' \item $signature extracted signatures,
#' \item $exposure inferred exposures,
#' } of \code{SomaticSignatures.NMF}, invisibly.
#'
#' @details Creates several
#'  files in \code{out.dir}. These are:
#'  TODO(Steve): list the files
#'
#'  TODO(Wuyang)
#'
#' @importFrom utils capture.output
#'
#' @export

RunSomaticSignatures <-
  function(input.catalog,
           out.dir,
           CPU.cores = NULL,
           seedNumber = 1,
           K.exact = NULL,
           K.range = NULL,
           nrun.est.K = 30,
           nrun.extract = 1,
           pConstant = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Check whether ONLY ONE of K.exact or K.range is specified.
    bool1 <- is.numeric(K.exact) & is.null(K.range)
    bool2 <- is.null(K.exact) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    ## Install SomaticSignatures, if not found in library
    if ("SomaticSignatures" %in% rownames(utils::installed.packages()) == FALSE)
      InstallSomaticSignatures()


    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                  strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]
    ## convSpectra: convert the ICAMS-formatted spectra catalog
    ## into a matrix which SomaticSignatures accepts:
    ## 1. Remove the catalog related attributes in convSpectra
    ## 2. Transpose the catalog
    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    dimnames(convSpectra) <- dimnames(spectra)
    sample.number <- dim(spectra)[2]
    ## Add pConstant to convSpectra.
    if(!is.null(pConstant)) convSpectra <- convSpectra + pConstant

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## CPU.cores specifies number of CPU cores to use.
    ## If CPU.cores is not specified, CPU.cores will
    ## be equal to the minimum of 30 or (total cores)/2
    if(is.null(CPU.cores)){
      CPU.cores = min(30,(parallel::detectCores())/2)
    } else {
      stopifnot(is.numeric(CPU.cores))
    }


    ## Before running NMF packge,
    ## Load it explicitly to prevent errors.
    requireNamespace("NMF")

    ## Run NMF using ICAMS-formatted spectra catalog
    ## Determine the best number of signatures (K.best).
    ## If K.exact is provided, use K.exact as the K.best.
    ## If K.range is provided, determine K.best by doing raw extraction.
    if(bool1){
      K.best <- K.exact
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){

      ## SomaticSignatures used approach in Hutchins et al.
      ## (2008) to estimate K.
      ## This requires calculation of second derivative of
      ## RSS at >2 integers, and thus requires at least 4 Ks
      ## to be assessed.
      if(max(K.range) - min(K.range) < 3)
        stop("To calculate second derivative, K.range should span at least 4 integers.\n")

      assess <- SomaticSignatures::assessNumberSignatures(
        convSpectra,
        nSigs = seq.int(min(K.range),max(K.range)),
        decomposition = SomaticSignatures::nmfDecomposition,
        seed = seedNumber,
        nrun = nrun.est.K)
      rownames(assess) <- as.character(assess[,"NumberSignatures"])

      ## Choose K.best as the smallest current.K
      ## which is an inflection point (changing sign of second derivative)
      ## of RSS value.
      ## For discrete function, we used forward formula to calculate
      ## first and second derivatives.

      ## memory is for storing the old second derivative
      memory <- NULL

      for(current.K in seq.int(K.range[1],K.range[2] - 2))
      {
        RSS.K <- assess[as.character(current.K),"RSS"]
        RSS.Kplus1 <- assess[as.character(current.K + 1),"RSS"]
        RSS.Kplus2 <- assess[as.character(current.K + 2),"RSS"]

        second.deriv <- (RSS.Kplus2 + RSS.K - RSS.Kplus1) / 2

        ## Choose the current.K if the second derivative equals to 0.
        if(second.deriv == 0)
          break

        ## Choose the current.K if the current second derivative has
        ## opposite sign.
        if(!is.null(memory)){
          if(memory * second.deriv < 0)
            break
        }
        memory <- second.deriv
      }
      K.best <- current.K
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))
    }
    ## Signature extraction
    res <- SomaticSignatures::identifySignatures(
      convSpectra,
      K.best,
      SomaticSignatures::nmfDecomposition,
      seed = seedNumber,
      nrun = nrun.extract)
    gc()
    gc()
    gc()
    ## un-normalized signature matrix
    sigsRaw <- res@signatures
    colnames(sigsRaw) <-
      paste("SS.NMF",1:ncol(sigsRaw),sep=".")
    extractedSignatures <- apply(sigsRaw,2,function(x) x/sum(x))   ## normalize each signature's sum to 1
    extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                             region = "unknown",
                                             catalog.type = "counts.signature")
    ## Output extracted signatures in ICAMS format
    ICAMS::WriteCatalog(extractedSignatures,
                        paste0(out.dir,"/extracted.signatures.csv"))


    ## Derive exposure count attribution results.
    rawExposures <- t(res@samples)
    rownames(rawExposures) <-
      paste("SS.NMF",1:nrow(rawExposures),sep=".")
    ## normalize exposure matrix
    exposureCounts <- apply(rawExposures,2,function(x) x/sum(x))
    ## Make exposureCounts real exposure counts.
    for (sample in seq(1,ncol(exposureCounts))){
      exposureCounts[,sample] <-
        colSums(spectra)[sample] * exposureCounts[,sample]
    }
    ## Write exposure counts in ICAMS and SynSig format.
    SynSigGen::WriteExposure(exposureCounts,
                             paste0(out.dir,"/inferred.exposures.csv"))


    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return a list of signatures and exposures
    invisible(list("signature" = extractedSignatures,
                   "exposure" = exposureCounts))
  }
