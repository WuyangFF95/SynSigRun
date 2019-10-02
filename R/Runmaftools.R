#' Install maftools from Bioconductor,
#' and its dependent package, NMF.
#'
#' @keywords internal
Installmaftools <- function(){
  message("Installing maftools from Bioconductor...")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    utils::install.packages("BiocManager")
  BiocManager::install("maftools")
  if ("NMF" %in% rownames(utils::installed.packages()) == FALSE)
    utils::install.packages("NMF")
}


#' Run maftools extraction ONLY on a spectra catalog file
#'
#' WARNING: maftools can only do signature extraction!
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param read.catalog.function Function to read a catalog
#' (can be spectra or signature catalog): it takes a file path as
#' its only argument and returning a catalog as a numeric matrix.
#'
#' @param write.catalog.function Function to write a catalog.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param CPU.cores Number of CPUs to use in running
#' maftools. For a server, 30 cores would be a good
#' choice; while for a PC, you may only choose 2-4 cores.
#' By default (CPU.cores = NULL), the CPU.cores would be equal
#' to \code{(parallel::detectCores())/2}, total number of CPUs
#' divided by 2.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run maftools. Setting seed can make the
#' attribution of maftools repeatable.
#' Default: 1.
#'
#' @param K,K.range \code{K} is the precise value for
#' the number of signatures active in spectra (K).
#' Specify \code{K} if you know precisely how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell maftools to search the best
#' signature number active in spectra, K, in this range of Ks.
#' Specify \code{K.range} if you don't know how many signatures
#' are active in the \code{input.catalog}.
#'
#' WARNING: You must specify only one of \code{K} or \code{K.range}!
#'
#' Default: NULL
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The extracted signatures of \code{maftools}, invisibly.
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

Runmaftools <-
  function(input.catalog,
           read.catalog.function,
           write.catalog.function,
           out.dir,
           CPU.cores = NULL,
           seedNumber = 1,
           K = NULL,
           K.range = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Check whether ONLY ONE of K or K.range is specified.
    bool1 <- is.numeric(K) & is.null(K.range)
    bool2 <- is.null(K) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    ## Install maftools, if not found in library
    if ("maftools" %in% rownames(utils::installed.packages()) == FALSE)
      Installmaftools()


    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- read.catalog.function(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]
    ## convSpectra: convert the ICAMS-formatted spectra catalog
    ## into a matrix which HDP accepts:
    ## 1. Remove the catalog related attributes in convSpectra
    ## 2. Transpose the catalog
    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    convSpectra <- list("nmf_matrix" = t(convSpectra))

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## CPU.cores specifies number of CPU cores to use.
    ## CPU.cores will be capped at 30.
    ## If CPU.cores is not specified, CPU.cores will
    ## be equal to the minimum of 30 or (total cores)/2
    if(is.null(CPU.cores)){
      CPU.cores = min(30,(parallel::detectCores())/2)
    } else {
      stopifnot(is.numeric(CPU.cores))
      if(CPU.cores > 30) CPU.cores = 30
    }

    ## Before running NMF packge,
    ## Load it explicitly to prevent errors.
    requireNamespace("NMF")

    ## Run NMF using ICAMS-formatted spectra catalog
    ## Determine the best number of signatures (K.best).
    ## If K is provided, use K as the K.best.
    ## If K.range is provided, determine K.best by doing raw extraction.
    if(bool1){
      grDevices::pdf(paste0(out.dir,"/maftools.plots.pdf"))
      gof_nmf <- maftools::extractSignatures(mat = convSpectra,
                          n = K,     ## n specifies number of signatures you want to assess
                          parallel  = paste0("p",CPU.cores))
      grDevices::dev.off()
      K.best <- K
      K.range <- K
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){
      grDevices::pdf(paste0(out.dir,"/maftools.plots.pdf"))
      sigs_nmf <- maftools::extractSignatures(convSpectra,
                          nTry = K.range[2],     ## nTry specifies maximal number of signatures you want to assess
                          parallel  = paste0("p",CPU.cores))
      grDevices::dev.off()

      K.best <- ncol(sigs_nmf$signatures) ## extractSignatures() will pick up K.best automatically.
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))
    }


    ## Generates a list contain extracted signatures
    extractedSignatures <- sigs_nmf$signatures   ## normalized signature matrix
    ## Add signature names for signature matrix extractedSignatures
    colnames(extractedSignatures) <-
      paste("maftools",1:ncol(extractedSignatures),sep=".")
    extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                             region = "unknown",
                                             catalog.type = "counts.signature")


    ## Output extracted signatures in ICAMS format
    write.catalog.function(extractedSignatures,
                           paste0(out.dir,"/extracted.signatures.csv"))


    ## Derive exposure count attribution results.


    ## exposure attributions (in mutation counts)
    exposureCounts <- (sigs_nmf$contributions)
    rownames(exposureCounts) <- paste("maftools",1:row(exposureCounts),sep=".")
    ## Write exposure counts in ICAMS and SynSig format.
    WriteExposure(exposureCounts,
                  paste0(out.dir,"/attributed.exposures.csv"))


    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return a list of signatures and exposures
    invisible(list("signature" = extractedSignatures,
                   "exposure" = exposureCounts))
  }
