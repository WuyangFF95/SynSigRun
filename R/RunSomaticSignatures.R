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



#' Run SomaticSignatures extraction and attribution on a spectra catalog file
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
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
#' of length 2 which tell SomaticSignatures to search the best
#' signature number active in spectra, K, in this range of Ks.
#' Specify \code{K.range} if you don't know how many signatures
#' are active in the \code{input.catalog}.
#'
#' WARNING: You must specify only one of \code{K.exact} or \code{K.range}!
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
#' @return The attributed exposure of \code{SomaticSignatures}, invisibly.
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
           seedNumber = 1,
           K.exact = NULL,
           K.range = NULL,
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

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## convSpectra: convert the ICAMS-formatted spectra catalog
    ## into a matrix which SomaticSignatures accepts:
    ## 1. Remove the catalog related attributes in convSpectra
    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    dimnames(convSpectra) <- dimnames(spectra)
    sample.number <- dim(spectra)[2]

    ## Signature extraction
    ## TODO(Wuyang)
    extractedSignatures <- NULL

    ## Write extracted signatures
    ICAMS::WriteCatalog(extractedSignatures,
                        paste0(out.dir,"/extracted.signatures.csv"))


    ## Derive exposure count attribution results.
    ## TODO(Wuyang)
    extractionObject <- NULL
    exposureCounts <- extractionObject$Ehat ## Unnormalized exposures
    ## Temporary workout
    K.best <- NULL
    rownames(exposureCounts) <- paste("SomaticSignatures",seq(1,K.best),sep = ".") ## Assign row names of exposure matrix as names of signatures
    colnames(exposureCounts) <- colnames(spectra) ## Assign column names of exposure matrix as names of tumors
    ## Normalize the attributed counts so that each column represents exposure of a signature
    for(ii in 1:ncol(exposureCounts)) {
      exposureCounts[,ii] <- exposureCounts[,ii] / sum(exposureCounts[,ii])
      exposureCounts[,ii] <- exposureCounts[,ii] * colSums(spectra)[ii]
    }
    ## Save exposure attribution results
    SynSigGen::WriteExposure(exposureCounts,
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
