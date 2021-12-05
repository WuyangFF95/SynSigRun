#' Install SparseSignatures from Bioconductor
#'
#' @keywords internal
InstallSparseSignatures <- function(){
  message("Installing SparseSignatures from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    utils::install.packages("BiocManager")
  BiocManager::install("SparseSignatures")
}



#' Run SparseSignatures extraction and attribution on a spectra catalog file
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run SparseSignatures. Setting seed can make the
#' attribution of SparseSignatures repeatable.
#' Default: 1.
#'
#' @param K.exact,K.range \code{K.exact} is the exact value for
#' the number of signatures active in spectra (K).
#' Specify \code{K.exact} if you know exactly how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell SparseSignatures to search the best
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
#' @return The inferred exposure of \code{SparseSignatures}, invisibly.
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

RunSparseSignatures <-
  function(input.catalog,
           out.dir,
           seedNumber = 1,
           K.exact = NULL,
           K.range = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    # Check whether ONLY ONE of K.exact or K.range is specified.
    bool1 <- is.numeric(K.exact) & is.null(K.range)
    bool2 <- is.null(K.exact) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    # Install SparseSignatures, if not found in library
    if ("SparseSignatures" %in% rownames(utils::installed.packages()) == FALSE)
      InstallSparseSignatures()


    # Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  # Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() # Save the random number generator (RNG) used


    # Read in spectra data from input.catalog file
    # spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]
    # convSpectra: convert the ICAMS-formatted spectra catalog
    # into a matrix which SparseSignatures accepts:
    # 1. Remove the catalog related attributes in convSpectra
    # 2. Transpose the catalog
    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    dimnames(convSpectra) <- dimnames(spectra)
    sample.number <- dim(spectra)[2]
    convSpectra <- t(convSpectra)

    # Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }



    # Determine the best number of signatures (K.best).
    # If K.exact is provided, use K.exact as the K.best.
    # If K.range is provided, determine K.best by doing raw extraction.
    if(bool1){
      # Determine hyperparameter (Lambda) when K is given.
      K.best <- K.exact
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
      LassoCVObject <- SparseSignatures::nmf.LassoCV(
        x = convSpectra,
        K = K.best,
        num_processes = NA,
        seed = seedNumber)
      res <- SparseSignatures::as.mean.squared.error(
        SparseSignatures::cv_example)$median
      resBest <- which(res==min(res),arr.ind=TRUE)
      # Record best Lambda number
      Lambda.best <- colnames(res)[resBest[2]]
      Lambda.best <- as.numeric(
        gsub(pattern="_.*",replacement = "",x = Lambda.best))

      print(paste0("The best number of signatures is found.",
               "It equals to: ",Lambda.best))
    }
    if(bool2){
      # Automatically determine best number of signatures (K)
      # and another hyperparameter, Lambda.
      LassoCVObject <- SparseSignatures::nmf.LassoCV(
        x = convSpectra,
        K = seq(K.range[1],K.range[2]),
        num_processes = NA,
        seed = seedNumber)
      res = SparseSignatures::as.mean.squared.error(SparseSignatures::cv_example)$median
      resBest = which(res==min(res),arr.ind=TRUE)
      # Record best number of signatures and Lambda number
      K.best = rownames(res)[resBest[1]]
      K.best = as.numeric(
      gsub(pattern="_.*",replacement = "",x = K.best))
      Lambda.best = colnames(res)[resBest[2]]
      Lambda.best = as.numeric(
        gsub(pattern="_.*",replacement = "",x = Lambda.best))
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",Lambda.best))
    }


    # Output extracted signatures in Duke-NUS format
    # Extract signatures using the best K and best Lambda.
    extractionObject <- SparseSignatures::nmf.LassoK(
      x = convSpectra,
      K = K.best,
      lambda_rate = Lambda.best,
      num_processes = NA,
      seed = seedNumber)
    extractedSignatures <- t(extractionObject$beta)
    rownames(extractedSignatures) <- rownames(spectra)
    colnames(extractedSignatures) <- paste0("SparseSignatures.",seq(1,ncol(extractedSignatures)))
    extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                             region = "unknown",
                                             catalog.type = "counts.signature")
    # Write extracted signatures
    ICAMS::WriteCatalog(extractedSignatures,
                           paste0(out.dir,"/extracted.signatures.csv"))


    # Derive exposure count attribution results.
    exposureCounts <- t(extractionObject$alpha)
    rownames(exposureCounts) <- paste0("SparseSignatures.",seq(1,nrow(exposureCounts)))
    colnames(exposureCounts) <- colnames(spectra) # Assign column names of exposure matrix as names of tumors

    # Save exposure attribution results
    SynSigGen::WriteExposure(exposureCounts,
                  paste0(out.dir,"/inferred.exposures.csv"))


    # Save seeds and session information
    # for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) # Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) # Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) # Save seed in use to a text file

    # Return a list of signatures and exposures
    invisible(list("signature" = extractedSignatures,
                   "exposure" = exposureCounts))
  }
