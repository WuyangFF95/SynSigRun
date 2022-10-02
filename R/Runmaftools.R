#' Install maftools from Bioconductor,
#' and its dependent package, NMF.
#'
#' @keywords internal
Installmaftools <- function(){
  message("Installing maftools from Bioconductor...")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    utils::install.packages("BiocManager")
  BiocManager::install("maftools")
  if (!requireNamespace("NMF", quietly = TRUE))
    utils::install.packages("NMF")
}


#' Run maftools extraction ONLY on a spectra catalog file
#'
#' WARNING: maftools can only do signature extraction!
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
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
#' @param K.exact,K.range \code{K.exact} is the exact value for
#' the number of signatures active in spectra (K).
#' Specify \code{K.exact} if you know exactly how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell maftools to search the best
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
#' NOTE: Unlike other NMF-based packages, parameter \code{nrun.extract} is
#' hard-coded as 1.
#'
#' @param pConstant A small positive value (a.k.a. pseudocount)
#' to add to every entry in the \code{input.catalog}.
#' Specify a value ONLY if an "non-conformable arrays error"
#' is raised.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#'
#' @param overwrite If TRUE, overwrite existing output.
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
           out.dir,
           CPU.cores = NULL,
           K.exact = NULL,
           K.range = NULL,
           nrun.est.K = 10,
           pConstant = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    # Check whether ONLY ONE of K.exact or K.range is specified.
    bool1 <- is.numeric(K.exact) & is.null(K.range)
    bool2 <- is.null(K.exact) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    # Install maftools, if failed to be loaded
    if (!requireNamespace("maftools", quietly = TRUE)) {
      Installmaftools()
    }

    # Set seed as 123456
    set.seed(123456)
    seedInUse <- .Random.seed  # Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() # Save the random number generator (RNG) used


    # Read in spectra data from input.catalog file
    # spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                  strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]
    # convSpectra: convert the ICAMS-formatted spectra catalog
    # into a matrix which maftools::estimateSignatures() and
    # matools::extractSignatures() accepts:
    # 1. Remove the catalog related attributes in convSpectra
    # 2. Transpose the catalog
    # The matools spectra can also be created by
    # matools::trinucleotideMatrix().
    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    # pConstant is added to convSpectra in
    # functions maftools::estimateSignatures() and
    # maftools::extractSignatures().
    convSpectra <- list("nmf_matrix" = t(convSpectra))

    # Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    # CPU.cores specifies number of CPU cores to use.
    # If CPU.cores is not specified, CPU.cores will
    # be equal to the minimum of 30 or (total cores)/2
    if(is.null(CPU.cores)){
      CPU.cores = min(30,(parallel::detectCores())/2)
    } else {
      stopifnot(is.numeric(CPU.cores))
    }

    # Before running NMF packge,
    # Load it explicitly to prevent errors.
    requireNamespace("NMF")

    # Run NMF using ICAMS-formatted spectra catalog
    # Determine the best number of signatures (K.best).
    # If K.exact is provided, use K.exact as the K.best.
    # If K.range is provided, determine K.best by doing raw extraction
    # with maftools::estimateSignature
    if(bool1){
      K.best <- K.exact
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){
      # Change K.range to a full vector
      # if it is already a full vector, just keep it.
      K.range <- seq.int(min(K.range),max(K.range))

      res <- maftools::estimateSignatures(
        mat = convSpectra,
        nMin = K.range[1],
        nTry = K.range[2],
        nrun = nrun.est.K,
        parallel = CPU.cores,
        pConstant = pConstant)

      gof_nmf <- res$nmfObj

      # Choose the best signature number (K.best) active in the spectra
      # catalog (input.catalog).
      ##
      # According to paper "A flexible R package for nonnegative matrix factorization"
      # (Gaujoux & Seoighe, 2010), the most common approach to choose number of
      # signature (K, a.k.a. rank in this paper) is to choose the smallest K for which
      # cophenetic correlation coefficient starts decreasing.
      for(current.K in K.range)
      {
        # Stop the cycle if current.K reaches the maximum.
        # At max(K.range), next.summary becomes meaningless.
        if(current.K == max(K.range))
          break

        current.summary <- NMF::summary(gof_nmf$fit[[as.character(current.K)]])
        current.cophenetic.coefficient <- current.summary["cophenetic"]

        next.summary <- NMF::summary(gof_nmf$fit[[as.character(current.K+1)]])
        next.cophenetic.coefficient <- next.summary["cophenetic"]

        if(current.cophenetic.coefficient > next.cophenetic.coefficient)
          break
      }
      K.best <- current.K # Choose K.best as the smallest current.K whose cophenetic
      # is greater than cophenetic from (current.K+1).
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))
    }

    # Extract signatures using maftools::extractSignatures
    sigs_nmf <- maftools::extractSignatures(
      convSpectra,
      n = K.best,
      parallel = CPU.cores,
      pConstant = pConstant)

    # Generates a list contain extracted signatures
    # sigs_nmf$signatures is already normalized.
    extractedSignatures <- sigs_nmf$signatures
    # Add signature names for signature matrix extractedSignatures
    colnames(extractedSignatures) <-
      paste("maftools",1:ncol(extractedSignatures),sep=".")
    extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                             region = "unknown",
                                             catalog.type = "counts.signature")
    # Output extracted signatures in ICAMS format
    ICAMS::WriteCatalog(extractedSignatures,
                        paste0(out.dir,"/extracted.signatures.csv"))


    # Derive exposure count attribution results.


    # exposure attributions (in percentage, normalized)
    exposureRaw <- (sigs_nmf$contributions)
    rownames(exposureRaw) <- paste("maftools",1:nrow(exposureRaw),sep=".")
    # convert exposure ratio to exposure counts
    exposureCounts <- exposureRaw
    for(sample in colnames(exposureCounts)){
      exposureCounts[,sample] <- exposureRaw[,sample] * sum(spectra[,sample])
    }


    # Write exposure counts in ICAMS and SynSig format.
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
