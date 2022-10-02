#' Install mutSignatures from github
InstallmutSignatures <- function(){
  message("Installing mutSignatures from github...\n")
  # Install dependent package: corpor from CRAN
  utils::install.packages("corpcor")
  # install mutSignatures
  remotes::install_github("dami82/mutSignatures")
}



#' Run mutSignatures extraction and attribution on a spectra catalog file
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param CPU.cores Number of CPUs to use in running
#' \code{\link[mutSignatures]{decipherMutationalProcesses}}.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run mutSignatures. Setting seed can make the
#' attribution of mutSignatures repeatable.
#' Default: 1.
#'
#' @param K.exact \code{K.exact} is the exact value for
#' the number of signatures active in spectra (K).
#' Specify \code{K.exact} if you know exactly how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' @param nrun.exact number of NMF runs for extracting signatures and inferring
#' exposures.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The inferred exposure of \code{mutSignatures}, invisibly.
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

RunmutSignatures <-
  function(input.catalog,
           out.dir,
           CPU.cores = NULL,
           seedNumber = 12345,
           K.exact = NULL,
           nrun.exact = 1000,
           test.only = FALSE,
           overwrite = FALSE) {

    # Check whether K.exact is specified as a numeric element.
    stopifnot(is.numeric(K.exact) & length(K.exact) == 1)

    # Install mutSignatures, if failed to be loaded
    if (!requireNamespace("mutSignatures", quietly = TRUE)) {
      InstallmutSignatures()
    }

    # Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  # Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() # Save the random number generator (RNG) used

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


    # Read in spectra data from input.catalog file
    # spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                  strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]

    # Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    # convSpectra: convert the ICAMS-formatted spectra catalog
    # into a matrix which mutSignatures accepts:
    # 1. Remove the catalog related attributes in convSpectra
    # 2. Transpose the catalog
    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    dimnames(convSpectra) <- dimnames(spectra)
    # convSpectra must be converted to data.frame
    # before converting to mutation counts object.
    convSpectra <- as.data.frame(convSpectra)
    sample.number <- dim(spectra)[2]
    convSpectra <- mutSignatures::as.mutation.counts(convSpectra)



    ### Extract signatures when K is specified.
    params.obj <-
      mutSignatures::setMutClusterParams(
        # num signatures to extract
        num_processesToExtract = K.exact,
        # Number of matrix decompositions, each with bootstrapping once
        # before running.
        # Same as "nrun.exact" in other R-NMF-Brunet packages.
        num_totIterations = nrun.exact,
        # total num of cores to use (parallelization)
        num_parallelCores = CPU.cores,
        # Set number of random seed.
        seed = seedNumber)

    print(paste0("Assuming there are ",K.exact," signatures active in input spectra."))

    grDevices::pdf(paste0(out.dir,"/Silhouette.pdf"))
    extractionObject <- mutSignatures::decipherMutationalProcesses(
      input = convSpectra,
      params = params.obj)
    grDevices::dev.off()




    # Output extracted signatures in ICAMS format
    # Normalize the extracted signatures so that frequencies of each signature sums up to 1
    signatureObj <- extractionObject$Results$signatures
    extractedSignatures <- as.matrix(signatureObj@mutationFreq)
    rownames(extractedSignatures) <- signatureObj@mutTypes[,1]
    colnames(extractedSignatures) <- signatureObj@signatureId[,1]
    extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                             region = "unknown",
                                             catalog.type = "counts.signature")
    # Write extracted signatures
    ICAMS::WriteCatalog(extractedSignatures,
                        paste0(out.dir,"/extracted.signatures.csv"))


    # Derive exposure count attribution results.
    exposureObj <- extractionObject$Results$exposures
    # Normalized exposures
    exposureCounts <- as.matrix(exposureObj@exposures)
    rownames(exposureCounts) <- exposureObj@signatureId[,1]
    colnames(exposureCounts) <- exposureObj@sampleId[,1]

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







#' Run mutSignatures attribution on a spectra catalog file
#' and known signatures.
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param gt.sigs.file File containing input mutational signatures.
#' Columns are signatures, rows are mutation types.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run mutSignatures. Setting seed can make the
#' attribution of mutSignatures repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The inferred exposure of \code{mutSignatures}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export

RunmutSignaturesAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    # Install mutSignatures, if failed to be loaded
    if (!requireNamespace("mutSignatures", quietly = TRUE)) {
      InstallmutSignatures()
    }

    # Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  # Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() # Save the random number generator (RNG) used


    # Read in spectra data from input.catalog file
    # spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                  strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]


    # Read in ground-truth signature file
    gtSignatures <- ICAMS::ReadCatalog(gt.sigs.file)

    # Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    # Convert ICAMS-formatted spectra and signatures
    # into mutSignatures format
    # Requires removal of redundant attributes.
    convSpectra <- spectra
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    class(convSpectra) <- "matrix"
    # convSpectra must be converted to data.frame
    # before converting to mutation counts object.
    convSpectra <- as.data.frame(convSpectra)
    convSpectra <- mutSignatures::as.mutation.counts(convSpectra)

    gtSignaturesMS <- gtSignatures
    attr(gtSignaturesMS,"catalog.type") <- NULL
    attr(gtSignaturesMS,"region") <- NULL
    class(gtSignaturesMS) <- "matrix"
    # gtSignaturesMS must be converted to data.frame
    # before converting to mutation counts object.
    gtSignaturesMS <- as.data.frame(gtSignaturesMS)
    gtSignaturesMS <- mutSignatures::as.mutation.signatures(gtSignaturesMS)

    # Obtain inferred exposures using resolveMutSignatures function
    run <- mutSignatures::resolveMutSignatures(mutCountData = convSpectra,
                                               signFreqData = gtSignaturesMS)
    # An S4 object storing exposures, names of signatures and samples.
    exposures <- run$results$count.result

    # Write exposure counts in ICAMS and SynSig format.
    exposureCounts <- as.matrix(exposures@exposures)
    rownames(exposureCounts) <- exposures@signatureId[,1]
    colnames(exposureCounts) <- exposures@sampleId[,1]

    SynSigGen::WriteExposure(exposureCounts,
                             paste0(out.dir,"/inferred.exposures.csv"))

    # Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    # Save seeds and session information
    # for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) # Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) # Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) # Save seed in use to a text file

    # Return inferred exposures
    invisible(exposureCounts)
  }
