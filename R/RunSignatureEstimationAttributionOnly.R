#' Install SignatureEstimation package from URL source.
#'
#' @keywords internal
InstallSignatureEstimation <- function(){
  message("Installing SignatureEstimation from URL source...\n")
  remotes::install_url("https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz")
}


#' Run SignatureEstimation Quadratic Programming (QP) attribution
#' on a spectra catalog file and known signatures.
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
#' used to run SignatureEstimation. Setting seed can make the
#' attribution of SignatureEstimation repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return Invisibly returns a list which contains: \itemize{
#' \item $exposuresCounts: the exposure counts inferred in ICAMSxtra format,
#' \item $exposureErrors: the MSE in ICAMSxtra format,
#' \item $SEoutput: A list which contains: \itemize{
#'   \item $exposures: exposure proportion in SignatureEstimation format,
#'   and errors invisibly.
#'   \item $errors: mean squared error (MSE) between normalized reconstructed
#'    spectra and normalized ground-truth mutational spectra.}
#' }
#'
#' @importFrom utils capture.output
#'
#' @export

RunSignatureEstimationQPAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    # Install SignatureEstimation from Bioconductor, if not found in library.
    if("SignatureEstimation" %in% rownames(utils::installed.packages()) == FALSE)
      InstallSignatureEstimation()


    # Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  # Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() # Save the random number generator (RNG) used


    # Read in spectra data from input.catalog file
    # spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                  strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]


    # Read in ground-truth signatures
    # gtSignatures: signature data.frame in ICAMS format
    gtSignatures <- ICAMS::ReadCatalog(gt.sigs.file)

    # Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    # Convert ICAMS-formatted spectra and signatures
    # into SignatureEstimation format
    # Requires removal of redundant attributes.
    convSpectra <- spectra
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    class(convSpectra) <- "matrix"

    gtSignaturesSE <- gtSignatures
    attr(gtSignaturesSE,"catalog.type") <- NULL
    attr(gtSignaturesSE,"region") <- NULL
    class(gtSignaturesSE) <- "matrix"

    # Obtain inferred exposures using decomposeQP/decomposeSA function


    SEoutput <- SignatureEstimation::findSigExposures(
        M = convSpectra,
        P = gtSignaturesSE,
        decomposition.method = SignatureEstimation::decomposeQP)


    # Obtain absolute exposure counts for tumos
    # from SEoutput$exposures which refers to relative exposures (exposure proportions)
    exposureCounts <- t(SEoutput$exposures)
    exposureCounts <- exposureCounts * colSums(convSpectra)
    exposureCounts <- t(exposureCounts)

    # Obtain absolute mean squared error (MSE) for tumors
    # from SEoutput$errors which refers to relative MSE
    exposureErrors <- t(SEoutput$errors)
    exposureErrors <- exposureErrors * colSums(convSpectra)
    exposureErrors <- t(exposureErrors)


    # Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    # Write inferred exposures into a SynSig formatted exposure file.
    SynSigGen::WriteExposure(exposureCounts,
                  paste0(out.dir,"/inferred.exposures.csv"))

    # Write inferred exposures into a SynSig formatted exposure file.
    SynSigGen::WriteExposure(exposureErrors,
                             paste0(out.dir,"/exposure.MSE.errors.csv"))

    # Save seeds and session information
    # for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) # Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) # Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) # Save seed in use to a text file

    # exposuresCounts: the exposure counts inferred in ICAMSxtra format,
    # exposureErrors: the MSE in ICAMSxtra format.
    # SEoutput:
    # $exposures: exposure proportion in SignatureEstimation format,
    # and errors invisibly.
    # $errors: mean squared error (MSE) between normalized reconstructed
    # spectra and normalized ground-truth mutational spectra.
    retval <- list(
      exposureCounts = exposureCounts,
      exposureErrors = exposureErrors,
      SEoutput = SEoutput
    )

    invisible(retval)
  }


#' Run SignatureEstimation Simulated Annealing (SA) attribution
#' on a spectra catalog file and known signatures.
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
#' used to run SignatureEstimation. Setting seed can make the
#' attribution of SignatureEstimation repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return Invisibly returns a list which contains: \itemize{
#' \item $exposuresCounts: the exposure counts inferred in ICAMSxtra format,
#' \item $exposureErrors: the MSE in ICAMSxtra format,
#' \item $SEoutput: A list which contains: \itemize{
#'   \item $exposures: exposure proportion in SignatureEstimation format,
#'   and errors invisibly.
#'   \item $errors: mean squared error (MSE) between normalized reconstructed
#'    spectra and normalized ground-truth mutational spectra.}
#' }
#'
#' @importFrom utils capture.output
#'
#' @export

RunSignatureEstimationSAAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    # Install SignatureEstimation from Bioconductor, if not found in library.
    if("SignatureEstimation" %in% rownames(utils::installed.packages()) == FALSE)
      InstallSignatureEstimation()


    # Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  # Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() # Save the random number generator (RNG) used


    # Read in spectra data from input.catalog file
    # spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                  strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]


    # Read in ground-truth signatures
    # gtSignatures: signature data.frame in ICAMS format
    gtSignatures <- ICAMS::ReadCatalog(gt.sigs.file)

    # Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    # Convert ICAMS-formatted spectra and signatures
    # into SignatureEstimation format
    # Requires removal of redundant attributes.
    convSpectra <- spectra
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    class(convSpectra) <- "matrix"

    gtSignaturesSE <- gtSignatures
    attr(gtSignaturesSE,"catalog.type") <- NULL
    attr(gtSignaturesSE,"region") <- NULL
    class(gtSignaturesSE) <- "matrix"

    # Obtain inferred exposures using decomposeQP/decomposeSA function


    SEoutput <- SignatureEstimation::findSigExposures(
        M = convSpectra,
        P = gtSignaturesSE,
        decomposition.method = SignatureEstimation::decomposeSA)


    # Obtain absolute exposure counts for tumos
    # from SEoutput$exposures which refers to relative exposures (exposure proportions)
    exposureCounts <- t(SEoutput$exposures)
    exposureCounts <- exposureCounts * colSums(convSpectra)
    exposureCounts <- t(exposureCounts)

    # Obtain absolute mean squared error (MSE) for tumors
    # from SEoutput$errors which refers to relative MSE
    exposureErrors <- t(SEoutput$errors)
    exposureErrors <- exposureErrors * colSums(convSpectra)
    exposureErrors <- t(exposureErrors)


    # Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    # Write inferred exposures into a SynSig formatted exposure file.
    SynSigGen::WriteExposure(exposureCounts,
                             paste0(out.dir,"/inferred.exposures.csv"))

    # Write inferred exposures into a SynSig formatted exposure file.
    SynSigGen::WriteExposure(exposureErrors,
                             paste0(out.dir,"/exposure.MSE.errors.csv"))

    # Save seeds and session information
    # for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) # Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) # Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) # Save seed in use to a text file

    # exposuresCounts: the exposure counts inferred in ICAMSxtra format,
    # exposureErrors: the MSE in ICAMSxtra format.
    # SEoutput:
    # $exposures: exposure proportion in SignatureEstimation format,
    # and errors invisibly.
    # $errors: mean squared error (MSE) between normalized reconstructed
    # spectra and normalized ground-truth mutational spectra.
    retval <- list(
      exposureCounts = exposureCounts,
      exposureErrors = exposureErrors,
      SEoutput = SEoutput
    )

    invisible(retval)
  }
