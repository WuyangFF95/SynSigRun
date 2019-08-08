#' Install SignatureEstimation package from URL source.
#'
#' @keywords internal
InstallSignatureEstimation <- function(){
  message("Installing SignatureEstimation from URL source...\n")

  devtools::install_url("https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz")
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
#' @param read.catalog.function Function to read a catalog
#' (can be spectra or signature catalog): it takes a file path as
#' its only argument and returning a catalog as a numeric matrix.
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
#' @return The attributed exposure of \code{SignatureEstimation}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export

RunSignatureEstimationQPAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           read.catalog.function,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Install SignatureEstimation from Bioconductor, if not found in library.
    if("SignatureEstimation" %in% rownames(installed.packages()) == FALSE)
      InstallSignatureEstimation()


    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- read.catalog.function(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]


    ## Read in ground-truth signatures
    ## gtSignatures: signature data.frame in ICAMS format
    gtSignatures <- read.catalog.function(gt.sigs.file)

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Convert ICAMS-formatted spectra and signatures
    ## into SignatureEstimation format
    ## Requires removal of redundant attributes.
    convSpectra <- spectra
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    attr(convSpectra,"class") <- NULL
    convSpectra <- as.matrix(convSpectra)
    dimnames(convSpectra) <- dimnames(spectra)

    gtSignaturesSE <- gtSignatures
    attr(gtSignaturesSE,"catalog.type") <- NULL
    attr(gtSignaturesSE,"region") <- NULL
    attr(gtSignaturesSE,"class") <- NULL
    gtSignaturesSE <- as.matrix(gtSignaturesSE)
    dimnames(gtSignaturesSE) <- dimnames(gtSignatures)

    ## Obtain attributed exposures using decomposeQP/decomposeSA function
    ## Note: SignatureEstimation::decomposeQP/decomposeSA() can only attribute ONE tumor at each run!
    num.tumors <- ncol(convSpectra)
    ## In each cycle, obtain attributed exposures for each tumor.
    exposureCounts <- data.frame()

    for(ii in 1:num.tumors){
      outputList <- SignatureEstimation::findSigExposures(
        M = convSpectra[, ii,drop = F],
        P = gtSignaturesSE,
        decomposition.method = decomposeQP)


      ## Obtain absolute exposure counts for current tumor
      exposuresOneTumor <- outputList$exposures ## Relative exposures (exposure proportions)
      exposuresOneTumor <- exposuresOneTumor * sum(convSpectra[,ii,drop = FALSE])

      ## Bind exposures for current tumor to exposure data.frame
      exposureCounts <- rbind(exposureCounts,t(exposuresOneTumor))
    }
    exposureCounts <- t(exposureCounts)

    ## Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    ## Write attributed exposures into a SynSig formatted exposure file.
    WriteExposure(exposureCounts,
                  paste0(out.dir,"/attributed.exposures.csv"))

    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return the exposures attributed, invisibly
    invisible(exposureCounts)
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
#' @param read.catalog.function Function to read a catalog
#' (can be spectra or signature catalog): it takes a file path as
#' its only argument and returning a catalog as a numeric matrix.
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
#' @return The attributed exposure of \code{SignatureEstimation}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export

RunSignatureEstimationSAAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           read.catalog.function,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Install SignatureEstimation from Bioconductor, if not found in library.
    if("SignatureEstimation" %in% rownames(installed.packages()) == FALSE)
      InstallSignatureEstimation()


    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- read.catalog.function(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]


    ## Read in ground-truth signatures
    ## gtSignatures: signature data.frame in ICAMS format
    gtSignatures <- read.catalog.function(gt.sigs.file)

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Convert ICAMS-formatted spectra and signatures
    ## into SignatureEstimation format
    ## Requires removal of redundant attributes.
    convSpectra <- spectra
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    attr(convSpectra,"class") <- NULL
    convSpectra <- as.matrix(convSpectra)
    dimnames(convSpectra) <- dimnames(spectra)

    gtSignaturesSE <- gtSignatures
    attr(gtSignaturesSE,"catalog.type") <- NULL
    attr(gtSignaturesSE,"region") <- NULL
    attr(gtSignaturesSE,"class") <- NULL
    gtSignaturesSE <- as.matrix(gtSignaturesSE)
    dimnames(gtSignaturesSE) <- dimnames(gtSignatures)

    ## Obtain attributed exposures using decomposeQP/decomposeSA function
    ## Note: SignatureEstimation::decomposeQP/decomposeSA() can only attribute ONE tumor at each run!
    num.tumors <- ncol(convSpectra)
    ## In each cycle, obtain attributed exposures for each tumor.
    exposureCounts <- data.frame()

    for(ii in 1:num.tumors){
      outputList <- SignatureEstimation::findSigExposures(
        M = convSpectra[, ii,drop = F],
        P = gtSignaturesSE,
        decomposition.method = decomposeSA)


      ## Obtain absolute exposure counts for current tumor
      exposuresOneTumor <- outputList$exposures ## Relative exposures (exposure proportions)
      exposuresOneTumor <- exposuresOneTumor * sum(convSpectra[,ii,drop = FALSE])

      ## Bind exposures for current tumor to exposure data.frame
      exposureCounts <- rbind(exposureCounts,t(exposuresOneTumor))
    }
    exposureCounts <- t(exposureCounts)


    ## Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    ## Write attributed exposures into a SynSig formatted exposure file.
    WriteExposure(exposureCounts,
                  paste0(out.dir,"/attributed.exposures.csv"))

    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return the exposures attributed, invisibly
    invisible(exposureCounts)
  }
