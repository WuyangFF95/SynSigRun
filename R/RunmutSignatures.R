#' Install mutSignatures from github
InstallmutSignatures <- function(){
  message("Installing mutSignatures from github...\n")
  # Install dependent package: corpor from CRAN
  utils::install.packages("corpcor")
  if(!"devtools" %in% rownames(utils::installed.packages()))
    utils::install.packages("devtools")
  # install mutSignatures
  devtools::install_github("dami82/mutSignatures")
}
#' Run deconstructSigs attribution on a spectra catalog file
#' and known signatures.
#'
#' @param input.catalog File containing input spectra catalog. Columns are
#' samples (tumors), rows are mutation types.
#'
#' @param gt.sigs.file File containing input mutational signatures. Columns are
#' signatures, rows are mutation types.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run deconstructSigs. Setting seed can make the
#' attribution of deconstructSigs repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The attributed exposure of \code{deconstructSigs}, invisibly.
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

    ## Install deconstructSigs, if not found in library.
    if("mutSignatures" %in% rownames(utils::installed.packages()) == FALSE)
      InstallmutSignatures()


    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                  strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]


    ## Read in ground-truth signature file
    ## gt.sigs: signature data.frame in ICAMS format
    gtSignatures <- ICAMS::ReadCatalog(gt.sigs.file)

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Convert ICAMS-formatted spectra and signatures
    ## into deconstructSigs format
    ## Requires removal of redundant attributes.
    convSpectra <- spectra
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    class(convSpectra) <- "matrix"

    gtSignaturesMS <- gtSignatures
    attr(gtSignaturesMS,"catalog.type") <- NULL
    attr(gtSignaturesMS,"region") <- NULL
    class(gtSignaturesMS) <- "matrix"

    ## Obtain attributed exposures using resolveMutSignatures function
    run <- mutSignatures::resolveMutSignatures(mutCountData = convSpectra,
                                               signFreqData = gtSignaturesMS)
    ## An S4 object storing exposures, names of signatures and samples.
    exposures <- run$results$count.result

    ## Write exposure counts in ICAMS and SynSig format.
    exposureCounts <- exposures@exposures
    rownames(exposureCounts) <- exposures@signatureId[,1]
    colnames(exposureCounts) <- exposures@sampleId[,1]

    WriteExposure(exposureCounts,
                  paste0(out.dir,"/attributed.exposures.csv"))

    ## Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return attributed exposures
    invisible(exposureCounts)
  }
