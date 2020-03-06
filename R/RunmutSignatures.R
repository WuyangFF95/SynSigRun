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
#' @return The attributed exposure of \code{mutSignatures}, invisibly.
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

    ## Install mutSignatures, if not found in library.
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
    gtSignatures <- ICAMS::ReadCatalog(gt.sigs.file)

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Convert ICAMS-formatted spectra and signatures
    ## into mutSignatures format
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



#' Run mutSignatures extraction and attribution on a spectra catalog file
#'
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param algorithm NMF implementation used to to extract signatures and
#' attribute exposures. Only "alexa", "brunet" or "lin" is valid.
#'
#' "alexa" or "brunet": Jean-Philippe Brunet's implementation.
#' This is the most widely used NMF implementation for signature extraction.
#' DOI: 10.1073/pnas.0308531101
#'
#' "lin": Chih-Jen Lin's implementation.
#' DOI:10.1109/TNN.2007.895831
#'
#' Default: "alexa".
#'
#' @param CPU.cores Number of CPUs to use in running
#' sigfit. For a server, 30 cores would be a good
#' choice; while for a PC, you may only choose 2-4 cores.
#' By default (CPU.cores = NULL), the CPU.cores would be equal
#' to \code{(parallel::detectCores())/2}, total number of CPUs
#' divided by 2.
#'
#' @param iterations Number of iterations in signature extraction.
#' Default: 1000.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run sigfit. Setting seed can make the
#' attribution of sigfit repeatable.
#' Default: 1.
#'
#' @param K,K.range \code{K} is the precise value for
#' the number of signatures active in spectra (K).
#' Specify \code{K} if you know precisely how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell sigfit to search the best
#' signature number active in spectra, K, in this range of Ks.
#' Specify \code{K.range} if you don't know how many signatures
#' are active in the \code{input.catalog}.
#' K.max - K.min >= 3, otherwise an error will be thrown.
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
#' @return The attributed exposure of \code{mutSignatures}, invisibly.
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
           algorithm = "alexa",
           CPU.cores = NULL,
           iterations = 1000,
           seedNumber = 1,
           K = NULL,
           K.range = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Check whether ONLY ONE of K or K.range is specified.
    bool1 <- is.numeric(K) & is.null(K.range)
    bool2 <- is.null(K) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    ## Check if algorithm parameter is correctly set
    stopifnot(algorithm %in% c("alexa","brunet","lin"))


    ## Install mutSignatures, if not found in library
    if ("mutSignatures" %in% rownames(utils::installed.packages()) == FALSE)
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

    ## convSpectra: convert the ICAMS-formatted spectra catalog
    ## into a matrix which mutSignatures accepts:
    ## 1. Remove the catalog related attributes in convSpectra
    ## 2. Convert the spectra matrix into an data.frame
    ## 3. Based on converted spectra, create an S4 object "spectraCounts"
    ## of class "mutationCounts"

    convSpectra <- spectra
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    class(convSpectra) <- "matrix"
    dimnames(convSpectra) <- dimnames(spectra)
    convSpectra <- as.data.frame(convSpectra)
    sample.number <- dim(spectra)[2]
    spectraCounts <- mutSignatures::as.mutation.counts(x = convSpectra) ## S4 object. Its class is "mutationCounts"

    ## Determine the best number of signatures (K.best).
    ## If K is provided, use K as the K.best.
    ## If K.range is provided, determine K.best by doing raw extraction
    ## with function ().
    ## Function can find the K.best with the Sihouette.
    if(bool1){
      K.best <- K
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){
      sigErrors <- mutSignatures::prelimProcessAssess(
        input = spectraCounts,
        plot = F,
        maxProcess = K.range[2],
        approach = "counts")

      ## Select the K empirically by SynSigRun:
      ## Choose the minimum K which let reconstrected spectra
      ## to have <0.1% difference to the ground-truth spectra.
      K.best <- min(sigErrors[which(sigErrors[,2] < 0.001),1])
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))
    }


    ## Precise extraction:


    ## 1. In Params Object, specifying number of signatures,
    ## and iterating more times to get more precise extraction.
    ## For faster computation, change num_parallelCores to be
    ## CPU.cores.
    extrParams <- mutSignatures::setMutClusterParams(
      num_processesToExtract = K.best,
      approach = "counts",
      num_totIterations = iterations,
      num_parallelCores = CPU.cores, # set to 1 to avoid parallelization
      debug = FALSE,
      algorithm = algorithm) ## Use Brunet NMF
    ## 2. Precise extraction,
    ## and report signatures and exposures from preciseExtr object.

    ## precise extraction will draw a Silhouette coefficient plot.
    ## Need to turn on the PDF device.
    grDevices::pdf(paste0(out.dir,"/Silhouette.coeffs.pdf"))
    preciseExtr <- mutSignatures::decipherMutationalProcesses(input = spectraCounts, params = extrParams)
    grDevices::dev.off()
    extractedSignatures <- preciseExtr$Results$signatures@mutationFreq
    rownames(extractedSignatures) <- preciseExtr$Results$signatures@mutTypes[,1]
    colnames(extractedSignatures) <- preciseExtr$Results$signatures@signatureId[,1]

    exposureCounts <- preciseExtr$Results$exposures@exposures
    rownames(exposureCounts) <- preciseExtr$Results$exposures@signatureId[,1]
    colnames(exposureCounts) <- preciseExtr$Results$exposures@sampleId[,1]

    extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                             region = "unknown",
                                             catalog.type = "counts.signature")

    ## Write extracted signatures into a ICAMS signature catalog file.
    ICAMS::WriteCatalog(extractedSignatures,
                        paste0(out.dir,"/extracted.signatures.csv"))

    ## Write attributed exposures into a SynSig formatted exposure file.
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
