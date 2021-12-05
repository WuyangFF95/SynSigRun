#' Install sigminer from Bioconductor,
#' also installing its dependent package, NMF.
#'
#' @keywords internal
Installsigminer <- function(){
  message("Installing sigminer from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    utils::install.packages("BiocManager")
  BiocManager::install("sigminer")
}



#' Run sigminer extraction and attribution on a spectra catalog file
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
#' sigminer. For a server, 30 cores would be a good
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
#' @param K.max \code{K.max} is the maximum number of signatures
#' users expect to active in \code{input.catalog}.
#' As this approach cannot specify \code{K.exact}, you can specify
#' \code{K.max} = \code{K.exact} If you know exactly how many signatures
#' are active in the \code{input.catalog}.
#' On the other hand, you may specify \code{max(K.range)} if you don't
#' know how many signatures are active in the \code{input.catalog}.
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
#' } of \code{sigminer}, invisibly.
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

Runsigminer <-
  function(input.catalog,
           out.dir,
           CPU.cores = NULL,
           seedNumber = 1,
           K.max = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    # Install sigminer, if not found in library
    if ("sigminer" %in% rownames(utils::installed.packages()) == FALSE)
      Installsigminer()


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
    # into a matrix which sigminer accepts:
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

    # Run sigminer using ICAMS-formatted spectra catalog
    print(paste0("Assuming there are at most ",K.max," signatures active in input spectra."))
    assess <- sigminer::sig_auto_extract(
      nmf_matrix = convSpectra,
      result_prefix = "BayesNMF",
      destdir = paste0(out.dir,"/sigminer.data/"),
      method = "L1W.L2H",
      strategy = "stable",
      K0 = K.max,
      nrun = 20,
      niter = 2e+05,
      tol = 1e-07,
      cores = CPU.cores,
      optimize = TRUE,
      skip = FALSE,
      recover = FALSE)
    gc()
    gc()
    gc()
    # normalized signature matrix
    extractedSignatures <- assess$Signature.norm
    colnames(extractedSignatures) <-
      paste("sigminer",1:ncol(extractedSignatures),sep=".")
    extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                             region = "unknown",
                                             catalog.type = "counts.signature")
    # Output extracted signatures in ICAMS format
    ICAMS::WriteCatalog(extractedSignatures,
                        paste0(out.dir,"/extracted.signatures.csv"))


    # Derive exposure count attribution results.
    # normalized exposure matrix
    exposureCounts <- assess$Exposure
    rownames(exposureCounts) <-
      paste("sigminer",1:nrow(exposureCounts),sep=".")
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
