#' Install mSigAct from GitHub
InstallmSigAct <- function(){
  message("Installing mSigAct from master...\n")
  devtools::install_github(repo = "steverozen/mSigAct",
                           ref = "master",
			   build_vignettes = TRUE)
}

#' Run mSigAct attribution on a spectra catalog file
#' and known signatures.
#'
#' @param input.catalog File containing input spectra catalog. Columns are
#' samples (tumors), rows are mutation types.
#'
#' @param gt.sigs.file File containing input mutational signatures. Columns are
#' signatures, rows are mutation types.
#'
#' @param read.catalog.function Function to read a catalog
#' (can be spectra or signature catalog): it takes a file path as
#' its only argument and returning a catalog as a numeric matrix.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param CPU.cores Number of CPUs to use in running
#' sigfit. For a server, 30 cores would be a good
#' choice; while for a PC, you may only choose 2-4 cores.
#' By default (CPU.cores = NULL), the CPU.cores would be equal
#' to \code{(parallel::detectCores())/2}, total number of CPUs
#' divided by 2.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run mSigAct. Setting seed can make the
#' attribution of mSigAct repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The attributed exposure of \code{mSigAct}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export

RunmSigActAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           read.catalog.function,
           out.dir,
           CPU.cores = NULL,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Install mSigAct, if not found in library.
    if("mSigAct" %in% rownames(installed.packages()) == FALSE)
      InstallmSigAct()


    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- read.catalog.function(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]


    ## Read in ground-truth signature file
    ## gtSignatures: signature data.frame in ICAMS format
    gtSignatures <- read.catalog.function(gt.sigs.file)

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
      CPU.cores = min(10,(parallel::detectCores())/2)
    } else {
      stopifnot(is.numeric(CPU.cores))
      if(CPU.cores > 10) CPU.cores = 10
    }

    ## mSigAct accepts ICAMS-formatted spectra and signature catalog.
    ## No need to convert. to convert catalog
    estimatedExposure <-
    SparseAssignActivity(spectra = spectra,
                         sigs = gtSignatures,
                         mc.cores = CPU.cores)
    exposureCounts <- estimatedExposure$exposure

    ## Write exposure counts in ICAMS and SynSig format.
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
