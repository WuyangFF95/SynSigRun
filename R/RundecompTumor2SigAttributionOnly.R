#' Install decompTumor2Sig from Biocondcutor
InstalldecompTumor2Sig <- function(){
  message("Installing decompTumor2Sig from Biocondcutor...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("decompTumor2Sig")

}
#' Run decompTumor2Sig attribution on a spectra catalog file
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

RundecompTumor2SigAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           read.catalog.function,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Install deconstructSigs, if not found in library.
    if("decompTumor2Sig" %in% rownames(installed.packages()) == FALSE)
      InstalldecompTumor2Sig()


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
    ## gt.sigs: signature data.frame in ICAMS format
    gtSignatures <- read.catalog.function(gt.sigs.file)

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
    ## To analyze Alexandrov-like spectra catalogs,
    ## decompoTumor2Sig requires signatures to be a LIST of
    ## probability vectors (sum equals to 1)
    convSpectraList <- list()
    G <- ncol(convSpectra)
    for (Gcurrent in 1:G){
      currentSpectrumName <- colnames(convSpectra)[Gcurrent]
      convSpectraList[[currentSpectrumName]] <- convSpectra[,Gcurrent]
      convSpectraList[[currentSpectrumName]] <-
        convSpectraList[[currentSpectrumName]] / sum(convSpectraList[[currentSpectrumName]])
    }

    gtSignaturesDT <- gtSignatures
    attr(gtSignaturesDT,"catalog.type") <- NULL
    attr(gtSignaturesDT,"region") <- NULL
    class(gtSignaturesDT) <- "matrix"

    ## To analyze Alexandrov-like signatures,
    ## decompoTumor2Sig requires signatures to be a LIST of
    ## probability vectors (sum equals to 1, tolerance is 1e-5!)
    gtSignaturesDTList <- list()
    K <- ncol(gtSignaturesDT)
    for (Kcurrent in 1:K){
      currentSigName <- colnames(gtSignaturesDT)[Kcurrent]
      gtSignaturesDTList[[currentSigName]] <- gtSignaturesDT[,Kcurrent]
      gtSignaturesDTList[[currentSigName]] <-
        gtSignaturesDTList[[currentSigName]] / sum(gtSignaturesDTList[[currentSigName]])
    }

    exposureList <- decompTumor2Sig::decomposeTumorGenomes(genomes = convSpectraList,
                                                       signatures = gtSignaturesDTList)
    ## Convert exposureList to exposureProb.
    exposureProb <- matrix(nrow = K, ncol = G)
    rownames(exposureProb) <- colnames(gtSignaturesDT)
    colnames(exposureProb) <- colnames(convSpectra)
    for(Gcurrent in 1:G){
      currentSpectrumName <- names(exposureList)[Gcurrent]
      exposureProb[,Gcurrent] <- exposureList[[Gcurrent]]
    }
    ## Convert exposureProb to exposureCounts
    exposureCounts <- exposureProb
    for(Gcurrent in 1:G){
      exposureCounts[,Gcurrent] <- exposureProb[,Gcurrent] * sum(convSpectra[,Gcurrent])
    }

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
