#' Install MutationalPatterns from Bioconductor,
#' and its dependent package, NMF.
#'
#' @keywords internal
InstallMutationalPatterns <- function(){
  message("Installing MutationalPatterns from Bioconductor...")

  if (!requireNamespace("BiocManager", quietly = TRUE))
    utils::install.packages("BiocManager")
  BiocManager::install("MutationalPatterns")

  if ("NMF" %in% rownames(utils::installed.packages()) == FALSE)
    utils::install.packages("NMF")
}


#' Run MutationalPatterns attribution on a spectra catalog file
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
#' used to run MutationalPatterns. Setting seed can make the
#' attribution of MutationalPatterns repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The inferred exposure of \code{MutationalPatterns}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export

RunMutationalPatternsAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Install MutationalPatterns, if not found in library
    if ("MutationalPatterns" %in% rownames(utils::installed.packages()) == FALSE)
      InstallMutationalPatterns()


    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <-ICAMS::ReadCatalog(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]

    ## Remove the catalog related attributes in convSpectra
    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    dimnames(convSpectra) <- dimnames(spectra)


    ## Read in ground-truth signature file
    ## gt.sigs: signature data.frame in ICAMS format
    gtSignatures <-ICAMS::ReadCatalog(gt.sigs.file)
    ## Remove the catalog related attributes in gtSignatures
    tmp <- dimnames(gtSignatures)
    class(gtSignatures) <- "matrix"
    attr(gtSignatures,"catalog.type") <- NULL
    attr(gtSignatures,"region") <- NULL
    dimnames(gtSignatures) <- tmp

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }


    ## Derive exposure count attribution results.
    ## WARNING: MutationalPatterns can only do exposure attribution
    ## using SBS96 spectra catalog and signature catalog!
    exposureObject <-
      MutationalPatterns::fit_to_signatures(mut_matrix = convSpectra,
                                            signatures = gtSignatures)
    ## exposure attributions (in mutation counts)
    exposureCounts <- exposureObject$contribution
    ## Write exposure counts in ICAMS and SynSig format.
    SynSigGen::WriteExposure(exposureCounts,
                  paste0(out.dir,"/inferred exposures.csv"))

    ## Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return inferred exposures
    invisible(exposureCounts)
  }



#' Run MutationalPatterns extraction and attribution on a spectra catalog file
#'
#' WARNING: MutationalPatterns can only do exposure attribution
#' using SBS96 spectra catalog and signature catalog!
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param CPU.cores Number of CPUs to use in running
#' MutationalPatterns. For a server, 30 cores would be a good
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
#' of length 2 which tell MutationalPatterns to search the best
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
#' @return The inferred exposure of \code{MutationalPatterns}, invisibly.
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

RunMutationalPatterns <-
  function(input.catalog,
           out.dir,
           CPU.cores = NULL,
           K.exact = NULL,
           K.range = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Check whether ONLY ONE of K.exact or K.range is specified.
    bool1 <- is.numeric(K.exact) & is.null(K.range)
    bool2 <- is.null(K.exact) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    ## Install MutationalPatterns, if not found in library
    if ("MutationalPatterns" %in% rownames(utils::installed.packages()) == FALSE)
      InstallMutationalPatterns()


    ## Set seed
    seedNumber <- 123456
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                  strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]
    ## convSpectra: convert the ICAMS-formatted spectra catalog
    ## into a matrix which HDP accepts:
    ## 1. Remove the catalog related attributes in convSpectra
    ## 2. Transpose the catalog
    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## CPU.cores specifies number of CPU cores to use.
    ## If CPU.cores is not specified, CPU.cores will
    ## be equal to the minimum of 30 or (total cores)/2
    if(is.null(CPU.cores)){
      CPU.cores = min(30,(parallel::detectCores())/2)
    } else {
      stopifnot(is.numeric(CPU.cores))
    }


    ## Before running NMF packge,
    ## Load it explicitly to prevent errors.
    requireNamespace("NMF")

    ## Run NMF using ICAMS-formatted spectra catalog
    ## Determine the best number of signatures (K.best).
    ## If K.exact is provided, use K.exact as the K.best.
    ## If K.range is provided, determine K.best by doing raw extraction.
    if(bool1){
      K.best <- K.exact
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){
      K.range <- seq.int(K.range[1],K.range[2]) ## Change K.range to a full vector
      gof_nmf <- NMF::nmf(convSpectra,
                          rank = K.range,     ## Rank specifies number of signatures you want to assess
                          nrun = 200,
                          method = "brunet",  ## "brunet" is the default NMF method in NMF package.
                          .options = paste0("p", CPU.cores),
                          seed = seedNumber)
      gc()
      gc()
      gc()


      ## Choose the best signature number (K.best) active in the spectra
      ## catalog (input.catalog).
      ##
      ## According to paper "A flexible R package for nonnegative matrix factorization"
      ## (Gaujoux & Seoighe, 2010), the most common approach to choose number of
      ## signature (K, a.k.a. rank in this paper) is to choose the smallest K for which
      ## cophenetic correlation coefficient starts decreasing.
      for(current.K in K.range)
      {
        current.summary <- NMF::summary(gof_nmf$fit[[as.character(current.K)]])
        current.cophenetic.coefficient <- current.summary["cophenetic"]

        next.summary <- NMF::summary(gof_nmf$fit[[as.character(current.K+1)]])
        next.cophenetic.coefficient <- next.summary["cophenetic"]

        if(current.cophenetic.coefficient > next.cophenetic.coefficient)
          break
      }
      K.best <- current.K ## Choose K.best as the smallest current.K whose cophenetic
                          ## is greater than cophenetic from (current.K+1).
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))
    }


    ## Generates a list contain extracted signatures
    sigs_nmf <- MutationalPatterns::extract_signatures(
      mut_matrix = convSpectra,
      rank = K.best,
      nrun = 200)
    gc()
    gc()
    gc()
    ## names(sigs_nmf)
    ## [1] "signatures"    "contribution"  "reconstructed"
    sigsRaw <- sigs_nmf$signatures ## un-normalized signature matrix
    extractedSignatures <- t(t(sigsRaw) / colSums(sigsRaw))   ## normalize each signature's sum to 1
    ## Add signature names for signature matrix extractedSignatures
    colnames(extractedSignatures) <-
      paste("MutationalPatterns",1:ncol(extractedSignatures),sep=".")
    extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                             region = "unknown",
                                             catalog.type = "counts.signature")


    ## Output extracted signatures in ICAMS format
    ICAMS::WriteCatalog(extractedSignatures,
                           paste0(out.dir,"/extracted.signatures.csv"))


    ## Derive exposure count attribution results.
    ## WARNING: MutationalPatterns can only do exposure attribution
    ## using SBS96 spectra catalog and signature catalog!
    exposureObject <- MutationalPatterns::fit_to_signatures(mut_matrix = convSpectra,
                                                            signatures = extractedSignatures)
    gc()
    gc()
    gc()
    ## exposure attributions (in mutation counts)
    exposureCounts <- (exposureObject$contribution)
    ## Write exposure counts in ICAMS and SynSig format.
    SynSigGen::WriteExposure(exposureCounts,
                  paste0(out.dir,"/inferred exposures.csv"))


    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return a list of signatures and exposures
    invisible(list("signature" = extractedSignatures,
                   "exposure" = exposureCounts))
  }
