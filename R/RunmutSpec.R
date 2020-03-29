

#' Run mutSpec extraction and attribution on a spectra catalog file
#'
#' NOTE: mutSpec can only do exposure attribution
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
#' mutSpec. For a server, 30 cores would be a good
#' choice; while for a PC, you may only choose 2-4 cores.
#' By default (CPU.cores = NULL), the CPU.cores would be equal
#' to \code{(parallel::detectCores())/2}, total number of CPUs
#' divided by 2.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run mutSpec. Setting seed can make the
#' attribution of mutSpec repeatable.
#' Default: 1.
#'
#' @param K,K.range \code{K} is the precise value for
#' the number of signatures active in spectra (K).
#' Specify \code{K} if you know precisely how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell mutSpec to search the best
#' signature number active in spectra, K, in this range of Ks.
#' Specify \code{K.range} if you don't know how many signatures
#' are active in the \code{input.catalog}.
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
#' @return The attributed exposure of \code{mutSpec}, invisibly.
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

RunmutSpec <-
  function(input.catalog,
           out.dir,
           CPU.cores = NULL,
           seedNumber = 1,
           K = NULL,
           K.range = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Check whether ONLY ONE of K or K.range is specified.
    bool1 <- is.numeric(K) & is.null(K.range)
    bool2 <- is.null(K) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    ## Install NMF, if not found in library
    if ("NMF" %in% rownames(utils::installed.packages()) == FALSE)
      install.packages("NMF")


    ## Set seed
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
    ## CPU.cores will be capped at 30.
    ## If CPU.cores is not specified, CPU.cores will
    ## be equal to the minimum of 30 or (total cores)/2
    if(is.null(CPU.cores)){
      CPU.cores = min(30,(parallel::detectCores())/2)
    } else {
      stopifnot(is.numeric(CPU.cores))
      if(CPU.cores > 30) CPU.cores = 30
    }
    ## "P" means that if the program cannot run parallelly, the NMF will abort.
    ## Therefore, we use "p" instead.
    nbCPU   <- paste0("vp", CPU.cores)

    ## Before running NMF packge,
    ## Load it explicitly to prevent errors.
    requireNamespace("NMF")

    ## Run NMF using ICAMS-formatted spectra catalog
    ## Determine the best number of signatures (K.best).
    ## If K is provided, use K as the K.best.
    ## If K.range is provided, determine K.best by doing raw extraction.
    if(bool1){
      K.best <- K
      K.range <- K
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){

      # Estimate the number of signatures with our data

      nbSign  <- seq.int(K.range[1],K.range[2]) ## Change K.range to a full vector
      # The minum number of signatures can't be lower than 2

      estim_r <- NMF::nmf(matrixNMF, method="brunet", nbSign, nrun=50, .opt=nbCPU)

      # Shuffle original data
      v_random <- NMF::randomize(matrixNMF)
      # Estimate quality measures from the shuffled data
      estim_r_random <- nmf(v_random, method="brunet", nbSign, nrun=50, .opt=nbCPU)

      ## Garbage collection
      gc()
      gc()
      gc()


      # Plot the estimation for our data and the random ones
      grDevices::graphics.off()
      options(bitmapType='cairo')
      grDevices::png(opt$output, width=3000, height=2000, res=300)
      plot(estim_r, estim_r_random)
      invisible( grDevices::dev.off() )



      ## Choose the best signature number (K.best) active in the spectra
      ## catalog (input.catalog).
      ##
      ## According to paper "A flexible R package for nonnegative matrix factorization"
      ## (Gaujoux & Seoighe, 2010), the most common approach to choose number of
      ## signature (K, a.k.a. rank in this paper) is to choose the smallest K for which
      ## cophenetic correlation coefficient starts decreasing.
      for(current.K in K.range)
      {
        current.summary <- NMF::summary(estim_r$fit[[as.character(current.K)]])
        current.cophenetic.coefficient <- current.summary["cophenetic"]

        next.summary <- NMF::summary(estim_r$fit[[as.character(current.K+1)]])
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
    res  <- NMF::nmf(matrixNMF, K.best, "brunet", nrun=200, .opt=nbCPU)

    # Recover the matrix W and H
    matrixW <- NMF::basis(res)
    matrixH <- NMF::coef(res) ## un-normalized signature matrix
    extractedSignatures <- t(t(matrixW) / colSums(matrixW))   ## normalize each signature's sum to 1
    ## Add signature names for signature matrix extractedSignatures
    colnames(extractedSignatures) <-
      paste("mutSpec",1:ncol(extractedSignatures),sep=".")
    extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                             region = "unknown",
                                             catalog.type = "counts.signature")


    ## Output extracted signatures in ICAMS format
    ICAMS::WriteCatalog(extractedSignatures,
                        paste0(out.dir,"/extracted.signatures.csv"))


    ## Derive exposure count attribution results.
    ## WARNING: mutSpec can only do exposure attribution
    ## using SBS96 spectra catalog and signature catalog!

    ## exposure attributions (in mutation counts)
    exposureCounts <- matrixH
    for(ii in 1:ncol(exposureCounts)){
      exposureCounts[,ii] <- exposureCounts[,ii] * colSums(convSpectra)[ii]
    }
    ## Add signature names for signature matrix extractedSignatures
    rownames(exposureCounts) <-
      paste("mutSpec",1:nrow(exposureCounts),sep=".")

    ## Write exposure counts in ICAMS and SynSig format.
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
