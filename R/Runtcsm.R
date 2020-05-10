#' tcsm's internal function to calculate likelihood for each \code{K}
#' using heldout method.
#'
#' This code requires the CRAN package stm https://cran.r-project.org/package=stm.
#'
#' tcsm advises to use "heldout" method to obtain the likelihood of number of
#' mutational signatures, \code{K}.
#' That is, to split the whole spectra catalogs dataset into 2 parts:
#'
#' 80% of spectra catalogs should be allocated into \code{train.mc.data}
#' for training \code{\link[stm]{stm}} model,
#' 20% of spectra catalogs should be allocated into \code{test.mc.data}
#' for estimating likelihood for given K.
#'
#' @param train.mc.data Transposed \code{\link[ICAMS]{ICAMS}} spectra catalog
#' containing mutations of training samples.
#'
#' @param test.mc.data Transposed \code{\link[ICAMS]{ICAMS}} spectra catalog
#' containing mutations of testing samples.
#'
#' @return heldout, which is a named list
#'   $missing
#'     $index
#'     $docs - list where each element is a list containing the mutations for a given sample
#'   $documents - list where each element is a list containing the mutations for a given sample
#'   $vocab - list the vocabularies, which refers to mutation types.
#'
#' @keywords internal
make.heldout.obj <- function(train.mc.data, test.mc.data){
  mc.data <- rbind(train.mc.data, test.mc.data)
  mat <- data.matrix(mc.data)
  corpus <- stm::readCorpus(mat, type="dtm")
  prep <- stm::prepDocuments(corpus$documents, corpus$vocab, lower.thresh=0)
  n_training_samples = dim(train.mc.data)[1]
  n_test_samples = dim(test.mc.data)[1]
  heldout <- list()
  heldout$documents <- prep$documents[1:n_training_samples]
  heldout$missing <- list()
  heldout$missing$docs <- prep$documents[(n_training_samples+1):(n_training_samples+n_test_samples)]
  heldout$missing$index <- seq(n_training_samples+1, n_training_samples+n_test_samples)
  for (i in heldout$missing$index){
    heldout$documents[[i]] <- matrix(, nrow=2, ncol=0)
  }
  heldout$vocab <- prep$vocab
  heldout
}

run.stm <- function(
  train.mc.data,
  test.mc.data,
  train.feature.file = NULL,
  test.feature.file = NULL,
  covariates = NULL,
  K,
  seed){
  heldout <- make.heldout.obj(train.mc.data, test.mc.data)
  if (is.null(covariates)){
    stm1 <- stm::stm(documents=heldout$documents, vocab=heldout$vocab, K=K, seed=seed,
                     max.em.its = 500, init.type = "Spectral")
  } else {
    train.feature.data <- read.delim(train.feature.file, sep = '\t', header = TRUE, row.names=1)
    test.feature.data <- read.delim(test.feature.file, sep = '\t', header = TRUE, row.names=1)
    feature.data <- rbind(train.feature.data, test.feature.data)
    covariate.formula <- stats::as.formula(paste0("~", covariates))
    # the heldout object
    stm1 <- stm::stm(documents=heldout$documents, vocab=heldout$vocab, K=K, seed=seed,
                     prevalence = covariate.formula, max.em.its = 500, data=feature.data,
                     init.type = "Spectral")
  }

  heldout.performance <- stm::eval.heldout(stm1, heldout$missing)
  heldout.likelihood <- heldout.performance$expected.heldout
  df <- data.frame("likelihood" = heldout.likelihood, K)
  return(df)
}



#' Run tcsm extraction and attribution on a spectra catalog file
#'
#' WARNING: tcsm can only do exposure attribution
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
#' tcsm. For a server, 30 cores would be a good
#' choice; while for a PC, you may only choose 2-4 cores.
#' By default (CPU.cores = NULL), the CPU.cores would be equal
#' to \code{(parallel::detectCores())/2}, total number of CPUs
#' divided by 2.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run tcsm. Setting seed can make the
#' attribution of tcsm repeatable.
#'
#' @param K.exact,K.range \code{K.exact} is the exact
#' value for the number of signatures active in spectra (K).
#' Specify \code{K.exact} if you know exact how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell tcsm to search the best
#' signature number active in spectra, K, in this range of Ks.
#' Specify \code{K.range} if you don't know how many signatures
#' are active in the \code{input.catalog}.
#' K.max - K.min >= 3, otherwise an error will be thrown.
#'
#' WARNING: You must specify only one of \code{K} or \code{K.range}!
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#'
#' @param overwrite If TRUE, overwrite existing output.
#'
#' @return The inferred exposure of \code{tcsm}, invisibly.
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

Runtcsm <-
  function(input.catalog,
           out.dir,
           CPU.cores = NULL,
           seedNumber = 1,
           K.exact = NULL,
           K.range = NULL,
		       covariates = NULL,
           test.only = FALSE,
           overwrite = FALSE,
		   feature.file = NULL,
		   effect.output.file = NULL,
		   sigma.output.file = NULL,
		   gamma.output.file = NULL
		   ) {

    ## Check whether ONLY ONE of K or K.range is specified.
    bool1 <- is.numeric(K.exact) & is.null(K.range)
    bool2 <- is.null(K.exact) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)



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
      dir.create(paste0(out.dir,"/other.outputs.by.tcsm"), recursive = T)
    }

    ## CPU.cores specifies number of CPU cores to use.
    ## If CPU.cores is not specified, CPU.cores will
    ## be equal to the minimum of 30 or (total cores)/2
    if(is.null(CPU.cores)){
      CPU.cores = min(30,(parallel::detectCores())/2)
    } else {
      stopifnot(is.numeric(CPU.cores))
    }

    ## For faster computation, enable parallel computing.
    ##

    ## Use CPU.cores cores for parallel computing
    options(mc.cores = CPU.cores)

    ## convSpectra: convert the ICAMS-formatted spectra catalog
    ## into a matrix which tcsm accepts:
    ## 1. Remove the catalog related attributes in convSpectra
    ## 2. Transpose the catalog

    convSpectra <- spectra
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    dimnames(convSpectra) <- dimnames(spectra)
    sample.number <- dim(spectra)[2]
    convSpectra <- t(convSpectra)

    ## Determine the best number of signatures (K.best).
    ## If K is provided, use K as the K.best.
    ## If K.range is provided, determine K.best by running .
    if(bool1){
      K.best <- K.exact
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){

      ## Choose the best signature number (K.best) active in the spectra
      ## catalog (input.catalog).
      ## Raw extraction: estimate most likely number of signatures
      ## (Nsig.max + 1) number of elements in mcmc_samples_extr
      ## The first Nsig.min number of elements are NULL elements
      ## Nsig.min+1 to Nsig.max elements are list elements of two elements: $data and $result
      ## The last element is the best signature number
      K.range <- seq.int(K.range[1],K.range[2])

      likelihoods <- numeric(0)
      for(K in K.range){

        currLikelihood <- data.frame()
        portion <- ceiling(sample.number/5)

        ## Do heldout training 5 times.
        ## In each run, take ~80% of samples to train,
        ## and take ~20% of samples to calculate likelihood of current K.
        for(ii in 1:5){

          if(ii == 5){
            test.samples <-
              seq((ii - 1)*portion + 1,sample.number)
          } else{
            test.samples <-
              seq((ii - 1)*portion + 1,ii*portion)
          }
          train.samples <- setdiff(seq.int(1,sample.number),test.samples)

          train.mc.data <- convSpectra[train.samples,]
          test.mc.data <- convSpectra[test.samples,]

          currLikelihood <- rbind(currLikelihood,
            run.stm(
            train.mc.data,
            test.mc.data,
            train.feature.file = NULL,
            test.feature.file = NULL,
            covariates = NULL,
            K,
            seed = seedNumber))
          gc()
        }

        likelihoods[as.character(K)] <- mean(currLikelihood[,1])
	  }

      if(TRUE){ ## debug
        ## Choose the K.best if likelihood(K.best + 1) - likelihood(K.best) < 0.01
        for(K in K.range){
          K.best <- K
          if(K == max(K.range))
            break
          if(likelihoods[as.character(K+1)] - likelihoods[as.character(K)] < 0.01)
            break
        }
      } else{
        K.best <- names(likelihoods)[which.max(likelihoods)] ## Choose K.best
        K.best <- as.integer(K.best)
      }
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))
    }


    ## Precise extraction:
    ## Specifying number of signatures, and iterating more times to get more precise extraction
    ## Return a list with two elements: $data and $result
    corpus <- stm::readCorpus(convSpectra, type="dtm")
    prep <- stm::prepDocuments(corpus$documents, corpus$vocab)

    if (is.null(covariates)){
      stm1 <- stm::stm(documents=prep$documents, vocab=prep$vocab, K=K.best, seed=seedNumber,
                  max.em.its = 500, init.type = "Spectral", sigma=0)
    } else {
      covariate.formula <- as.formula(paste0("~", covariates))
      feature.data <- read.delim(feature.file, sep = '\t', header = TRUE, row.names=1)
      stm1 <- stm::stm(documents=prep$documents, vocab=prep$vocab, K=K.best, seed=seedNumber,
                  prevalence = covariate.formula, max.em.its = 500, data=feature.data,
                  init.type = "Spectral", sigma=0)
      effect <- estimateEffect(covariate.formula, stm1, metadata=feature.data)
      effect.summary <- summary(effect)
      effect.tables <- effect.summary$tables
      results <- lapply(effect.summary$tables, function(x) x[, "Estimate"])
      effect.frame <- as.data.frame(do.call(rbind, results))
      utils::write.table(effect.frame, file=paste0(out.dir,"/other.outputs.by.tcsm/effect.frame.tsv"), sep="\t")
      covariate.list <- strsplit(covariates, "\\+")[[1]]
      gamma <- stm1$mu$gamma
      rownames(gamma) <- c("default", covariate.list)
      utils::write.table(gamma, file=paste0(out.dir,"/other.outputs.by.tcsm/gamma.tsv"), sep="\t")
  }

  utils::write.table(stm1$sigma, file=paste0(out.dir,"/other.outputs.by.tcsm/sigma.tsv"), sep="\t")
  # process the signatures and save them
  # the K-by-V matrix logbeta contains the natural log of the probability of seeing each word conditional on the topic
  mat <- stm1$beta$logbeta[[1]]
  signatures <- apply(mat, 1:2, exp)
  ## Mutation types
  colnames(signatures) <- stm1$vocab
  ## Get Signature names from exposure object
  dt <- stm::make.dt(stm1)
  rownames(signatures) <- colnames(dt)[-1]
  ## signatures should be transposed to ICAMS format.
  extractedSignatures <- t(signatures)

  ## Change signature names for signature matrix extractedSignatures:
  ## E.g., replace "Signature A" with "tcsm.A".
  colnames(extractedSignatures) <-
    gsub(pattern = "Topic",replacement = "tcsm.",colnames(extractedSignatures))
  extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                           region = "unknown",
                                           catalog.type = "counts.signature")

  ## Write extracted signatures into a ICAMS signature catalog file.
  ICAMS::WriteCatalog(extractedSignatures,
                      paste0(out.dir,"/extracted.signatures.csv"))


  # Get raw exposure object, dt.
  # dt records exposure proportion in each tumor spectrum.
  # That is, the signature exposure of each tumor sums to 1.
  dt <- stm::make.dt(stm1)
  rawExposures <- as.matrix(dt[,-1])
  rownames(rawExposures) <- colnames(spectra)
  colnames(rawExposures) <-
    gsub(pattern = "Topic",replacement = "tcsm.",
	colnames(rawExposures))
  ## Calculate exposureCounts from rawExposures.
  exposureCounts <- rawExposures * colSums(spectra)
  exposureCounts <- t(exposureCounts)

  ## Write inferred exposures into a SynSig formatted exposure file.

  SynSigGen::WriteExposure(
    exposureCounts,
    paste0(out.dir,"/inferred.exposures.csv"))


    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return a list of signatures and exposures
    invisible(list("signature" = extractedSignatures,
                   "exposure" = exposureCounts))
  }
