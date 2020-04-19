#' Run hdp extraction and attribution on a spectra catalog file
#'
#' @param input.catalog A catalog of spectra catalog
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#' @param CPU.cores Number of CPUs to use in running
#'    \code{\link[hdp]{hdp_posterior}}.
#'
#' @param seedNumber Specify the random seed for repeatable results.
#'
#' @param K.guess Suggested initial value of the number of
#' signatures, passed to \code{\link[hdp]{dp_activate}} as
#' \code{initcc}.
#'
#' @param multi.types A logical scalar or
#' a character vector.
#' If \code{FALSE}, hdp will regard all input spectra as one tumor type,
#' and will allocate them to one single dirichlet process node.
#'
#' If \code{TRUE}, hdp will infer tumor types based on the string before "::" in their names.
#' e.g. Tumor type for "SA.Syn.Ovary-AdenoCA::S.500" would be "SA.Syn.Ovary-AdenoCA"
#'
#' If it is a character vector, it should be a vector of case-sensitive tumor
#' types.
#' e.g. \code{c("SA.Syn.Ovary-AdenoCA", "SA.Syn.Ovary-AdenoCA", "SA.Syn.Kidney-RCC")}.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @param num.posterior Number of posterior sampling chains; can set to
#'   1 for testing.
#'
#' @param posterior.verbosity Pass to \code{\link[hdp]{hdp_posterior}}
#'      \code{verbosity}.
#'
#' @return The attributed exposure of \code{hdp}, invisibly.
#'
#' @export

RunhdpInternal <-
  function(input.catalog,
           # out.dir,
           CPU.cores           = 1,
           seedNumber          = 1,
           K.guess,
           multi.types         = FALSE,
           verbose             = TRUE,
           num.posterior       = 4,
           posterior.verbosity = 0) {

    if (!exists("stir.closure", envir = .GlobalEnv)) {
      assign("stir.closure", hdp::make.stirling(), envir = .GlobalEnv)
    }

    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  # To document the seed used
    RNGInUse <- RNGkind()      # To document the random number generator (RNG) used

    # input.catalog into a matrix that accepts.
    # hdp gets confused if the class in not matrix.
    convSpectra <- input.catalog
    class(convSpectra) <- "matrix"
    convSpectra <- t(convSpectra)

    number.channels <- nrow(input.catalog)
    number.samples  <- ncol(input.catalog)

    if (verbose) {
      message("Guessed number of signatures ",
              "(= Dirichlet process data clusters) = ", K.guess)
    }

    # Initialize hdp object
    # Allocate process index for hdp initialization.
    # Each different index number refers to a dirichlet process
    # for one tumor type. TODO Wuyang -- this comment does not seem correct
    if(multi.types == FALSE){ ## All tumors belong to one tumor type (default)
      num.tumor.types <- 1
      process.index <- c(0,1,rep(2,number.samples))
    } else if(multi.types == TRUE){
      ## There are multiple tumors in the sample.
      ## Tumor type will be inferred by the string before "::" in the column names.
      ## e.g. Tumor type for "SA.Syn.Ovary-AdenoCA::S.500" would be "SA.Syn.Ovary-AdenoCA"
      tumor.types <- sapply(
        colnames(input.catalog),
        function(x) {strsplit(x,split = "::",fixed = T)[[1]][1]})
      num.tumor.types <- length(unique(tumor.types))
      ## 0 refers to the grandparent DP node. All signatures are drawn from this node.
      ## Signature of each tumor type is drawn from a parent DP node (level 1).
      ## If a dataset has X tumor types, then we need to specify X level-1 nodes.
      process.index <- c(0, rep(1,num.tumor.types))
      ## For every tumor of the 1st/2nd/3rd/... tumor type,
      ## we need to specify a level 2/3/4/... DP node for the tumor.
      process.index <- c(process.index,1 + as.numeric(as.factor(tumor.types)))
      # Something like
      # c(0, 1, 1, 2, 2, 2, 3, 3)
      # 0 is grandparent
      # 1 is a parent of one type (there are 2 types)
      # 2 indcates tumors of the first type
      # 3 indicates tumors of second type
    } else if (is.character(multi.types)){ ## multi.types is a character vector recording tumor types
      num.tumor.types <- length(unique(multi.types))
      process.index <- c(0, rep(1,num.tumor.types))
      process.index <- c(process.index, 1 + as.numeric(as.factor(multi.types)))
    } else {
      stop("Error. multi.types should be TRUE, FALSE, or a vector of tumor types for each tumor sample.\n")
    }

    ## Specify ppindex as process.index,
    ## and cpindex (concentration parameter) as 1 + process.index
    ppindex <- process.index
    cpindex <- 1 + process.index

    ## Calculate the number of levels in the DP node tree.
    dp.levels <- length(unique(ppindex))

    al <- rep(1,dp.levels)

    ## initialise hdp
    if (verbose) message("calling hdp_init")
    hdpObject <- hdp::hdp_init(ppindex = ppindex,
                               cpindex = cpindex,
                               hh = rep(1,number.channels),
                               alphaa = al,
                               alphab = al)

    # num.process is the number of samples plus number of cancer types plus 1 (grandparent)
    num.process <- hdp::numdp(hdpObject)

    if (verbose) message("calling hdp_setdata")
    hdpObject <- hdp::hdp_setdata(
      hdpObject,
      (1 + num.tumor.types + 1):num.process,
      convSpectra)

    if (verbose) message("calling dp_activate")

    # dp_activate requires that stir.closure exists in .GlobalEnv
    hdpObject <- hdp::dp_activate(hdpObject,
                                  1:num.process,
                                  initcc = K.guess,
                                  seed = seedNumber)

    # Run num.posterior independent sampling chains
    f_posterior <- function(seed) {
      if (verbose) message("calling hdp_posterior")
      retval <- hdp::hdp_posterior (
        hdp       = hdpObject,
        verbosity = posterior.verbosity,
        # The remaining values, except seed, are from the vignette
        burnin    = 4000,
        n         = 50,
        space     = 50,
        cpiter    = 3,
        seed      = seed)
      return(retval)
    }

    chlist <- parallel::mclapply(
      # Must choose a different seed for each of the chains
      X = (seedNumber + 1:num.posterior * 10^6) %% 10^7 ,
      FUN = f_posterior,
      mc.cores = CPU.cores)

    # Generate the original multi_chain for the sample
    if (verbose) message("calling hdp_multi_chain")
    mut_example_multi <- hdp::hdp_multi_chain(chlist)

    if (verbose) message("calling hdp_extract_components")
    # Extract components(here is signature) with cosine.merge = 0.90 (default) TODO: Wu Yang, this is clustering of signatures?
    mut_example_multi_extracted <-
      hdp::hdp_extract_components(mut_example_multi)

    if (verbose) message("calling hdp::comp_categ_distn")
    extractedSignatures <-
      t(hdp::comp_categ_distn(mut_example_multi_extracted)$mean)

    rownames(extractedSignatures) <- rownames(input.catalog)
    # Set signature names to "hdp.0","hdp.1","hdp.2", ...
    colnames(extractedSignatures) <-
      paste("hdp", colnames(extractedSignatures), sep = ".")

    ## Calculate the exposure probability of each signature (component) for each
    ## tumor sample (posterior sample corresponding to a dirichlet process node).
    ## This is the probability distribution of signatures (components) for all
    ## tumor samples (DP nodes); exposureProbs is the normalized
    ## signature exposure all tumor samples # TODO Wuyang, what do you mean
    # by normalize?

    if (verbose) message("Calling hdp::comp_dp_distn to generate exposure probability")
    exposureProbs <- hdp::comp_dp_distn(mut_example_multi_extracted)$mean

    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.
    # Transpose so it conforms to SynSigEval format
    exposureProbs <- t(exposureProbs[-(1:(num.tumor.types + 1)), ])
    # Now in exposureProb rows are signatures, columns are samples

    # Calculate exposure counts from exposure probabilies and total mutation
    # counts
    exposureCounts <- exposureProbs %*% diag(rowSums(convSpectra))
    colnames(exposureCounts) <- colnames(input.catalog)
    rownames(exposureCounts) <- colnames(extractedSignatures)

    invisible(list(signature       = extractedSignatures,
                   exposure        = exposureCounts,
                   exposure.p      = exposureProbs,
                   multi.chains    = mut_example_multi,
                   ex.multi.chains = mut_example_multi_extracted,
                   seedInUse       = seedInUse,
                   RNGInUse        = RNGInUse))
  }
