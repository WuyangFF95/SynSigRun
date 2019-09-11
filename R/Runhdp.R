#' Install hdp from GitHub.
#'
#' @keywords internal
Installhdp <- function(){
  message("Installing hdp from GitHub nicolaroberts/hdp ...\n")
  devtools::install_github("nicolaroberts/hdp", build_vignettes = TRUE)
}



#' Run hdp extraction and attribution on a spectra catalog file
#'
#' WARNING: hdp can only do exposure attribution
#' using SBS96 spectra catalog and signature catalog!
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param read.catalog.function Function to read a catalog
#' (can be spectra or signature catalog): it takes a file path as
#' its only argument and returning a catalog as a numeric matrix.
#'
#' @param write.catalog.function Function to write a catalog.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run hdp. Setting seed can make the
#' attribution of hdp repeatable.
#' Default: 1.
#'
#' @param K,K.range \code{K} is the precise value for
#' the number of signatures active in spectra (K).
#' Specify \code{K} if you know precisely how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell hdp to search the best
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
#' @return The attributed exposure of \code{hdp}, invisibly.
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

Runhdp <-
  function(input.catalog,
           read.catalog.function,
           write.catalog.function,
           out.dir,
           seedNumber = 1,
           K = NULL,
           K.range = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Check whether ONLY ONE of K or K.range is specified.
    bool1 <- is.numeric(K) & is.null(K.range)
    bool2 <- is.null(K) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    ## Install hdp, if not found in library
    if ("hdp" %in% rownames(utils::installed.packages()) == FALSE)
      Installhdp()


    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- read.catalog.function(input.catalog,
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
    dimnames(convSpectra) <- dimnames(spectra)
    convSpectra <- t(convSpectra)

    number.channels <- dim(spectra)[1]
    number.samples <- dim(spectra)[2]

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Determine the best number of signatures (K.best).
    ##
    ## hdp accepts an initial guess of number of signatures (K.initial), and later
    ## determine the best number of signatures (K.best)
    ##
    ## If K is provided, use K as the K.initial.
    ## If K.range is provided, use the largest value as the K.initial.
    ##
    ##
    if(bool1)
      K.initial <- K
    if(bool2)
      K.initial <- max(K.range)

    print(paste0("Assuming there are ",K.initial," signatures active in input spectra."))
    print(paste0("But the final number of signatures may not equal to ",K.initial,"\n."))

    ## Run hdp main program.
    ## Step 1: initialize hdp object
    {
      ## initialise hdp
      ppindex <- c(0, 1, rep(2,number.samples))
      cpindex <- c(1, 2, rep(3,number.samples))

      hdpObject <- hdp::hdp_init(ppindex = ppindex,
                      cpindex = cpindex,
                      hh = rep(1,number.channels),
                      alphaa = rep(1,3), alphab = rep(1,3))
      num.process <- hdp::numdp(hdpObject)

      ## Add data
      hdpObject <- hdp::hdp_setdata(hdpObject, 3:num.process, convSpectra)

      hdp::numdp(hdpObject)

      ## hdp also has to enter number of signatures in advance, but the final result doesn't necessarily equal to the initial input value
      ## When the number of sample is too large, hdp may eat up all your RAMs!
      hdpObject <- hdp::dp_activate(hdpObject,
                                    1:num.process,
                                    K.initial,
                                    seed=seedNumber)	## K.initial initial components(start with 30 signatures)

      hdpObject

      ## Release the occupied RAM by dp_activate
      gc()
      gc()
      gc()
    }

    ## Step 2: run 4 independent sampling chains
    {
      ## Run four independent posterior sampling chains
      chlist <- vector("list", 4)	#4 is too much here!

      i <- 1 ## Should execute manually rather than for cycle, because it is too slow!
      for (i in 1:4){
        chlist[[i]] <- hdp::hdp_posterior(hdpObject,
                                     burnin=4000,	## 4000 is too large, but necessary
                                     n=50,	## n here refers to the times of posterior sampling after burnin. To be faster, n can be set to 50.
                                     space=50,	## space is the time of iterations between two samplings. In this case, I need to iterate 9000 times.
                                     cpiter=3,
                                     seed= (seedNumber + i * 10^6) %% (10^7) ) ## Cannot choose the same seed for 4 chains!
      }

      ## Generate the original multi_chain for the sample
      mut_example_multi <- hdp::hdp_multi_chain(chlist)
      mut_example_multi
    }

    ## Step 3: Plot the diagnostics of sampling chains.
    {
      ## Plotting using hdp functions
      ## Plotting device on the server does not work
      ## Need to plot the file into a pdf



      ## Draw the DP oscillation plot for mut_example_multi(original_sample)
      {
        grDevices::pdf(file = paste0(out.dir,"/original_sample.pdf"))

        graphics::par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
        p1 <- lapply(hdp::chains(mut_example_multi),
                     hdp::plot_lik, bty="L", start=500)
        p2 <- lapply(hdp::chains(mut_example_multi),
                     hdp::plot_numcluster, bty="L")

        grDevices::dev.off()
      }

      ## Extract components(here is signature) with cosine.merge = 0.90 (default)
      mut_example_multi_extracted <- hdp::hdp_extract_components(mut_example_multi)
      mut_example_multi_extracted


      ## Generate a pdf for mut_example_multi_extracted
      {
        grDevices::pdf(file = paste0(out.dir,"/signature_hdp_embedded_func.pdf"))
        ## Draw the DP oscillation plot for mut_example_multi_extracted
        graphics::par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
        p1 <- lapply(hdp::chains(mut_example_multi_extracted),
                     hdp::plot_lik, bty="L", start=500)
        p2 <- lapply(hdp::chains(mut_example_multi_extracted),
                     hdp::plot_numcluster, bty="L")

        ## Draw the computation size plot
        graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
        hdp::plot_comp_size(mut_example_multi_extracted, bty="L")

        ## Close the PDF device so that the plots are exported to PDF
        grDevices::dev.off()
      }
    }

    ## Step 4: Using hdp samples to extract signatures
      {

        ## Calculate the mutation composition in each signature:
        extractedSignatures <- hdp::comp_categ_distn(mut_example_multi_extracted)$mean
        dim(extractedSignatures)
        ## Add base context for extractedSignatures
        extractedSignatures <- t(extractedSignatures)
        rownames(extractedSignatures) <- rownames(spectra)
        ## Change signature names in extractedSignatures
        ## from "0","1","2" to "hdp.0","hdp.1","hdp.2"
        colnames(extractedSignatures) <-
          paste("hdp", colnames(extractedSignatures), sep = ".")
        ## Remove "hdp.0" (noise signature) if it is a null signature or NA signature
        flagRemoveHDP0 <- FALSE
        if(is.null(extractedSignatures[1,"hdp.0"]) | is.na(extractedSignatures[1,"hdp.0"]))
          flagRemoveHDP0 <- TRUE
        if(flagRemoveHDP0){
            sigToBeRemoved <- which(colnames(extractedSignatures) == "hdp.0")
            extractedSignatures <- extractedSignatures[,-(sigToBeRemoved),drop = FALSE]
        }


        ## Convert extractedSignatures to ICAMS-formatted catalog.
        extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                                 region = "unknown",
                                                 catalog.type = "counts.signature")
        ## Output the signatures extracted
        write.catalog.function(extractedSignatures,
                               paste0(out.dir,"/extracted.signatures.csv"))
      }


    ## Step 5: Using hdp samples to attribute exposure counts.
    {

      ## Calculate the exposure probability of each signature(component) for each tumor sample(posterior sample corresponding to a dirichlet process node):
      ## This is the probability distribution of signatures(components) for all tumor samples(DP nodes)
      ## exposureProbs proves to be the normalized signature exposure all 100 tumor samples

      exposureProbs <- hdp::comp_dp_distn(mut_example_multi_extracted)$mean
      dim(exposureProbs)
      exposureProbs <- exposureProbs[3:dim(exposureProbs)[1],]
      rownames(exposureProbs) <- rownames(convSpectra)[1:dim(exposureProbs)[1]]
      ## Remove NA or NULL "hdp.0" signature in exposureProbs matrix.
      if(flagRemoveHDP0)
        exposureProbs <- exposureProbs[,-(sigToBeRemoved),drop = FALSE]
      ## Change signature names in exposureCounts
      ## from "0","1","2" to "hdp.0","hdp.1","hdp.2"
      colnames(exposureProbs) <- colnames(extractedSignatures)
      dim(exposureProbs)


      ## Calculate signature exposure counts from signature exposure probability
      ## Unnormalized exposure counts = Normalized exposure probability * Total mutation count in a sample
      sample_mutation_count <- apply(convSpectra,1,sum)

      exposureCounts <- matrix(nrow = dim(exposureProbs)[1], ncol = dim(exposureProbs)[2])
      dimnames(exposureCounts) <- dimnames(exposureProbs)
      for (sample in seq(1,dim(exposureProbs)[1]))
        exposureCounts[sample,] <- sample_mutation_count[[sample]] * exposureProbs[sample,]

      ## Change exposure count matrix to SynSigEval format.
      exposureCounts <- t(exposureCounts)

      ## Next, write the exposureCounts to a file
      ## Write attributed exposures into a SynSig formatted exposure file.
      WriteExposure(exposureCounts,
                    paste0(out.dir,"/attributed.exposures.csv"))
    }



      ## Save seeds and session information
      ## for better reproducibility
      capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
      write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
      write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

      ## Return a list of signatures and exposures
      invisible(list("signature" = extractedSignatures,
                     "exposure" = exposureCounts))
  }
