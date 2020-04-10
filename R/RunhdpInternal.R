#' Run hdp extraction and attribution on a spectra catalog file
#'
#' @param input.catalog A catalog of spectra catalog
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param CPU.cores Number of CPUs to use in running
#'    \code{\link[hdp]{hdp_posterior}}.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run hdp. Setting seed can make the
#' attribution of hdp repeatable.
#' Default: 1.
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
#' @param remove.noise Whether to remove noise signature "hdp.0"? In normal cases scenarios,
#' only few mutations will be assigned to noise signature.
#'
#' For result visualization and assessment of \code{hdp} package, select \code{TRUE};
#' for diagnostic purposes, select \code{FALSE}.
#'
#' @param overwrite If TRUE, overwrite existing output.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @param num.posterior Number of posterior sampling chains; can set to
#'   1 for testing.
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

RunhdpInternal <-
  function(input.catalog,
           out.dir,
           CPU.cores     = 1,
           seedNumber    = 1,
           K.guess,
           multi.types   = FALSE,
           remove.noise  = FALSE,
           overwrite     = FALSE,
           verbose       = TRUE,
           num.posterior = 4) {

    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  # To document the seed used
    RNGInUse <- RNGkind()      # To document the random number generator (RNG) used

    # Read in spectra data from input.catalog file
    spectra <- input.catalog

    ## CPU.cores specifies number of CPU cores to use.
    ## CPU.cores will be capped at 30.
    ## If CPU.cores is not specified, CPU.cores will
    ## be equal to the minimum of 30 or (total cores)/2

    ## convSpectra: convert the ICAMS-formatted spectra catalog
    ## into a matrix which HDP accepts:
    ## 1. Remove the catalog related attributes in convSpectra
    ## 2. Transpose the catalog
    convSpectra <- spectra

    # hdp get confused if there are extra attributes
    # or classes.
    class(convSpectra) <- "matrix"
    attr(convSpectra,"catalog.type") <- NULL
    attr(convSpectra,"region") <- NULL
    dimnames(convSpectra) <- dimnames(spectra) # check if really necessary
    convSpectra <- t(convSpectra)

    number.channels <- nrow(spectra)
    number.samples  <- ncol(spectra)

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    if (verbose) {
      message("number of Dirichlet process data clusters = ", K.guess)
    }
    ## Run hdp main program.
    ## Step 1: initialize hdp object
    {
      ## Allocate process index for hdp initialization.
      ## Each different index number refers to a dirichlet process
      ## for one tumor type.
      if(multi.types == FALSE){ ## All tumors belong to one tumor type (default)
        num.tumor.types <- 1
        process.index <- c(0,1,rep(2,number.samples))
      } else if(multi.types == TRUE){
        ## There are multiple tumors in the sample.
        ## Tumor type will be inferred by the string before "::" in the column names.
        ## e.g. Tumor type for "SA.Syn.Ovary-AdenoCA::S.500" would be "SA.Syn.Ovary-AdenoCA"
        tumor.types <- sapply(
          colnames(spectra),
          function(x) {strsplit(x,split = "::",fixed = T)[[1]][1]})
        num.tumor.types <- length(unique(tumor.types))
        ## 0 refers to top grandparent DP node. All signatures are drawn from this node.
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
      ## DP node of each level share two Dirichlet Hyperparameters:
      ## shape (alphaa) and rate (alphab).
      ## For mutational signature analysis purpose,
      ## alphaa and alphab are set as 1 for each level.
      alphaa <- rep(1,dp.levels)
      alphab <- rep(1,dp.levels)

      ## initialise hdp
      if (verbose) message("calling hdp_init")
      hdpObject <- hdp::hdp_init(ppindex = ppindex,
                                 cpindex = cpindex,
                                 hh = rep(1,number.channels),
                                 alphaa = alphaa,
                                 alphab = alphab)

      # num.process is the number of samples plus number cancer types plus 1 (grandparent)
      num.process <- hdp::numdp(hdpObject)

      if (verbose) message("calling hdp_setdata")
      hdpObject <- hdp::hdp_setdata(
        hdpObject,
        (1 + num.tumor.types + 1):num.process,
        convSpectra)

      if (verbose) message("calling dp_activate")

      # dp_activate calls hdp:::stirling, which when
      # called on a large number (e.g. > 200000), requires a great deal of
      # memory (e.g. estimated 2 Terabyte for stirling(200000)).
      hdpObject <- hdp::dp_activate(hdpObject,
                                    1:num.process,
                                    initcc = K.guess,
                                    seed = seedNumber)

      ## Release the occupied RAM by dp_activate
      gc()
    }

    ## Step 2: run 4 independent sampling chains
    {
      if(CPU.cores == 1){ ## debug
        ## Run four independent posterior sampling chains
        chlist <- vector("list", 4)	#4 is too much here!

        for (i in 1:4) {

          if (verbose) message("calling hdp_posterior ", i)
          chlist[[i]] <-
            hdp::hdp_posterior(
              hdpObject,
              # The remaining values, except seed, are from the vignette; there
              # are no defaults.
              burnin = 4000,
              n      = 50,
              space  = 50,
              cpiter = 3,
              # Must choose a different seed for each of the 4 chains:
              seed   = (seedNumber + i * 10^6) %% (10^7) )
        }
      } else {
        f_posterior <- function(seed,hdpObject) {
          if (verbose) message("calling hdp_posterior")
            retval <- hdp::hdp_posterior(
              hdpObject,
              # The remaining values, except seed, are from the vignette; there
              # are no defaults.
              burnin = 4000,
              n      = 50,
              space  = 50,
              cpiter = 3,
              # Must choose a different seed for each of the 4 chains:
              seed   = seed %% (10^7) )
            return(retval)
        }

        chlist <- parallel::mclapply(
          X = seedNumber + 10^6 * 1:4,
          FUN = f_posterior,
          hdpObject = hdpObject,
          mc.cores = CPU.cores
        )
      }

      ## Generate the original multi_chain for the sample
      if (verbose) message("calling hdp_multi_chain")
      mut_example_multi <- hdp::hdp_multi_chain(chlist)
    }


    ## Step 3: Plot the diagnostics of sampling chains.
    {
      ## Plotting using hdp functions
      ## Plotting device on the server does not work
      ## Need to plot the file into a pdf


      if (verbose) message("plotting to original_sample.pdf")
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

      if (verbose) message("calling hdp_extract_components")
      ## Extract components(here is signature) with cosine.merge = 0.90 (default)
      mut_example_multi_extracted <- hdp::hdp_extract_components(mut_example_multi)
      mut_example_multi_extracted


      ## Generate a pdf for mut_example_multi_extracted
      {
        if (verbose) message("plotting to signature_hdp_embedded_func.pdf")
        grDevices::pdf(
          file = paste0(out.dir,"/signature_hdp_embedded_func.pdf"))
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
      if (verbose) message("calling hdp::comp_categ_distn")
      ## Calculate the mutation composition in each signature:
      extractedSignatures <-
        hdp::comp_categ_distn(mut_example_multi_extracted)$mean
      dim(extractedSignatures)
      ## Add base context for extractedSignatures
      extractedSignatures <- t(extractedSignatures)
      rownames(extractedSignatures) <- rownames(spectra)
      ## Change signature names in extractedSignatures
      ## from "0","1","2" to "hdp.0","hdp.1","hdp.2"
      colnames(extractedSignatures) <-
        paste("hdp", colnames(extractedSignatures), sep = ".")
      ## Remove "hdp.0" (noise signature) if remove.noise == TRUE
      flagRemoveHDP0 <- FALSE
      if(FALSE){ ## debug
        ## Remove "hdp.0" (noise signature) if it is a null signature or NA signature
        if(is.null(extractedSignatures[1,"hdp.0"]) | is.na(extractedSignatures[1,"hdp.0"]))
          flagRemoveHDP0 <- TRUE
      }
      if(flagRemoveHDP0){
        sigToBeRemoved <- which(colnames(extractedSignatures) == "hdp.0")
        if(length(sigToBeRemoved) == 1)
          extractedSignatures <- extractedSignatures[,-(sigToBeRemoved),drop = FALSE]
      }


      ## Convert extractedSignatures to ICAMS-formatted catalog.
      extractedSignatures <- ICAMS::as.catalog(extractedSignatures,
                                               region = "unknown",
                                               catalog.type = "counts.signature")

      if (verbose) message("calling ICAMS::WriteCatalog")
      ICAMS::WriteCatalog(extractedSignatures,
                          paste0(out.dir,"/extracted.signatures.csv"))
    }


    ## Step 5: Using hdp samples to attribute exposure counts.
    {

      ## Calculate the exposure probability of each signature(component) for each tumor sample(posterior sample corresponding to a dirichlet process node):
      ## This is the probability distribution of signatures(components) for all tumor samples(DP nodes)
      ## exposureProbs proves to be the normalized signature exposure all tumor samples

      if (verbose) message("Calling hdp::comp_dp_distn to generate exposure probability")
      exposureProbs <- hdp::comp_dp_distn(mut_example_multi_extracted)$mean
      dim(exposureProbs)
      exposureProbs <- exposureProbs[(num.tumor.types + 2):dim(exposureProbs)[1],]
      rownames(exposureProbs) <- rownames(convSpectra)[1:dim(exposureProbs)[1]]
      ## Remove NA or NULL "hdp.0" signature in exposureProbs matrix.
      if(flagRemoveHDP0)
        exposureProbs <- exposureProbs[,-(sigToBeRemoved),drop = FALSE]
      ## Change signature names in exposureCounts
      ## from "0","1","2" to "hdp.0","hdp.1","hdp.2"
      colnames(exposureProbs) <- colnames(extractedSignatures)
      dim(exposureProbs)
      ## Transpose exposureProbs so that it conforms to SynSigEval format.
      exposureProbs <- t(exposureProbs)


      ## Calculate signature exposure counts from signature exposure probability
      ## Unnormalized exposure counts = Normalized exposure probability * Total mutation count in a sample
      sample_mutation_count <- apply(convSpectra,1,sum)

      exposureCounts <- matrix(nrow = nrow(exposureProbs), ncol = ncol(exposureProbs))
      dimnames(exposureCounts) <- dimnames(exposureProbs)
      for (sample in seq(1,ncol(exposureProbs)))
        exposureCounts[,sample] <- sample_mutation_count[[sample]] * exposureProbs[,sample]

      if (verbose) message("Calling WriteExposure() to write exposure counts")
      WriteExposure(exposureProbs,
                    paste0(out.dir,"/exposure.probs.csv"))
      WriteExposure(exposureCounts,
                    paste0(out.dir,"/inferred.exposures.csv"))
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
