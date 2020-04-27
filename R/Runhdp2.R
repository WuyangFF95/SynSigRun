#' Run hdp extraction and attribution on a spectra catalog file
#'
#' @inheritParams RunhdpInternal
#'
#' @param input.catalog.file File containing a spectra catalog
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param remove.noise Deprecated; ignored
#'
#' For result visualization and assessment of \code{hdp} package, select \code{TRUE};
#' for diagnostic purposes, select \code{FALSE}.
#'
#' @param test.only If > 0, only analyze the first \code{test.only} columns
#'  in \code{input.catalog.file}.
#'
#' @return The same list as returned by \code{\link{RunhdpInternal}}.
#'
#' @details Creates several files in \code{out.dir}. These are:
#'  TODO(Steve): list the files
#'
#' @export

Runhdp2 <-
  function(input.catalog.file,
           out.dir,
           CPU.cores           = 1,
           seedNumber          = 1,
           K.guess,
           multi.types         = FALSE,
           remove.noise        = FALSE,
           test.only           = 0,
           overwrite           = FALSE,
           verbose             = TRUE,
           num.posterior       = 4,
           post.burnin         = 4000,
           post.n              = 50,
           post.space          = 50,
           post.cpiter         = 3,
           post.verbosity      = 0,
           cos.merge           = 0.9,
           min.sample          = 1) {

    if (verbose) message("Reading input catalog file ", input.catalog.file)
    spectra <- ICAMS::ReadCatalog(input.catalog.file, strict = FALSE)
    if (test.only > 0) spectra <- spectra[ , 1:test.only]

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
      if (verbose) message("Using existing out.dir ", out.dir)
    } else {
      dir.create(out.dir, recursive = T)
      if (verbose) message("Created new out.dir", out.dir)
    }

    retval <- RunhdpInternal(
      input.catalog   = spectra,
      CPU.cores       = CPU.cores,
      seedNumber      = seedNumber,
      K.guess         = K.guess,
      multi.types     = multi.types,
      num.posterior   = num.posterior,
      verbose         = verbose,
      post.burnin     = post.burnin,
      post.n          = post.n,
      post.space      = post.space,
      post.cpiter     = post.cpiter,
      post.verbosity  = post.verbosity,
      cos.merge       = cos.merge,
      min.sample      = min.sample
    ) # 14 Arguments

    multi <- retval[["multi.chains"]]
    chains <- hdp::chains(multi) # Gibbs sampling chains
    # cat("Runhdp2, class of chains ", class(chains), "\n")
    # cat("Runhdp2, class of multi ",  class(multi),  "\n")

    # Plot the diagnostics of sampling chains.
    if (verbose) message("Writing hdp.diagnostics.pdf")
    grDevices::pdf(file = paste0(out.dir,"/hdp.diagnostics.pdf"))
    graphics::par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
    p1 <- lapply(chains, hdp::plot_lik, bty="L")
    p2 <- lapply(chains, hdp::plot_numcluster, bty="L")

    graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
    hdp::plot_comp_size(multi, bty="L")

    graphics::par(mfrow=c(8, 1), mar = c(1, 1, 1, 1))
    hdp::plot_comp_distn(multi)

    # TODO, need argument dpinices and col_comp;
    # Need to return the hdp object (perhaps) from RunhdpInternal
    # to get the required values.
    # hdp::plot_dp_comp_exposure(multi)

    grDevices::dev.off()

    if (verbose) message("Writing signatures")
    extractedSignatures <- ICAMS::as.catalog(retval$signature,
                                             region       = "unknown",
                                             catalog.type = "counts.signature")
    ICAMS::WriteCatalog(extractedSignatures,
                        paste0(out.dir,"/extracted.signatures.csv"))

    if (verbose) message("Writing exposures")
    WriteExposure(retval$exposure.p,
                  paste0(out.dir,"/exposure.probs.csv"))
    WriteExposure(retval$exposure,
                  paste0(out.dir,"/inferred.exposures.csv"))

    if (verbose) message("Writting additonal information")
    utils::capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt"))

    invisible(retval)
  }
