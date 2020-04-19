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
           posterior.verbosity = 0) {

    if (verbose) message("Reading catalog file ", input.catalog.file)
    spectra <- ICAMS::ReadCatalog(input.catalog.file, strict = FALSE)
    if (test.only > 0) spectra <- spectra[ , 1:test.only]

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
      message("Using existing out.dir ", out.dir)
    } else {
      dir.create(out.dir, recursive = T)
      message("Created new out.dir", out.dir)
    }

    retval <- RunhdpInternal(
      input.catalog       = spectra,
      CPU.cores           = CPU.cores,
      seedNumber          = seedNumber,
      K.guess             = K.guess,
      multi.types         = multi.types,
      verbose             = verbose,
      num.posterior       = num.posterior,
      posterior.verbosity = posterior.verbosity
    )

    # Plot the diagnostics of sampling chains.
    if (verbose) message("Writing sampling chain information")

    ## Draw the DP oscillation plot for mut_example_multi(original_sample)
    chains <- hdp::chains(retval$multi.chains) # Gibbs sampling chains
    grDevices::pdf(file = paste0(out.dir,"/original_sample.pdf"))
    graphics::par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
    p1 <- lapply(chains, hdp::plot_lik, bty="L", start=500)
    p2 <- lapply(chains, hdp::plot_numcluster, bty="L")
    grDevices::dev.off()


    ex.chains <- hdp::chains(retval$ex.multi.chains) # extracted chains
    grDevices::pdf(
      file = paste0(out.dir,"/signature_hdp_embedded_func.pdf"))
    ## Draw the DP oscillation plot for mut_example_multi_extracted
    graphics::par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
    p1 <- lapply(ex.chains, hdp::plot_lik, bty="L", start=500)
    p2 <- lapply(ex.chains, hdp::plot_numcluster, bty="L")
    ## Draw the computation size plot
    graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
    hdp::plot_comp_size(mut_example_multi_extracted, bty="L")
    grDevices::dev.off()

    if (verbose) message("Writing signatures")
    extractedSignatures <- ICAMS::as.catalog(retval$signture,
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
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt"))
    write(x = retval$seedInUse, file = paste0(out.dir,"/seedInUse.txt"))
    write(x = retval$RNGInUse, file = paste0(out.dir,"/RNGInUse.txt"))

    invisible(retval)
  }
