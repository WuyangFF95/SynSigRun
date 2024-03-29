#' Install YAPSA package from Bioconductor.
#'
#' @keywords internal
InstallYAPSA <- function(){
  message("Installing YAPSA from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    utils::install.packages("BiocManager")
  BiocManager::install("YAPSA")
}


#' Run YAPSA attribution on a spectra catalog file
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
#' used to run YAPSA. Setting seed can make the
#' attribution of YAPSA repeatable.
#' Default: 1.
#'
#' @param signature.cutoff A numeric vector of values less than 1.
#' Signatures from within W with an overall exposure
#' less than the respective value in \code{in_cutoff_vector}
#' will be discarded.
#' Default: vector length of number of sigs with all zeros
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The inferred exposure of \code{YAPSA}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export
#'
RunYAPSAAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           out.dir,
           seedNumber = 1,
           signature.cutoff = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    # Install YAPSA from Bioconductor, if failed to be loaded
    if (!requireNamespace("YAPSA", quietly = TRUE)) {
      InstallYAPSA()
    }

    # Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  # Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() # Save the random number generator (RNG) used


    # Read in spectra data from input.catalog file
    # spectra: spectra data.frame in ICAMS format
    spectra <- ICAMS::ReadCatalog(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]


    # Read in ground-truth signatures
    # gtSignatures: signature data.frame in ICAMS format
    gtSignatures <- ICAMS::ReadCatalog(gt.sigs.file)

    # Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    # If signature.cutoff is NULL (by default),
    # set it to all zeros of length K (number of signatures)
    if(is.null(signature.cutoff))
      signature.cutoff = rep(0,times = ncol(gtSignatures))

    # Derive exposure count attribution results.

    # Known signature matrix
    in_signatures_df <- gtSignatures
    class(in_signatures_df) <- "matrix"
    attr(in_signatures_df,"catalog.type") <- NULL
    attr(in_signatures_df,"region") <- NULL

    ### Tumor spectra matrix and related parameters
    in_mutation_catalogue_df <- spectra # Converted spectra matrix
    size <- colSums(in_mutation_catalogue_df) # Total mutation count of each spectrum
    class(in_mutation_catalogue_df) <- "matrix"
    attr(in_mutation_catalogue_df,"catalog.type") <- NULL
    attr(in_mutation_catalogue_df,"region") <- NULL

    # Plotting parameter - maximum height in the plot
    ymax <- rep(0.4,ncol(in_mutation_catalogue_df))
    names(ymax) <- colnames(in_mutation_catalogue_df)

    # Using Linear Combination Decomposition to attribute exposures
    # YAPSA::LCD() is not recommended. The author recommended YAPSA::LCD_complex_cutoff(),
    # which is a wrapper of it.
    # YAPSA also supports different presence cutoff for different signatures,
    # this is done by providing different values of cutoff in LCD_complex_cutoff function.
    # Authors suggest to use YAPSA::LCD_complex_cutoff() rather than YAPSA::LCD() in most cases.
    LCD_complex_object <- YAPSA::LCD_complex_cutoff(in_mutation_catalogue_df,
                                                    in_signatures_df,
                                                    in_cutoff_vector = signature.cutoff, # If there are 2 signatures in the spectra,
                                                    # you must provide a
                                                    in_rescale = TRUE)  # Rescale signature exposures so that the sum of exposure for each tumor
    # equals to the exposure sum in original spectra
    # This prevents the difference between original spectra and observed spectra
    class(LCD_complex_object) # [1] "list"
    names(LCD_complex_object) # For detail, see YAPSA user manual
    ##[1] "exposures"                   "norm_exposures"
    ##[3] "signatures"                  "choice"
    ##[5] "order"                       "residual_catalogue"
    ##[7] "rss"                         "cosDist_fit_orig_per_matrix"
    ##[9] "cosDist_fit_orig_per_col"    "sum_ind"
    ##[11] "out_sig_ind_df"              "aggregate_exposures_list"

    # Exposures generated by LCD_complex_object()
    # does not equal to exposures generated by LCD()
    # Because by default, LCD_complex_object normalizes the counts.
    if(FALSE){
      dim(LCD_complex_object$exposures) == dim(LCD_object) # [1] TRUE
      LCD_complex_object$exposures == LCD_object # [1] FALSE
    }

    # For each tumor spectrum, $exposures (the exposure counts inferred by LCD_complex_object())
    # sums up to the total mutation counts in 500 tumors in the dataset.
    # But $norm_exposures (relative exposure probs inferred by LCD_complex_object())
    # sums up to number of tumors only.
    sum(LCD_complex_object$exposures) == sum(spectra) # [1] TRUE
    sum(LCD_complex_object$norm_exposures) # [1] (Number of tumors in spectra)

    # For each tumor spectrum, sum of normalized inferred exposures by LCD_complex_cutoff()
    # does not equal to the sum of ground-truth exposures.
    all( colSums(LCD_complex_object$norm_exposures) == colSums(spectra) ) # [1] FALSE

    # Export inferred exposure probs
    LCD_exposure_prob <- LCD_complex_object$norm_exposures
    # Export inferred exposure counts
    exposureCounts <- LCD_complex_object$exposures # Export exposure probs

    # Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    # Write inferred exposures into a SynSig formatted exposure file.
    SynSigGen::WriteExposure(exposureCounts,
                  paste0(out.dir,"/inferred.exposures.csv"))

    # Save seeds and session information
    # for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) # Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) # Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) # Save seed in use to a text file

    # Return the exposures inferred, invisibly
    invisible(exposureCounts)
  }
