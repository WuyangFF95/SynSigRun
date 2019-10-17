# Script to run SignatureAnalyzer on the exome-subset of PCAWG
# ID mutational spectra

library(ICAMS)
library(SynSigEval)

num.runs                   <- 20 # 2 for debugging
# signatureanalyzer.code.dir <- "SignatureAnalzyer.052418" # for debugging on Laptop
signatureanalyzer.code.dir <- "/home/gmssgr/bin/SignatureAnalzyer.052418/" # for monster
input.catalog              <- "sa.ID.exome.subset/pcawg-as-exome-ID.csv"
out.dir                    <- "sa.ID.exome.subset/"
maxK                       <- 30
test.only                  <- FALSE
delete.tmp.files           <- TRUE
overwrite                  <- TRUE
mc.cores                   <- 20 # Will be overidden and set to 1 on Windows
verbose                    <- TRUE

# TODO(Steve): move RNGking and set.seed into SAMultiRunOneCatalog
cat("\nRunning, maxK is ", maxK, "\n\n", sep = "")
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(897)

sa.ID.exome.subset.res <-
  SAMultiRunOneCatalog(
    num.runs                   = num.runs,
    signatureanalyzer.code.dir = signatureanalyzer.code.dir,
    input.catalog              = input.catalog,
    out.dir                    = out.dir,
    maxK                       = maxK,
    tol                        = 1e-7,
    test.only                  = test.only,
    delete.tmp.files           = delete.tmp.files,
    overwrite                  = overwrite,
    mc.cores                   = mc.cores,
    verbose                    = verbose)
