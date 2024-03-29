% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SignatureAnalyzerInteraction.R
\name{SignatureAnalyzerOneRun}
\alias{SignatureAnalyzerOneRun}
\title{Source SignatureAnalyzer and run it once on a single data set
and put results in specified location.}
\usage{
SignatureAnalyzerOneRun(
  signatureanalyzer.code.dir,
  input.catalog,
  out.dir,
  seedNumber = NULL,
  input.exposures = NULL,
  maxK = 30,
  tol = 1e-07,
  test.only = FALSE,
  delete.tmp.files = TRUE,
  verbose = 0,
  overwrite = FALSE
)
}
\arguments{
\item{signatureanalyzer.code.dir}{The directory holding the
SignatureAnalyzer code.}

\item{input.catalog}{File containing input catalog.  Columns are
samples (tumors), rows are signatures.  SignatureAnalyzer does
not care about the row names (I think) TODO(Steve): check this.}

\item{out.dir}{Directory that will be created for the output;
abort if it already exits.  Log files will be in
\code{paste0(out.dir, "/tmp")}.}

\item{seedNumber}{Specify the pseudo-random seed number
used to run SignatureAnalyzer. Setting seed can make the
attribution of SignatureAnalyzer repeatable.
If NULL, this function will not specify seed number.
Default: NULL.}

\item{input.exposures}{A file with the synthetic exposures used to generate
\code{input.catalog}; if provided here,
this is copied over to the output directory
for downstream analysis.}

\item{maxK}{The maximum number of signatures to consider
extracting.}

\item{tol}{Controls when SignatureAnalyzer will terminate
its search; \code{tol} was 1.e-05 for the PCAWG7 analysis.}

\item{test.only}{If TRUE, only analyze the first 10 columns
read in from \code{input.catalog}.}

\item{delete.tmp.files}{If TRUE delete the many temporary
files generated by SignatureAnalyzer.}

\item{verbose}{If \code{TRUE}, then print various messages.}

\item{overwrite}{If TRUE, overwrite existing output}
}
\value{
A list with the following elements:
\enumerate{
 \item \code{signatures.W} The raw signature matrix, *including*
       columns of all zeros.
 \item \code{exposures.H} The raw exposure matrix, *excluding*
       rows of all zeros. The matrix
       product of the non-zero columns of \code{signatures.w}
       and \code{exposures.H} approximates the input spectrum
       matrix.
 \item \code{likelihood} The likelihood as returned by
 SignatureAnalyzer.
 \item \code{evidence} -1 * the posterior probability
 as returned by SignatureAnalyzer.
 \item \code{relevance} One for each column of the \code{signatures.W},
 as returned by SignatureAnalyzer.
 \item \code{error} A measure of reconstruction error (?) as
 returned by SignatureAnalyzer
 \item \code{normalized.sigs} The non-0 columns of \code{signatures.W}
 normalized so that each column sum is 1.
}
}
\description{
Source SignatureAnalyzer and run it once on a single data set
and put results in specified location.
}
\details{
Creates several files in \code{out.dir}:
\enumerate{
\item \code{sa.output.sigs.csv} Normalized signatures (no all-0 signatures,
column sums all 0)
\item \code{sa.output.raw.exp.csv} Raw exposures (attributions)
\item \code{sa.output.exp.csv} Same as \code{sa.output.raw.exp.csv}
\item \code{sa.output.other.data.csv}, contains a summary of important
information, including the number of signatures extracted.
\item \code{input.syn.exp.csv} Optional, a copy of \code{input.exposures},
if it was provided.
}
}
