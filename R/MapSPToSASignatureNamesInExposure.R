#' With the signatures represented in a matrix of exposures, find the nearest
#' SignatureAnalyzer exposure.
#'
#' @param sp.exposures The exposures
#'
#' @param sa.sig.names.to.consider A subset of the colnames
#' of \code{\link{sa.96.sigs}}
#'
#' @details IMPORTANT: uses the package global
#' variables \code{\link{sa.96.sigs}}
#' and \code{\link{sp.sigs}}.
#'
#' @return A list with
#'
#' \enumerate{
#' \item \code{exp2} Copy of \code{sp.exposures} with the
#' rownames(signature names) updated according to the
#' match.
#'
#' \item \code{sp.to.sa.sig.match}
#'
#' \item \code{sa.to.sp.sig.match} Best matches in the opposite direction
#' }
#'
#' @export

MapSPToSASignatureNamesInExposure <-
  function(sp.exposures,
           sa.sig.names.to.consider = colnames(sa.96.sigs))
  {
    SP.SA.mappings <-
      MatchSigs2Directions(
        sp.sigs[ , rownames(sp.exposures)],
        sa.96.sigs[ , sa.sig.names.to.consider, drop = FALSE])
    new.syn.exp.rownames <- SP.SA.mappings$match1[rownames(sp.exposures), "to"]

    if(any(duplicated(new.syn.exp.rownames))) {
      print(SP.SA.mappings$match1)
      stop("\n\nCannot handle this situation:\n",
           "> 1 SigProfiler signature maps to\n",
           "the same SignatureAnalyzer signature\n\n")
      # Presumably this would die a few lines below when we
      # set rownames(exp2)
    }

    exp2 <- sp.exposures
    rownames(exp2) <- new.syn.exp.rownames
    return(list(exp2               = exp2,
                sp.to.sa.sig.match = SP.SA.mappings$match1,
                sa.to.sp.sig.match = SP.SA.mappings$match2))
  }
