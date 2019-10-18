#' @keywords internal
MustCreateDir <- function(dir, overwrite = FALSE) {
  if (dir.exists(dir)) {
    if (!overwrite) stop(dir, " exists and overwrite is FALSE")
    return(NULL)
  }
  if (!dir.create(dir, recursive = TRUE)) {
    stop("Unable to create dir ", dir )
  }
  return(NULL)
}
