context("Test summarizing extraction results in one directory")

test_that("signature.analyzer.sa.sa.COMPOSITE.out", {
  load("signature.analyzer.sa.sa.COMPOSITE.out")
  expect_equal(
    SummarizeSigOneSASubdir(
      run.dir = "./test.sa.COMPOSITE/sa.results/",
      ground.truth.exposure.dir = "./test.sa.COMPOSITE/",
      which.run = "1.run",
      overwrite = TRUE),
    signature.analyzer.sa.sa.COMPOSITE.out)
  if (TRUE) # Set this to FALSE to inspect the file and directory output
  {
    # Clean up
    unlink("./test.sa.sa.COMPOSITE/sa.results/summary",
           recursive = TRUE, force = TRUE)
  }
})

test_that("signature.analyzer.sa.sa.96.out", {
  load("signature.analyzer.sa.sa.96.out.Rdata")
  expect_equal(
    SummarizeSigOneSASubdir(
      run.dir = "./test.sa.sa.96/sa.results/",
      ground.truth.exposure.dir = "./test.sa.sa.96/",
      which.run = "1.run.sa.sa.96",
      overwrite = TRUE),
    signature.analyzer.sa.sa.96.out)
  if (TRUE) # Set this to FALSE to inspect the file and directory output
  {
    # Clean up
    file.to.unlink <- "./test.sa.sa.96/sa.results/summary"
    if (!file.exists(file.to.unlink)) cat("Coding error,", file.to.unlink,
                                          "does not exist\n")
    res <- unlink(file.to.unlink,
           recursive = TRUE, force = TRUE)
    if (res != 0) cat("Failed to unlink", file.to.unlink, "\n")
  }
})

test_that("SummarizeSigOneSigProExtractorSubdir (SBS96)", {
  load("./test.sigprofiler.sp.sp.out.Rdata")

  # Warning, do not change t1 to a longer name,
  # or file2 below will be too long for portable zip'ing.
  tdir.res <- "./"
  expect_equal(
    SummarizeSigOneSigProExtractorSubdir(
      run.dir = tdir.res,
      overwrite = T,
      ground.truth.exposure.dir = "tests/testthat/sp.96.ground.truth/"),
    sigprofiler.sp.sp.out)
  if (TRUE) # Set this to FALSE to inspect the file and directory output
  {
    # Clean up
    file1 <- paste0(tdir.res, "summary")
    if (!file.exists(file1)) cat("coding error, wrong file\n")
    res <- unlink(file1, recursive = TRUE, force = TRUE)
    if (res != 0) cat("Failed to unlink ", file1, "\n")
    file2 <- paste0(tdir.res,
                    "SBS96/Suggested_Solution/",
                    "De_Novo_Solution/signatures.PCAWG.format.csv")
    res <- unlink(file2)
    if (res != 0) cat("Failed to unlink ", file2, "\n")
    file3 <- paste0(tdir.res,
                    "SBS96/Suggested_Solution/",
                    "De_Novo_Solution/attributed.exposures.csv")
    res <- unlink(file3)
    if (res != 0) cat("Failed to unlink ", file3, "\n")
  }
})
