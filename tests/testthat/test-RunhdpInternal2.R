
test_that("RunhdpInternal2", {

  input.catalog <-
    ICAMS::ReadCatalog(
      system.file(
        "tests/SBS96.ground.truth/ground.truth.syn.catalog.csv",
        package = "SynSigRun",
        mustWork = TRUE))

  load(
    system.file("tests/RunhdpInternal.testdata/test.RunhdpInternal.2.Rdata",
                package = "SynSigRun",
                mustWork = TRUE))

  retvalx <- RunhdpInternal(
    input.catalog = input.catalog[1:10 , 1:15],
    CPU.cores = 1,
    seedNumber = 44,
    K.guess = 5,
    multi.types = FALSE,
    verbose       = TRUE,
    num.posterior = 1
  )

  testthat::expect_equal(retvalx$signature, test.RunhdpInternal.2$signature)
  testthat::expect_equal(retvalx$exposure, test.RunhdpInternal.2$exposure)
})
