
test_that("RunhdpInternal", {
input.catalog <-
  ICAMS::ReadCatalog(
    system.file(
      "tests/SBS96.ground.truth/ground.truth.syn.catalog.csv",
      package = "SynSigRun",
      mustWork = TRUE))

load(
  system.file("tests/RunhdpInternal.testdata/t2.out.Rdata",
              package = "SynSigRun",
              mustWork = TRUE))

retvalx <- RunhdpInternal(
  input.catalog = input.catalog[ , 1:15],
  # out.dir = "tests/test_RunhdpInternal_out_dir",
  CPU.cores = 1,
  seedNumber = 44,
  K.guess = 5,
  multi.types = FALSE,
  remove.noise = FALSE,
  overwrite = TRUE,
  verbose       = TRUE,
  num.posterior = 1
)

testthat::expect_equal(retvalx, t2.out)
})
