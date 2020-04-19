
test_that("Runhdp2", {
  input.catalog.file <-
    system.file(
      "tests/SBS96.ground.truth/ground.truth.syn.catalog.csv",
      package = "SynSigRun",
      mustWork = TRUE)

load(
  system.file("tests/RunhdpInternal.testdata/t2.out.Rdata",
              package = "SynSigRun",
              mustWork = TRUE))

out.dir.root <- system.file("tests",
                       package = "SynSigRun",
                       mustWork = TRUE)

retvalx <- Runhdp2(
  input.catalog.file = input.catalog.file,
  out.dir            = file.path(out.dir.root, "test_Runhdp2_out_dir"),
  CPU.cores          = 1,
  seedNumber         = 44,
  K.guess            = 5,
  multi.types        = FALSE,
  overwrite          = TRUE,
  verbose            = TRUE,
  num.posterior      = 1,
  test.only          = 15
)

foo <- t2.out$signature
class(foo) <- "matrix"
testthat::expect_equal(retvalx$signature, foo, check.attributes = FALSE)
testthat::expect_equal(retvalx$exposure, t2.out$exposure)
})
