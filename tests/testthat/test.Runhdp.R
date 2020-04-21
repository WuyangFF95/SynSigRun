
input.catalog <-
  ICAMS::ReadCatalog("tests/SBS96.ground.truth/ground.truth.syn.catalog.csv")

stir.closure <- hdp::make.stirling()

retval <- SynSigRun:::RunhdpInternal(
  input.catalog = input.catalog[ , 1:15],
  out.dir = "tests/test_Hdprun_out_dir",
  CPU.cores = 1,
  seedNumber = 44,
  K.guess = 5,
  multi.types = FALSE,
  remove.noise = FALSE,
  overwrite = TRUE,
  verbose       = TRUE,
  num.posterior = 1
)
