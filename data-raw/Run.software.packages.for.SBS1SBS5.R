
## Load required packages
library(ICAMS)
library(SynSigEval)



## Specify slopes and Rsqs for the datasets
slopes <- c(0.1,0.5,1,2,10)
Rsqs <- c(0.1,0.2,0.3,0.6)
datasetNames <- c()
for(slope in slopes)
  for(Rsq in Rsqs)
    datasetNames <- c(datasetNames, paste0("S.",slope,".Rsq.",Rsq))

## Specify 20 seeds used in software running
seedsInUse <- c(1, 691, 1999, 3511, 8009,
                9902, 10163, 10509, 14476,
                20897, 27847, 34637, 49081, 75679,
                103333, 145879, 200437, 310111, 528401)



## Run R-based Extraction and attribution packages.
## sigproextractor (Python package) and MultiModalMuSig (Julia package)
## needs to be run with external script.
extrAttrPackages <- c("signeR","sigfit","MutationalPatterns",
                      "SomaticSignatures","hdp")
for(extrAttrPackage in extrAttrPackages){
  func <- get(paste0("Run",extrAttrPackage))
  for(seedInUse in seedsInUse){
    for(datasetName in datasetNames){
      func(input.catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
           read.catalog.function = ICAMS::ReadCatalog,
           write.catalog.function = ICAMS::WriteCatalog,
           out.dir = paste0(datasetName,"/sp.sp/ExtrAttr/",extrAttrPackage,".results/seed.",seedInUse),
           seedNumber = seedInUse,
           K.range = c(1,10),
           overwrite = T)
    }
  }
}


## Run SignatureAnalyzer
for(seedInUse in seedsInUse){
  for(datasetName in datasetNames){
    SAMultiRunOneCatalog(
      num.runs = 20,
      signatureanalyzer.code.dir = "~/softwares/SynSigEval/data-raw/SignatureAnalzyer.052418/",
      input.catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
      out.dir = paste0(datasetName,"/sp.sp/ExtrAttr/SignatureAnalyzer.results/seed.",seedInUse),
      maxK = 10,
      overwrite = T,
      mc.cores = 10,
      seed = seedInUse)
  }
}


## Run Attribution-only packages (deconstructSigs, YAPSA, SignatureEstimation)
## Note: SignatureEstimation has two methods:
attrOnlyPackages <- c("deconstructSigs","YAPSA")

for(attrOnlyPackage in attrOnlyPackages){
  func <- get(paste0("Run",attrOnlyPackage,"AttributeOnly"))
  for(datasetName in datasetNames){
    for(seedInUse in seedsInUse){
    func(input.catalog = paste0(datasetName,"/sp.sp/ground.truth.syn.catalog.csv"),
         gt.sigs.file = paste0(datasetName,"/sp.sp/ground.truth.syn.sigs.csv"),
         read.catalog.function = ICAMS::ReadCatalog,
         out.dir = paste0(datasetName,"/sp.sp/Attr/",attrOnlyPackage,".results/seed.",seedInUse),
         seedNumber = seedInUse,
         overwrite = T)
    }
  }
}
