#################################################################################################
###### load prerequisites
#################################################################################################
import sys,subprocess,os
import random ## Required to designate random seed
import gc ## Garbage collector


#### Naming the seeds
seedNumbers = (1, 691, 1999, 3511, 8009,
    9902, 10163, 10509, 14476, 20897,
    27847, 34637, 49081, 75679, 103333,
    145879, 200437, 310111, 528401, 1076753)


#### Naming the datasets for cycling
slopes = (0.1,0.5,1,2,10)
Rsqs = (0.1,0.2,0.3,0.6)
datasetNames = tuple()
for slope in slopes:
    for Rsq in Rsqs:
        datasetNames = datasetNames + ("S."+str(slope)+".Rsq."+str(Rsq),)



#### Set GLOBAL working directory
workingDir = os.getcwd()


#### Import sigproextractor: de novo extraction+attribution module
import sigproextractor
from sigproextractor import sigpro as sig
toolName = "sigproextractor"
toolAcronym = "sp"

#### Import sigproSS: attribution-only module
import sigproSS
from sigproSS import spss, spss_pcwag

## Set global parameters
minK = 1
maxK = 10
totalIterations = 100
numcpus = 10

## Run sigproextractor for each dataset for each seed.
for seedNumber in seedNumbers:
    for datasetName in datasetNames:
        inputDir = "/".join([datasetName,"sp.sp/"])
        inputCatalog = "/".join([inputDir,"ground.truth.syn.catalog.csv"])
        outputDir = "/".join([inputDir,"ExtrAttr/sigproextractor.results/seed."+str(seedNumber)])
        ## Set seed
        random.seed(seedNumber) ## Set seed to seedNumber
        ## extract signatures and attribute exposures
        sig.sigProfilerExtractor("csv", outputDir, inputCatalog, refgen="GRCh37", 
        startProcess = minK, endProcess = maxK, totalIterations = numiters, 
        cpu = numcpus, hierarchy = False, mtype = ["default"], exome = False)


#### Run attribution-only script.
for seedNumber in seedNumbers:
    for datasetName in datasetNames:
        inputDir = "/".join([datasetName,"sp.sp/"])
        inputCatalog = "/".join([inputDir,"ground.truth.syn.catalog.csv"])
        inputSigs = "/".join([inputDir,"ground.truth.syn.sigs.csv"])
        outputDir = "/".join([inputDir,"Attr/sigproextractor.results/seed."+str(seedNumber)])
        spss_pcwag.single_sample_pcwag(samples=inputCatalog, output=outputDir, 
        sigbase=inputSigs, n_cpu = numcpus)
		
		





