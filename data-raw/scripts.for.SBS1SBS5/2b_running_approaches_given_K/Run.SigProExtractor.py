
# SigProExtractor accepts ICAMS-formatted catalogs.
# However, SigProExtractor results need to be summarized by
# its own summary function SynSigEval::SummarizeSigOneSigProExtractorSubdir()


#################################################################################################
###### load prerequisites
#################################################################################################
import sys,subprocess,os
import os.path
import random ## Required to designate random seed
import gc ## Garbage collector

#### Read old working directory
oldWorkingDir = os.getcwd()


# Set working directory to "<SynSigRun Home>/data-raw/scripts.for.SBS1SBS5"
# before running this script.
#### Set GLOBAL working directory
os.chdir("../")
workingDir = os.getcwd()


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

## Set global parameters
exactK = 2
numiters = 1000
numcpus = 10

## Run sigproextractor for each dataset for each seed.
for seedNumber in seedNumbers:
    for datasetName in datasetNames:
        inputDir = "/".join([datasetName,"sp.sp/"])
        inputCatalog = "/".join([inputDir,"ground.truth.syn.catalog.csv"])
        outputDir = "/".join([inputDir,"ExtrAttrExact/SigProExtractor.results/seed."+str(seedNumber)])
        ## Set seed
        random.seed(seedNumber) ## Set seed to seedNumber
        ## extract signatures and attribute exposures
        print("\n\n#####################")
        print("Start running catalog "+str(inputCatalog)+" without knowing K using seed "+str(seedNumber)+"...\n")
        print("#####################\n")
        if os.path.exists(outputDir+"/SBS96/Suggested_Solution/Decomposed_Solution"):
            continue
        sig.sigProfilerExtractor("csv", outputDir, inputCatalog, refgen="GRCh37", 
        startProcess = exactK, endProcess = exactK, totalIterations = numiters, 
        cpu = numcpus, hierarchy = False, mtype = ["default"], exome = False)
		
## Restore old working directory		
os.chdir(oldWorkingDir)




