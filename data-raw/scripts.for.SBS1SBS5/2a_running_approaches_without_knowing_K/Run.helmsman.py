
# To convert ICAMS-formatted catalog to 
# catalog format acceptable by helmsman,
# Prep.For.Running.helmsman.R needs to be run
# before running this script.

#################################################################################################
###### load prerequisites
#################################################################################################
import os,sys,subprocess



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


for seedNumber in seedNumbers:
    for datasetName in datasetNames:
        inputPath = "/".join([datasetName,"sp.sp/ExtrAttr/helmsman.NMF.results"])
        inputCatalog = "/".join([inputPath,"ground.truth.syn.catalog.tsv"])
        outputPath = "/".join([inputPath,"seed."+str(seedNumber)])
        arguments = ['python3','~/practice/3_Signature_Challenge/helmsman/helmsman.py',
            '--seed',str(seedNumber),'--mode','agg',
            '--input',inputCatalog,
            '--projectdir',outputPath,
            '--verbose',
            '--decomp','nmf','--rank','0']
        process = subprocess.Popen(" ".join(arguments),shell = True)
        process.wait()

## Restore old working directory		
os.chdir(oldWorkingDir)