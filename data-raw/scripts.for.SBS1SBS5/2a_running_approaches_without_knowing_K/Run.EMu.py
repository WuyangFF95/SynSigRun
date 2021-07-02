
# To convert ICAMS-formatted catalog to 
# catalog format acceptable by EMu,
# Prep.For.Running.EMu.R needs to be run
# before running this script.

#################################################################################################
###### load prerequisites
#################################################################################################

import os,sys,subprocess

#### Read old working directory
oldWorkingDir = os.getcwd()

# Set working directory to the folder which contains results of
# computational approaches on SBS1-SBS5-correlated data sets
# before running this script.
#
# PATH = "<path_to_results_on_SBS1-SBS5-correlated_datasets>"
#
# os.setcwd(PATH)
topLevelFolder4Run = "./2a.Full_output_K_unspecified"

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




for nrun in range(1,21):
    for datasetName in datasetNames:
        inputPath = topLevelFolder4Run+"/EMu.results/"+datasetName+"/run."+str(nrun)+"/"
        inputCatalog = "/".join([inputPath,"ground.truth.syn.catalog.tsv"])
        outputPath = inputPath
        print("\n\n======================================\n")
        print(str(nrun)+"-th running EMu for dataset "+str(datasetName)+" ...")
        print("\n\n======================================\n")
        ## The first argument should be replaced by the locatation of the compiled EMu in your machine.
        arguments = ['/home/wuyang/practice/3_Signature_Challenge/EMu/EMu/build/EMu',
            '--mut',inputCatalog,
            '--opp','human-genome',
            '--pre',outputPath]
        process = subprocess.Popen(" ".join(arguments),shell = True)
        process.wait()


## Restore old working directory		
os.chdir(oldWorkingDir)
