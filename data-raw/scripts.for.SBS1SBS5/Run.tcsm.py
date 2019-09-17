import os,sys,subprocess

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


for seedInUse in seedsInUse:
    for datasetName in datasetNames:
        inputPath = "/".join([datasetName,"sp.sp/ExtrAttr/EMu.results/run."+str(index)+"/"])
        inputCatalog = "/".join([inputPath,"ground.truth.syn.catalog.tsv"])
        outputPath = inputPath
        arguments = ['~/practice/3_Signature_Challenge/EMu/EMu/build/EMu',
            '--mut',inputCatalog,
            '--opp','human-genome',
            '--pre',outputPath]
        process = subprocess.Popen(" ".join(arguments),shell = True)
        process.wait()
