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


Kexact = 2

for index in range(1,21):
    for datasetName in datasetNames:
        inputPath = "/".join([datasetName,"sp.sp/ExtrAttrExact/EMu.results/run."+str(index)+"/"])
        inputCatalog = "/".join([inputPath,"ground.truth.syn.catalog.tsv"])
        outputPath = inputPath
        print("\n\n======================================\n")
        print(str(index)+"-th running EMu for dataset "+str(datasetName)+" ...")
        print("\n\n======================================\n")
        arguments = ['~/practice/3_Signature_Challenge/EMu/EMu/build/EMu',
            '--force',str(Kexact),
            '--mut',inputCatalog,
            '--opp','human-genome',
            '--pre',outputPath]
        process = subprocess.Popen(" ".join(arguments),shell = True)
        process.wait()
