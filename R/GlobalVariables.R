# GlobalVariables.R

# Define variables used by SignatureAnalyzer-related functions
# in an environment called "envSA"
# The parent of this environment is "parent.frame()", in the
# package, it will point to "<environment: namespace:SynSigRun>"
envSA <- new.env()
## Define variables in envSA.
for(varName in c("INPUT","TEMPORARY"))
  assign(varName,NULL,envir = envSA)


# Dump variables used by tests into environment envTest.
envTest <- new.env()

