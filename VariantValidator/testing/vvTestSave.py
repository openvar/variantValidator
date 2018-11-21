#Saving script

import vvTestFunctions as fn
import sys

sys.stdout=open("vvTestSaveOutput.txt","w")

inputVariants="inputVariants.txt"
saveDirectory="testOutputs"

fn.generateTestFolder(saveDirectory,inputVariants)

