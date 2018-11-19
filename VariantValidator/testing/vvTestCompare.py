#Saving script

import vvTestFunctions as fn

masterDirectory="testOutputsMaster101"
testDirectories=["testOutputs246","testOutputsITS"]

for d in testDirectories:
    print("Comparing "+masterDirectory+" and "+d)
    fn.compareBatches(masterDirectory,d)

