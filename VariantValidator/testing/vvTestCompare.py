#Saving script

import vvTestFunctions as fn

masterDirectory="testOutputsMasterITS"
testDirectories=["testOutputs","testOutputs246"]

for d in testDirectories:
    print("Comparing "+masterDirectory+" and "+d)
    fn.compareBatches(masterDirectory,d)

