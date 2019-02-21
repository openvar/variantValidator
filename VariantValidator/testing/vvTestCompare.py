#Saving script

import vvTestFunctions as fn

masterDirectory="testOutputsMasterITS	"
testDirectories=["testOutputs"]

for d in testDirectories:
    print(("Comparing "+masterDirectory+" and "+d))
    fn.compareBatches(masterDirectory,d)

