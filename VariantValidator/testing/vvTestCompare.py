#Saving script

import vvTestFunctions as fn

masterDirectory="outputs101"
testDirectory1="outputs246"
testDirectory2="outputsITS"
testDirectory2="outputs249"

print("Comparing "+masterDirectory+" and "+testDirectory1)
fn.compareBatches(masterDirectory,testDirectory1)
print("Comparing "+masterDirectory+" and "+testDirectory2)
fn.compareBatches(masterDirectory,testDirectory2)
print("Comparing "+masterDirectory+" and "+testDirectory3)
fn.compareBatches(masterDirectory,testDirectory2)

