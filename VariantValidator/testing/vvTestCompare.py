#Saving script

import vvTestFunctions as fn

masterDirectory="outputs101"
testDirectory1="outputs246"
testDirectory2="outputsITS"

fn.compareBatches(masterDirectory,testDirectory1)

fn.compareBatches(masterDirectory,testDirectory2)

