#Saving script

import vvTestFunctions as fn
#from VariantValidator import Validator
import VariantValidator as vv
import os

val=vv.Validator()
os.environ["ADD_LOGS"]="True"
fn.generateTestFolder("testOutputs","inputVariants.txt",val)
