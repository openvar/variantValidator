#Saving script

import vvTestFunctions as fn
from VariantValidator import Validator
import os

val=Validator()
os.environ["ADD_LOGS"]="True"
fn.generateTestFolder("testOutputsReworked","inputVariants.txt",val)
