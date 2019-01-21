#Saving script

import vvTestFunctions as fn
from vvObjects import Validator
import os

val=Validator()
os.environ["ADD_LOGS"]="True"
fn.generateTestFolder("testOutputsReworked","inputVariants.txt",val)
