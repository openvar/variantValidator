#Saving script

import vvTestFunctions as fn
from vvObjects import Validator
import os

val=Validator()

fn.generateTestFolder("testOutputsReworked","inputVariants.txt",val)
