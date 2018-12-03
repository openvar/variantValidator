#Saving script

import vvTestFunctions as fn
import sys
from StringIO import StringIO

sysOut=StringIO()

#sys.stdout=sysOut

inputVariants="inputVariants.txt"
saveOut="testJSON.json"

fn.generateTestJSON(saveOut,inputVariants,sysOut)

