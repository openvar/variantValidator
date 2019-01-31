#Saving script

import vvTestFunctions as fn
import sys
from StringIO import StringIO

from VariantValidator import variantValidator as vv

vv.my_config()

sysOut=StringIO()

#sys.stdout=sysOut

inputVariants="inputVariants.txt"
#saveOut="testJSON.json"

#fn.generateTestJSON(saveOut,inputVariants,sysOut)
fn.generateTestFolder("testOutputs",inputVariants)
