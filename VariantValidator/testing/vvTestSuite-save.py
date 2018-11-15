#PJDP testing suite for variant validator

#Run this test to validate all variants and SAVE the results for comparison with a different version.
#The input variants file should contain a bunch of variants on each line in quotes. Anything outside the
#quotes is discarded.
inputVariants="inputVariants.txt"
saveDirectory="referenceOutputs"

import unittest
import os
from VariantValidator import variantValidator as vv


def loadDatabases():
    seqrepo_current_version='2018-08-21'
    HGVS_SEQREPO_DIR='/home/pjdp2/seqrepo/'+seqrepo_current_version
    os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
    uta_current_version='uta_20180821'
    UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
    os.environ['UTA_DB_URL']=UTA_DB_URL

def config():
    vv.my_config()
    print "Configured"
    
def saveTestResults():
    #Saves the results of running inputVariants to a folder given in saveDirectory.
    if not os.path.isdir(saveDirectory):
        os.mkdir(saveDirectory)
    variantArray=loadVariantList(inputVariants)
    #Go through the variant array, validating, and save the results.
    batch=validateBatch(variantArray)
    #Save copy of the resulting dictionary
    for i,v in enumerate(batch):
        with open(os.path.join(saveDirectory,"variant"+str(i)+".txt") ,"w") as f:
            pickle.dump(out,f)

def loadVariantList(path):
    out=[]
    #Load up the input variant file, should be passed in path.txt. Extra space, commas and quotes will be stripped.
    with open(path) as f:
        for l in f.readlines():
            l=l.strip()
            if len(l)>3:
                if l[-1]==",":
                    l=l[:-1]
                if l[-1]=='"':
                    l=l[:-1]
                if l[0]=='"':
                    l=l[1:]
                out.append(l)
    return out

def validateBatch(variantArray):
    #Returns an array of validations (themselves dictionary objects).
    out=[]
    selectTranscripts='all'
    selectedAssembly='GRCh37'
    for i,v in enumerate(variantArray[:3]):
        print("VALIDATING",str(i)+"/"+str(len(variantArray)),v)
        out.append(vv.validator(v,selectedAssembly,selectTranscripts))
    return out

#Main chain
loadDatabases()
config()
saveTestResults()

#variant='NG_005905.2:g.172252G>A'
#validation=vv.validator(variant,selected_assembly,select_transcripts)
#print json.dumps(validation, sort_keys=True, indent=4, separators=(',',":"))
