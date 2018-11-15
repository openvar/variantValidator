#PJDP testing suite for variant validator

#Run this test to validate all variants and SAVE the results for comparison with a different version.
#The input variants file should contain a bunch of variants on each line in quotes. Anything outside the
#quotes is discarded.
inputVariants="inputVariants.txt"

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
    with open(inputVariants) as f:
        for l in f.readlines():
            print(l)



#variant='NG_005905.2:g.172252G>A'
#select_transcripts='all'
#selected_assembly='GRCh37'
#validation=vv.validator(variant,selected_assembly,select_transcripts)
#print json.dumps(validation, sort_keys=True, indent=4, separators=(',',":"))
