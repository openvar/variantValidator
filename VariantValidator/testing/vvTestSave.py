#Saving script

import vvTestFunctions as fn
import sys
from StringIO import StringIO
import sqlite3
import os

try:
    print("Configuring for lamp")
    seqrepo_current_version='2018-08-21'
    HGVS_SEQREPO_DIR='/local/seqrepo/'+seqrepo_current_version
    os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
    uta_current_version='uta_20180821'
    UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
    os.environ['UTA_DB_URL']=UTA_DB_URL
    from VariantValidator import variantValidator as vv
    vv.my_config()
except sqlite3.OperationalError:
    print("Configuring for VM")
    seqrepo_current_version='2018-08-21'
    HGVS_SEQREPO_DIR='/home/pjdp2/seqrepo/'+seqrepo_current_version
    os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
    uta_current_version='uta_20180821'
    UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
    os.environ['UTA_DB_URL']=UTA_DB_URL
    try:
        from VariantValidator import variantValidator as vv
        vv.my_config()
    except sqlite3.OperationalError:
        print("Configuring for VM")
        seqrepo_current_version = '2018-08-21'
        HGVS_SEQREPO_DIR = '/Users/pjf9/variant_validator_data/seqrepo/' + seqrepo_current_version
        os.environ['HGVS_SEQREPO_DIR'] = HGVS_SEQREPO_DIR
        uta_current_version = 'uta_20180821'
        UTA_DB_URL = 'postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
        os.environ['UTA_DB_URL'] = UTA_DB_URL
        os.environ['PYLIFTOVER_DIR'] = '/Users/pjf9/variant_validator_data/pyLiftover/'
        from VariantValidator import variantValidator as vv
        vv.my_config()
except OSError:
    print("Configuring for personal linux")
    seqrepo_current_version='2018-08-21'
    HGVS_SEQREPO_DIR='/home/buran/documents/workspace/ITS/seqrepo/'+seqrepo_current_version
    os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
    uta_current_version='uta_20180821'
    UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
    os.environ['UTA_DB_URL']=UTA_DB_URL
    from VariantValidator import variantValidator as vv
    vv.my_config()


sysOut=StringIO()

#sys.stdout=sysOut

inputVariants="inputVariants.txt"
#saveOut="testJSON.json"

#fn.generateTestJSON(saveOut,inputVariants,sysOut)
fn.generateTestFolder("testOutputs",inputVariants)
