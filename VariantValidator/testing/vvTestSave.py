#Saving script

import vvTestFunctions as fn
import sys
from StringIO import StringIO
import sqlite3
import os

class vvHub():
    #Variant validator configuration hub object
    def __init__(self):
        seqrepo_current_version='2018-08-21'
        HGVS_SEQREPO_DIR='/home/buran/documents/workspace/ITS/seqrepo/'+seqrepo_current_version
        os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
        self.hvgsSeqrepoPath=HGVS_SEQREPO_DIR
        uta_current_version='uta_20180821'
        UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
        os.environ['UTA_DB_URL']=UTA_DB_URL
        self.utaPath=UTA_DB_URL
        import VariantValidator.variantanalyser.vvLogging as vvLogging
        self.logger=vvLogging.logger
        from VariantValidator import variantValidator as vv
        self.vv=vv
        self.vv.my_config()


hub=vvHub()
fn.generateTestFolder("testOutputs","inputVariants.txt",hub)