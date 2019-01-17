#Why isn't this discovered I wonder.
import os
import pytest
import VariantValidator.testing.vvTestFunctions as fn

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

inputVariants=fn.loadVariantFile("VariantValidator/testing/inputVariants.txt")

'''
print("Configuring for personal linux")
seqrepo_current_version='2018-08-21'
HGVS_SEQREPO_DIR='/home/buran/documents/workspace/ITS/seqrepo/'+seqrepo_current_version
os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
uta_current_version='uta_20180821'
UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
os.environ['UTA_DB_URL']=UTA_DB_URL
from VariantValidator import variantValidator as vv
vv.my_config()
'''

def constructHub():
    hub=vvHub()
    return hub

@pytest.fixture(params=inputVariants[:])
def constructValidation(request):
    hub=constructHub()
#    print request.param
    selectTranscripts='all'
    selectedAssembly='GRCh37'
    os.environ["ADD_LOGS"]="True"
    return hub,hub.vv.validator(request.param,selectedAssembly,selectTranscripts)

def test_validation_output(constructValidation):
    hub,v=constructValidation
    assert v!=None

def test_validation_errors(constructValidation):
    hub,v=constructValidation
    logs=v["metadata"]["logs"].split("\n")
    e=0
    for l in logs:
        if "ERROR:" in l:
            e+=1
    assert e==0

def test_validation_criticals(constructValidation):
    hub,v=constructValidation
    logs=v["metadata"]["logs"].split("\n")
    c=0
    for l in logs:
        if "CRIT:" in l:
            c+=1
    assert c==0

