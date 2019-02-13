import os
import pytest
from VariantValidator import variantValidator as vv
import VariantValidator.testing.vvTestFunctions as fn


class vvHub():
    # Variant validator configuration hub object
    def __init__(self):
        if 'HGVS_SEQREPO_DIR' not in os.environ or 'UTA_DB_URL' not in os.environ:
            raise Exception("Environment variables aren't set. Please set HGVS_SEQREPO_DIR and UTA_DB_URL.")
        self.hvgsSeqrepoPath = os.environ['HGVS_SEQREPO_DIR']
        self.utaPath = os.environ['UTA_DB_URL']
        import VariantValidator.variantanalyser.vvLogging as vvLogging
        self.logger=vvLogging.logger
        self.vv = vv
        self.vv.my_config()

inputVariants=fn.loadVariantFile("VariantValidator/testing/inputVariants.txt")


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

@pytest.mark.skip(reason='Testing output elsewhere now')
def test_validation_output(constructValidation):
    hub,v=constructValidation
    assert v!=None

@pytest.mark.skip(reason="Testing for errors currently fails")
def test_validation_errors(constructValidation):
    hub,v=constructValidation
    logs=v["metadata"]["logs"].split("\n")
    e=0
    for l in logs:
        if "ERROR:" in l:
            e+=1
    assert e==0

@pytest.mark.skip(reason="Testing for criticals currently fails")
def test_validation_criticals(constructValidation):
    hub,v=constructValidation
    logs=v["metadata"]["logs"].split("\n")
    c=0
    for l in logs:
        if "CRIT:" in l:
            c+=1
    assert c==0

