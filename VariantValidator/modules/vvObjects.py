import os
from configparser import ConfigParser,RawConfigParser
import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.dataproviders.seqfetcher
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.sequencevariant
import hgvs.validator
import hgvs.exceptions
import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.normalizer
import re
#import io
from vvDatabase import vvDatabase
from vvLogging import logger

# Custom Exceptions
class variantValidatorError(Exception):
    pass

'''
This file contains the validator object, which is instantiated in order to perform validator functions.
The validator contains configuration information and permanent copies of database links and the like.
Much of the validator's inner workings are stored in special one-off function container objects:
validator.db : The validator's MySQL database access functions

The validator configuration is stored in ~/.config/VariantValidator/config.ini . This is loaded
when the validator object is initialized.

Running variant validator should hopefully be as simple as writing a script like this:
import VariantValidator

val=Validator()
val.validate("some kind of gene situation","the transcripts to use")

'''

'''        
        Renaming of variables :
        'seqrepo_directory': HGVS_SEQREPO_DIR,           #self.seqrepoPath
        'uta_url': UTA_DB_URL,                           #self.utaPath
        'py_liftover_directory': PYLIFTOVER_DIR,         #self.liftoverPath
        'variantvalidator_data_url': VALIDATOR_DB_URL,   #self.db.path
        'entrez_id': ENTREZ_ID,                          #self.entrezID
        'variantvalidator_version': VERSION,             #self.version
        'variantvalidator_hgvs_version': hgvs_version,   #self.hgvsVersion
        'uta_schema': str(hdp.data_version()),           #self.uta_schema
        'seqrepo_db': HGVS_SEQREPO_DIR.split('/')[-1]    #self.seqrepoVersion
'''



class Validator():
    # This object contains configuration options.
    def __init__(self):
        # First load from the configuration file, if it exists.
        configName="config.ini"
        homePath=os.path.expanduser("~")
        configPath=os.path.join(homePath,".config","VariantValidator")
        if not os.path.isdir(configPath):
            os.makedirs(configPath)
        # Now configpath points to the config file itself.
        configPath=os.path.join(configPath,configName)
        # Does the file exist?
        if not os.path.exists(configPath):
            self.createConfig(configPath)

        # Load the configuration file.
        with open(configPath) as file:
            lines=file.read()
        config=RawConfigParser(allow_no_value=True)
        config.read(configPath)
        # The custom vvLogging module will set itself up using the VALDIATOR_DEBUG environment variable.
        logString = config["logging"]['string']
        os.environ["VALIDATOR_DEBUG"] = logString

        # Handle databases
        self.entrezID=config["EntrezID"]["entrezID"]
        if config["seqrepo"]["location"]!=None:
            self.seqrepoVersion=config["seqrepo"]["version"]
            self.seqrepoPath=config["seqrepo"]["location"]+self.seqrepoVersion
            os.environ['HGVS_SEQREPO_DIR']=self.seqrepoPath
        else:
            raise ValueError("The seqrepo location has not been set in ~/.config/VariantValidator/config.ini")
        os.environ['UTA_DB_URL']=config["uta"]["location"]+config["uta"]["version"]
        self.utaPath=config["uta"]["location"]+config["uta"]["version"]
        self.dbConfig = {
            'user':    config["mysql"]["user"],
            'password':config["mysql"]["password"],
            'host':    config["mysql"]["host"],
            'database':config["mysql"]["database"],
    	    'raise_on_warnings': True
        }
        self.db=vvDatabase(self,self.dbConfig)
        # Set up versions
        __version__ = config["variantValidator"]['version']
        self.version=__version__
        if re.match('^\d+\.\d+\.\d+$', __version__) is not None:
            self.releasedVersion=True
            _is_released_version = True
        else:
            self.releasedVersion=False
        self.hgvsVersion=hgvs.__version__

        # Set up other configuration variables
        self.liftoverPath=config["liftover"]["location"]
        if not self.liftoverPath==None:
            os.environ['PYLIFTOVER_DIR']=self.liftoverPath
        self.entrezID=config["EntrezID"]['entrezid']

        # Set up HGVS
        # Configure hgvs package global settings
        hgvs.global_config.uta.pool_max = 25
        hgvs.global_config.formatting.max_ref_length = 1000000
        # Create HGVS objects
        self.hdp = hgvs.dataproviders.uta.connect(pooling=True)
        self.hp = hgvs.parser.Parser() #P arser
        self.vr = hgvs.validator.Validator(self.hdp) # Validator
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp) # Variant mapper
        # Create a lose vm instance
        self.lose_vm = hgvs.variantmapper.VariantMapper(self.hdp,
                                                   replace_reference=True,
                                                   prevalidation_level=None
                                                   )
        self.nr_vm = hgvs.variantmapper.VariantMapper(self.hdp, replace_reference=False)
        self.sf = hgvs.dataproviders.seqfetcher.SeqFetcher() # Seqfetcher
        # Set standard genome builds
        self.genome_builds = ['GRCh37', 'hg19', 'GRCh38']
        self.uta_schema = str(self.hdp.data_version())

        #Transfer function handle from other file.
        self.validate=vvCore.validate

    def validate(self):
        pass

    def createConfig(self,outPath):
        # This function reads from the default configuration file stored in the same folder as this module,
        # and transfers it to outPath.
        # Outpath should include a filename.
        lines=[]
        inPath=os.path.join(os.path.dirname(os.path.realpath(__file__)),"defaultConfig.ini")
#        print(os.path.join(inPath,"defaultConfig.ini"))
        with open(inPath) as file:
            for l in file:
                lines.append(l)
        with open(outPath, "w") as file:
            for l in lines:
                file.write(l)
                

class Validation():
    #Validation objects contain a number of variant interpretations
    pass

class ValOutput():
    #This object contains a single possible interpretation of a variant
    pass


