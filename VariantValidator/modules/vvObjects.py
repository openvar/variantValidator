import os
from configparser import ConfigParser,RawConfigParser
import io

#        uta_current_version='uta_20180821'
#        UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
#        seqrepo_current_version='2018-08-21'
#        HGVS_SEQREPO_DIR='/home/buran/documents/workspace/ITS/seqrepo/'+seqrepo_current_version

class Validator():
    #This object contains configuration options.
    def __init__(self,hgvsPath=None,utaPath=None):
        #First load from the configuration file, if it exists.
        configName="config.ini"
        homePath=os.path.expanduser("~")
        configPath=os.path.join(homePath,".config","VariantValidator")
        if not os.path.isdir(configPath):
            os.makedirs(configPath)
        #Now configpath points to the config file itself.
        configPath=os.path.join(configPath,configName)
        #Does the file exist?
        if not os.path.exists(configPath):
            self.createConfig(configPath)

        #Load the configuration file.
        with open(configPath) as file:
            lines=file.read()
        config=RawConfigParser(allow_no_value=True)
        #print(configPath)
        config.read(configPath)
        #print config.sections()
        print config["seqrepo"]["location"]
        '''
        #Load hgvs
        if hgvsPath!=None:
            os.environ['HGVS_SEQREPO_DIR']=hgvsPath
            self.hgvsPath=hgvsPath
        else:
            self.hgvsPath=hgvsPath
        if utaPath!=None:
            os.environ['UTA_DB_URL']=utaPath
            self.utaPath=utaPath
        else:
            self.utaPath=utaPath
            seqrepo_current_version='2018-08-21'
            HGVS_SEQREPO_DIR='/home/buran/documents/workspace/ITS/seqrepo/'+seqrepo_current_version
            #HGVS_SEQREPO_DIR='/local/seqrepo/'+seqrepo_current_version
            os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
            uta_current_version='uta_20180821'
            UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
            #export postgresql://uta_admin:uta_admin@127.0.0.1/uta/uta_20180821
            os.environ['UTA_DB_URL']=UTA_DB_URL
            #from VariantValidator import variantValidator as vv

        '''

    def validate(self):
        pass

    def createConfig(self,outPath):
        #This function reads from the default configuration file stored in the same folder as this module.
        #Outpath should include a filename.
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


