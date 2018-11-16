#PJDP testing suite for variant validator

#Run this test to validate all variants and SAVE the results for comparison with a different version.
#The input variants file should contain a bunch of variants on each line in quotes. Anything outside the
#quotes is discarded.

import unittest
import os
import pickle

import sqlite3

try:
    seqrepo_current_version='2018-08-21'
    #HGVS_SEQREPO_DIR='/home/pjdp2/seqrepo/'+seqrepo_current_version
    HGVS_SEQREPO_DIR='/local/seqrepo/'+seqrepo_current_version
    os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
    uta_current_version='uta_20180821'
    UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
    os.environ['UTA_DB_URL']=UTA_DB_URL
    #from VariantValidator import variantValidator as vv
    print "Databases loaded"
    from VariantValidator import variantValidator as vv
    vv.my_config()
    print("Configured for LAMP")
except sqlite3.OperationalError:
    seqrepo_current_version='2018-08-21'
    HGVS_SEQREPO_DIR='/home/pjdp2/seqrepo/'+seqrepo_current_version
    #HGVS_SEQREPO_DIR='/local/seqrepo/'+seqrepo_current_version
    os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
    uta_current_version='uta_20180821'
    UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
    os.environ['UTA_DB_URL']=UTA_DB_URL
    #from VariantValidator import variantValidator as vv
    print "Databases loaded"
    from VariantValidator import variantValidator as vv
    vv.my_config()
    print("Configured for VM")

def saveValidations(path,inputVariants):
    #Saves the results of running inputVariants to a folder given in saveDirectory.
    if not os.path.isdir(path):
        os.mkdir(path)
    variantArray=loadVariantList(inputVariants)
    #Go through the variant array, validating, and save the results.
    batch=validateBatch(variantArray)
    #Save copy of the resulting dictionary
    for i,v in enumerate(batch):
        with open(os.path.join(path,"variant"+str(i)+".txt") ,"w") as f:
            pickle.dump(v,f)

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

def loadValidations(path):
    #Saves the results of running inputVariants to a folder given in saveDirectory.
    out=[]
    for paths,dirs,files in os.walk(path):
        for filePath in files:
            with open(os.path.join(paths,filePath)) as f:
                out.append(pickle.load(f))
                #print(type(out[-1]))
    return out

def validateBatch(variantArray):
    #Returns an array of validations (themselves dictionary objects).
    out=[]
    selectTranscripts='all'
    selectedAssembly='GRCh37'
    for i,v in enumerate(variantArray):
        print("VALIDATING Variant"+str(i)+" "+str(i+1)+"/"+str(len(variantArray))+" "+str(v))
        try:
            out.append(vv.validator(v,selectedAssembly,selectTranscripts))
        except Exception as e:
            print("FATAL error processing variant: "+str(e))
            out.append({"ERROR":str(e)})
    return out

def compareValidations(v1,v2):
    #print(v1,v2)
    for vk in v1.keys():
        if not (vk in v2.keys()):
            print(vk,"not found in second variant")
            return False
    for vk in v2.keys():
        if not (vk in v1.keys()):
            print(vk,"not found in first variant")
            return False
    for vk in v1.keys():
        if not (v1[vk]==v2[vk]):
            print("Different variant values - "+str(vk)+":"+str(v1[vk])+"vs."+str(vk)+":"+str(v2[vk]))
            return False
    return True

def compareBatches(v1path,v2path):
    #Loads all files in validations folder and compares them
    outFlags=[]
    passScore=0
    v1batch=loadValidations(v1path)
    v2batch=loadValidations(v2path)
    for i,v in enumerate(v1batch):
        outFlags.append(compareValidations(v1batch[i],v2batch[i]))
        if outFlags[-1]:
            passScore+=1
    print("Passed "+str(passScore)+"/"+str(len(v1batch)))
    if passScore==len(v1batch):
        #Test passed.
        print("Validation sets are identical")
        return True
    else:
        print("Validation sets are not identical - differences are:")
        for i,v in enumerate(v1batch):
            if not outFlags[i]:
                print("Mismatch in validation "+str(i))
                print(v1batch[i])
                print("Verses")
                print(v2batch[i])
        return False


    

