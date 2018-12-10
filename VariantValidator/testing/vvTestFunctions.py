#PJDP testing suite for variant validator

#Run this test to validate all variants and SAVE the results for comparison with a different version.
#The input variants file should contain a bunch of variants on each line in quotes. Anything outside the
#quotes is discarded.

import os
import pickle
import json
import sys

import sqlite3

try:
    print("Configuring for lamp")
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
    print("Configuring for VM")
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
except OSError:
    print("Configuring for personal linux")
    seqrepo_current_version='2018-08-21'
    HGVS_SEQREPO_DIR='/home/buran/documents/workspace/ITS/seqrepo/'+seqrepo_current_version
    #HGVS_SEQREPO_DIR='/local/seqrepo/'+seqrepo_current_version
    os.environ['HGVS_SEQREPO_DIR']=HGVS_SEQREPO_DIR
    uta_current_version='uta_20180821'
    UTA_DB_URL='postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
    #export postgresql://uta_admin:uta_admin@127.0.0.1/uta/uta_20180821
    os.environ['UTA_DB_URL']=UTA_DB_URL
    #from VariantValidator import variantValidator as vv
    print "Databases loaded"
    from VariantValidator import variantValidator as vv
    vv.my_config()
    print("Configured for VM")


def generateTestFolder(path, inputVariants):
    #Saves the results of running inputVariants to a folder given in saveDirectory.
    if not os.path.isdir(path):
        os.mkdir(path)
    variantArray=loadVariantFile(inputVariants)
    #Go through the variant array, validating, and save the results.
    batch=validateBatch(variantArray)
    #Save copy of the resulting dictionary
    saveValidationsAsFolder(path,batch)

def generateTestJSON(path, inputVariants,sysOut):
    variantArray=loadVariantFile(inputVariants)
    #Go through the variant array, validating, and save the results.
    batch=validateBatch(variantArray)
    #batch.append(sysOut.getvalue())
    #Save copy of the resulting dictionary
    saveValidationsAsJSON(path,batch)

def saveValidationsAsFolder(path, validations):
    #Pickles validation dictionaries into the given folder.
    for i,v in enumerate(validations):
        with open(os.path.join(path,"variant"+str(i)+".txt") ,"w") as f:
            pickle.dump(v,f)

def saveValidationsAsJSON(path,validations):
    #Saves a set of validations (v is a list of dictionaries) or a bunch of validations (v is a list of dictionaries)
    #as the json given in path. The name of the file will be that of the input variant string.
    jOut=json.dumps(validations)
    with open(path,"w") as f:
        f.write(jOut)
    print("JSON saved to "+path)

def loadVariantFile(path):
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

def saveVariantFile(path, variants):
    #Saves a variant input array (a bunch of strings) into a new text file given by path.
    with open(path,"w") as f:
        for v in variants:
            f.write(v+"\n")

def mergeVariantList(variants1,variants2):
    #Merges two lists of variants, avoiding duplicants.
    out=[]
    for v in variants1:
        if not v in out:
            out.append(v)
    for v in variants2:
        if not v in out:
            out.append(v)
    return out

def loadValidations(path):
    #Loads a set of validations from the folder given in path.
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
        except KeyboardInterrupt:
            print("Exiting...")
            sys.exit()
        except Exception as e:
            print("FATAL error processing variant: "+str(e))
            out.append({"ERROR":str(e)})
    return out

def retrieveVariant(validation):
    #Returns the variant string (if possible) from a validation.
    out=None
    for v in validation.values():
        try:
            if type(v)==type({}) and "submitted_variant" in v.keys():
                out=v["submitted_variant"]
                return out
        except (KeyError, TypeError, AttributeError):
            pass
    raise AttributeError("Validation does not contain the original variant string")

def compareValidations(v1,v2,id):
    #print(v1,v2)
    for vk in v1.keys():
        if not (vk in v2.keys()):
#            print("tag "+vk+" : "+str(v1[vk])+" not found in second variant")
            print("Variant "+str(id)+": Tag "+vk+" not found in second variant")
            return False
    for vk in v2.keys():
        if not (vk in v1.keys()):
#            print("tag "+vk+" : "+str(v2[vk])+" not found in first variant")
            print("Variant "+str(id)+": Tag "+vk+" not found in first variant")
            return False
    for vk in v1.keys():
        if not (v1[vk]==v2[vk]):
            if type(v1[vk])==type(dict()) or type(v2[vk])==type(dict()):
                print("Variant " + str(id) + ": Different tag values for key " + str(vk))
            else:
                print("Variant "+str(id)+": Different tag values - "+str(vk)+" : "+str(v1[vk])+" vs. "+str(vk)+" : "+str(v2[vk]))
            return False
    return True

def compareBatches(v1path,v2path):
    #Loads all files in validations folder and compares them
    outFlags=[]
    passScore=0
    v1batch=loadValidations(v1path)
    v2batch=loadValidations(v2path)
    print("Comparing validation sets...")
    for i,v in enumerate(v1batch):
#        print("Comparing validation "+str(i))
        outFlags.append(compareValidations(v1batch[i],v2batch[i],i))
        if outFlags[-1]:
            passScore+=1
    if passScore==len(v1batch):
        #Test passed.
        print("Validation sets are identical, "+str(passScore)+" passed")
        return True
    else:
        print("Validation sets are NOT identical, passed " + str(passScore) + "/" + str(len(v1batch)))
        #for i,v in enumerate(v1batch):
            #if not outFlags[i]:
                #print("Mismatch in validation "+str(i))
                #print(v1batch[i])
                #print("Verses")
                #print(v2batch[i])
        return False


    

