"""
Exon-numbering

Authors: Katie Williams and Katherine Wingfield

This code will ultimately aim to provide exon numbering information for VariantValidator

"""

# Import the relevant packages/functions
import requests #this is needed to talk to the API
import json #this is needed to format the output data
import re 


# Function 1 
# Call the VV endpoint of the API for a BRCA1 transcript variant (originated from Clinvar)

input = "BRCA1"

# Function to request BRCA1 data from the gene2transcripts VariantValidator API
# This function will request info from https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/BRCA1 
base_url_VV = "https://rest.variantvalidator.org/"
server_G2T = "/VariantValidator/tools/gene2transcripts/"
gene_query = "BRCA1"


parameters = '?content-type=application/json'

def request_sequence(base_url, server, gene_name, parameters):
    url = base_url + server + gene_name + parameters
    print(url)
    #make the request and pass the object to the function
    response = requests.get(url) #this is the code that actually queries the API
    print("Querying " + url)
    return response

#request the sequence.
response = request_sequence(base_url_VV, server_G2T, gene_query, parameters) #added transcript_id in place of 'ENST00000297261'

#Convert response (JSON) to python dictionary
response_dictionary = response.json()

rd_type = type(response_dictionary) #troubleshooting checks - is this outputting a dictionary?
print(rd_type)

#Print the response
print(json.dumps(response_dictionary, sort_keys=True, indent=4, separators=(',', ': ')))