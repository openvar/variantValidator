"""
Exon_numbering

Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)

This code will ultimately aim to provide exon numbering information for VariantValidator

"""

# Import the relevant packages/functions
import requests #this is needed to talk to the API
import json #this is needed to format the output data
import re 


####### Code to request BRCA1 data from the gene2transcripts VariantValidator API #########
# This code will request info from https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/BRCA1 

# Define the URL information as strings
base_url_VV = "https://rest.variantvalidator.org/"
server_G2T = "/VariantValidator/tools/gene2transcripts/"
gene_query = "BRCA1"

#use ensembl rest api
server_variant = 'variant_recoder/human/'
base_url_ensembl = 'https://rest.ensembl.org/'

# Define the parameter for retrieving in a JSON format
parameters = '?content-type=application/json'

# Create a function that will call the API and retrieve the information
def request_sequence(base_url, server, gene_name, parameters):
    url = base_url + server + gene_name + parameters
    print(url)
    #make the request and pass the object to the function
    response = requests.get(url) #this is the code that actually queries the API
    print("Querying " + url)
    return response

# Request the information (about BRCA1)
response = request_sequence(base_url_VV, server_G2T, gene_query, parameters)

#Convert response (JSON) to python dictionary
response_dictionary = response.json()

rd_type = type(response_dictionary) #troubleshooting checks - is this outputting a dictionary?
#print(rd_type)

#Print the response
#print(json.dumps(response_dictionary, sort_keys=True, indent=4, separators=(',', ': ')))


############ Use a variant ID, and call VV API ##############

# Pre-determine the variant
variant_id = "NM_007294.4:c.1067A>G"

# Maybe write code to split this up?
#query the API with the variant
variant_response = request_sequence(base_url_ensembl, server_variant, variant_id, parameters)

#convert the response (JSON) to a python dictionary
variant_response_dictionary = variant_response.json()

#print response
#print(json.dumps(variant_response_dictionary, sort_keys=True, indent=4, separators=(',', ': ')))

#find genome co-ordinates for the variant:
genomic_coordinates = [variant_response_dictionary[0]['spdi'][0].split(':')[1]]
print(genomic_coordinates)

#check the genome coordinates are what you expect, genome_coordinates is a list so need to use [0] to get contents of list
if genomic_coordinates[0] == '43094463':
    print("yes this works")

print(response_dictionary["genomic_spans"])
    # print("yes this works (again)")

def lookup(key, dictionary=response_dictionary):
    if key in dictionary: return dictionary[key]
    for value in dictionary.values():
        if isinstance(value, dict):
            a = lookup(key, value)
            if a is not None:
                return a
    return None

print(lookup('NC_000017.11'))

#create dictionary of all intron exon positions for transcripts
exon_start_end_positions = {}
for transcript_id in variant_response_dictionary:
    exon_start_end_positions[transcript_id] = {"start_exon": start_exon, "end_exon": end_exon}

'''
define function to parse through exon_start_end_positions and return the intron
and exon numbering for the start and end of the query variant
'''
def finds_variant_in_dict(reference):
    for transcript_id in exon_start_end_positions:
        if transcript_id == reference:
            print("The exon numbering for " + transcript_id + " starts in " + start_exon + " ends in " + end_exon)
