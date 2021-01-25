"""
Exon_numbering

Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)

This code will ultimately aim to provide exon numbering information for VariantValidator

"""

# Import the relevant packages/functions
import requests #this is needed to talk to the API
import json #this is needed to format the output data
import re 


# Define all the URL information as strings
base_url_VV = "https://rest.variantvalidator.org/"
server_G2T = "VariantValidator/tools/gene2transcripts/"
gene_query = "BRCA1"
server_variant = 'VariantValidator/variantvalidator/'
genome_build = "GRCh38"

# Define the parameter for retrieving in a JSON format
parameters = '?content-type=application/json'

# Create a function that will call an API and retrieve the information
def request_sequence(base_url, server, gene_name, parameters):
    url = base_url + server + gene_name + parameters
    # print(url)
    # make the request and pass the object to the function
    response = requests.get(url) #this is the code that actually queries the API
    print("Querying " + url)
    return response


############ FUNCTION 1 ######################################################
# Use a variant ID, and call VV API

# Pre-determine the variant
variant_id = "NM_007294.3:c.1067A>G"

# Query the VV API with the variant
variant_response = request_sequence(base_url_VV, server_variant, genome_build + '/' + variant_id + '/all', parameters)

# Convert the response (JSON) to a python dictionary
variant_response_dictionary = variant_response.json()

# Print response
#print(json.dumps(variant_response_dictionary, sort_keys=True, indent=4, separators=(',', ': ')))


####### FUNCTION 2 #########################################################
# Code to request BRCA1 data from the gene2transcripts VariantValidator API 

# First, do using gene_query as "BRCA1"
# response = request_sequence(base_url_VV, server_G2T, gene_query, parameters)
# response_dictionary = response.json()

# Print the response
# print(json.dumps(response_dictionary, sort_keys=True, indent=4, separators=(',', ': ')))

# Now, rather than the gene symbol BRCA1, it pass a transcript ID
transcript_id = variant_id.split(":")[0]
response = request_sequence(base_url_VV, server_G2T, transcript_id, parameters)
response_dictionary = response.json()
# Print the response
#print(json.dumps(response_dictionary, sort_keys=True, indent=4, separators=(',', ': ')))

# Note, function 2 will pull back the exon structures for all the transcripts
# so will need to filter out the transcript you are interested in based on the transcript ID.

# Open an empty dictionary
dict1 = {}
# Add our variant
dict1['Variant HGVS'] = variant_id
# Format the ouput dictionary nicely and print
output_dict = json.dumps(dict1, sort_keys=True, indent=4, separators=(',', ': '))
#print the output dictionary
# print(output_dict)


#print(type(response_dictionary["transcripts"][0]))
#print(response_dictionary["transcripts"][0].keys())
# print(response_dictionary["transcripts"][3]["reference"])
# print(len(response_dictionary["transcripts"])) # there are 7 transcripts

num_transcripts = len(response_dictionary["transcripts"])
for i in range(num_transcripts):
    if response_dictionary["transcripts"][i]["reference"] == transcript_id:
        transcipt_accession = i
        brca_exon_structure = response_dictionary["transcripts"][i]["genomic_spans"]
        break

print(brca_exon_structure.keys())
#print(brca_exon_structure)

####### FUNCTION 3 #########################################################
# Works out the exon/intron for the transcript variant for each aligned chromosomal or gene reference sequence
# Set up output dictionary  
exon_start_end_positions = {}
# This dictionary will have the keys as the aligned chromosomal and gene reference sequences
# And the values of these keys will be another dictionary
# With keys, start_exon and end_exon
# With respective values relating the the position of variant in the reference seqeuence
# e.g. {NC_000: {"start_exon": "1", "end_exon": "1i"}, NC_0000 .... 

# #find genome co-ordinates for the variant:
#print(variant_response_dictionary.keys())
genomic_coordinates = variant_response_dictionary[variant_id]["primary_assembly_loci"][genome_build.lower()]['vcf']['pos']
print(type(genomic_coordinates)) #genomic coordinates = 43094464

#print(brca_exon_structure['NC_000017.11'])

def finds_exon_number(coordinates, exon_structure_dictionary=brca_exon_structure):
    for exon in brca_exon_structure['NC_000017.11']['exon_structure']:
        print(exon)
        if coordinates >= exon['genomic_start'] and coordinates <= exon['genomic_end']:
            print("Exon number is " + str(exon['exon_number']))
            return exon['exon_number']

print(finds_exon_number(int(genomic_coordinates)))

# def lookup(key, dictionary=brca_exon_structure):
#     if key in dictionary: return dictionary[key]
#     for value in dictionary.values():
#         if isinstance(value, dict):
#             a = lookup(key, value)
#             if a is not None:
#                 return a
#     return None

#print(lookup('NC_000017.11'))

#create dictionary of all intron exon positions for transcripts
# for transcript_id in variant_response_dictionary:
#     exon_start_end_positions[transcript_id] = {"start_exon": start_exon, "end_exon": end_exon}

'''
define function to parse through exon_start_end_positions and return the intron
and exon numbering for the start and end of the query variant
'''
# def finds_variant_in_dict(reference):
#     for transcript_id in exon_start_end_positions:
#         if transcript_id == reference:
#             print("The exon numbering for " + transcript_id + " starts in " + start_exon + " ends in " + end_exon)

#close