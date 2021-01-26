"""
Exon_numbering

Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)

This code will ultimately aim to provide exon numbering information for VariantValidator

"""

# Import the relevant packages/functions
import requests #this is needed to talk to the API
import json #this is needed to format the output data
import re #this is used to manipulate the variant nomenclature

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

    # make the request and pass the object to the function
    response = requests.get(url)  # this is the code that actually queries the API
    print("Querying " + url)
    return response


# function to find exon numbering for a given variant
def finds_exon_number(variant):
    
    # extract the transcript ID from the variant nomenclature
    transcript_id = variant.split(":")[0]

    # request variant data from the gene2transcripts VariantValidator API
    response = request_sequence(base_url_VV, server_G2T, transcript_id, parameters)
    
    # Convert the response (JSON) to a python dictionary
    response_dictionary = response.json()

    # Note, function 2 will pull back the exon structures for all the transcripts
    # so will need to filter out the transcript you are interested in based on the transcript ID.
    #this for loop finds the exon structure for the given transcript
    num_transcripts = len(response_dictionary["transcripts"])
    for i in range(num_transcripts):
        if response_dictionary["transcripts"][i]["reference"] == transcript_id:
            
            #returns an exon structure dictionary
            exon_structure_dict = response_dictionary["transcripts"][i]["genomic_spans"]

            #identifies the coding start:
            coding_start = response_dictionary['transcripts'][i]["coding_start"]
            break
    
    #find the variant position from the variant nomenclature
    coordinates = variant.split(":")[1].split(".")[1]
    coordinates = re.sub('[^0-9, +, \-, _]','', coordinates) #removes A,G,C,T from HGVS nomenclature
    
    #identify start and end of variant from input coordinates
    if '_' in coordinates:
        split = coordinates.split('_')
        start_position = split[0]
        end_position = split[1]
    else:
        start_position = coordinates
        end_position = coordinates

    # #take into account coding start and end
    # adjusted_start_position = coding_start + (start_position) -1 #640+1 split()
    # adjusted_end_position = coding_start + int(end_position) - 1

    #create empty output dictionary
    exon_start_end_positions = {}

    # Works out the exon/intron for the transcript variant for each aligned chromosomal or gene reference sequence
    # This dictionary will have the keys as the aligned chromosomal and gene reference sequences
    # And the values of these keys will be another dictionary
    # With keys, start_exon and end_exon
    # With respective values relating the the position of variant in the reference seqeuence
    # e.g. {NC_000: {"start_exon": "1", "end_exon": "1i"}, NC_0000 ...
    for transcript in exon_structure_dict.keys():
 
        for exon in exon_structure_dict[transcript]['exon_structure']:

            #runs to identify which exon the variant is in 
            #start position
            if '+' not in str(start_position) and '-' not in str(start_position):
                adj_start_position = int(start_position) + coding_start - 1
                if adj_start_position >= exon['transcript_start'] and adj_start_position <= exon['transcript_end']:
                    start_exon = str(exon['exon_number'])
            
            elif '+' in start_position:
                nearest_exon_boundary = start_position.split('+')[0]
                adj_nearest_exon_boundary = int(nearest_exon_boundary) + coding_start - 1
                if adj_nearest_exon_boundary == exon['transcript_end']:
                    start_exon = str(exon['exon_number']) + 'i'
            
            elif '-' in start_position:
                nearest_exon_boundary = start_position.split('-')[0]
                adj_nearest_exon_boundary = int(nearest_exon_boundary) + coding_start - 1
                if adj_nearest_exon_boundary == exon['transcript_start']:
                    start_exon = str(exon['exon_number'] - 1)+ 'i'
            #end position
            if  '+' not in str(end_position) and '-' not in str(end_position):
                adj_end_position = int(end_position) + coding_start - 1
                if adj_end_position >= exon['transcript_start'] and adj_end_position <= exon['transcript_end']:
                    end_exon = str(exon['exon_number'])
                 
            elif '+' in end_position:
                nearest_exon_boundary = end_position.split('+')[0]
                adj_nearest_exon_boundary = int(nearest_exon_boundary) + coding_start - 1
                if adj_nearest_exon_boundary == exon['transcript_end']:
                    end_exon =  str(exon['exon_number'])+ 'i'
            
            elif '-' in end_position:
                nearest_exon_boundary = start_position.split('-')[0]
                adj_nearest_exon_boundary = int(nearest_exon_boundary) + coding_start - 1
                if adj_nearest_exon_boundary == exon['transcript_start']:
                    end_exon = str(exon['exon_number'] - 1) + 'i'
            
        exon_start_end_positions[transcript] = {"start_exon": start_exon, "end_exon": end_exon}
    return exon_start_end_positions

# <LICENSE>
# Copyright (C) 2021 VariantValidator Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
