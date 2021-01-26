"""
Exon_numbering.py Module

Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)

See exon_numbering.md markdown for a full description on how this
module operates.
*** Please note this Markdown is not finished yet ***

This code will ultimately aim to provide exon numbering information for
VariantValidator

"""

import requests  # This is needed to talk to the API
import json      # This is needed to format the output data
import re        # This is used to manipulate the variant nomenclature

# Define all the URL information as strings
BASE_URL_VV = "https://rest.variantvalidator.org/"
SERVER_G2T = "VariantValidator/tools/gene2transcripts/"

# Define the parameter for retrieving in a JSON format
PARAMETERS = '?content-type=application/json'


# Create a function that will call an API and retrieve the information
def request_sequence(base_url, server, variant_or_transcript, parameters):
    url = base_url + server + variant_or_transcript + parameters

    # Query the API and pass the object to the function
    response = requests.get(url)
    print("Querying " + url)
    return response


# Function to find and output exon numbering for a given variant
def finds_exon_number(variant):

    # Extract the transcript ID from the variant nomenclature
    transcript_id = variant.split(":")[0]

    # Request variant data from the gene2transcripts VariantValidator API
    response = request_sequence(BASE_URL_VV, SERVER_G2T, transcript_id,
                                PARAMETERS)

    # Convert the response (JSON) to a python dictionary
    response_dictionary = response.json()

    # Filter out the response_disctionary for the variant transcript
    # This will find the exon structure dictionary for the given transcript
    # And select the coding start position number
    for i in range(len(response_dictionary["transcripts"])):
        if response_dictionary["transcripts"][i]["reference"] == transcript_id:

            # Returns an exon structure dictionary
            exon_structure_dict = response_dictionary["transcripts"][i]["genomic_spans"]

            # Returns the start of coding (will this need to correct the position)
            coding_start = response_dictionary['transcripts'][i]["coding_start"]
            break

    # Find the variant position from the variant nomenclature
    coordinates = variant.split(":")[1].split(".")[1]
    # Remove A,G,C,T and variant description type from HGVS nomenclature
    # leaves only numbers, +, -, and _
    coordinates = re.sub('[^0-9, +, \-, _]', '', coordinates)

    # Identify start and end of variant from input coordinates
    if '_' in coordinates:
        start_position, end_position = coordinates.split('_')
    else:
        # If SNV, then start = end position
        start_position = coordinates
        end_position = coordinates

    # Create empty output dictionary
    exon_start_end_positions = {}

    # Works out the exon/intron for the transcript variant for each aligned chromosomal or gene reference sequence
    # This dictionary will have the keys as the aligned chromosomal and gene reference sequences
    # And the values of these keys will be another dictionary
    # With keys, start_exon and end_exon
    # With respective values relating the the position of variant in the reference seqeuence
    # e.g. {NC_000: {"start_exon": "1", "end_exon": "1i"}, NC_0000 ...
    for transcript in exon_structure_dict:
 
        for exon in exon_structure_dict[transcript]['exon_structure']:

            # For loop that runs to identify which exon/inton the variant is in
            # 'i' denotes introns, i.e. exon 2i is intron 2
            # Separated by start and end position of the variant as they may be different
            # if the variant is not a SNP.

            # Start position
            if '+' not in str(start_position) and '-' not in str(start_position):
                # This works for positions in exons
                adj_start_position = int(start_position) + coding_start - 1
                if adj_start_position >= exon['transcript_start'] and adj_start_position <= exon['transcript_end']:
                    start_exon = str(exon['exon_number'])
            
            elif '+' in start_position:
                # This works for positions that are + the exon boundary
                nearest_exon_boundary = start_position.split('+')[0]
                adj_nearest_exon_boundary = int(nearest_exon_boundary) + coding_start - 1
                if adj_nearest_exon_boundary == exon['transcript_end']:
                    start_exon = str(exon['exon_number']) + 'i'
           
            elif '-' in start_position:
                # This works for positions that are - the exon boundary
                nearest_exon_boundary = start_position.split('-')[0]
                adj_nearest_exon_boundary = int(nearest_exon_boundary) + coding_start - 1
                if adj_nearest_exon_boundary == exon['transcript_start']:
                    start_exon = str(exon['exon_number'] - 1) + 'i'
            # End position
            if  '+' not in str(end_position) and '-' not in str(end_position):
                # This works for positions in exons
                adj_end_position = int(end_position) + coding_start - 1
                if adj_end_position >= exon['transcript_start'] and adj_end_position <= exon['transcript_end']:
                    end_exon = str(exon['exon_number'])
                
            elif '+' in end_position:
                # This works for positions that are + the exon boundary
                nearest_exon_boundary = end_position.split('+')[0]
                adj_nearest_exon_boundary = int(nearest_exon_boundary) + coding_start - 1
                if adj_nearest_exon_boundary == exon['transcript_end']:
                    end_exon = str(exon['exon_number']) + 'i'
            
            elif '-' in end_position:
                # This works for positions that are - the exon boundary
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
