# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 17:13:08 2021

@author: naomi
"""

"""
Protein variant format:
    NP_000079.2:p.(Gly197Cys)
    NP_000079.2:p.(G197C)
"""

#Import modules
import re #needed to split the string with multiple delimiters
import json #needed to create json object

#Define data
#Will want to replace the variant_accession with a VV input in the long term
variant_accession = "NP_000079.2:p.(Ter97Gly)"
print(variant_accession)

#Split string to get amino acid information
#Note this code would also work to get just the nucleotide variant
variant_accession_split = re.split('[()]', variant_accession)
print(variant_accession_split)

#define the protein variant
protein_variant = variant_accession_split[1]
print(protein_variant)

#Use re to split the variant into numbers and letters
number_letter = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
protein_variant_split = number_letter.match(protein_variant).groups()
print(protein_variant_split)

#Use logic to determine variant type
#This currently works for three letter aa codes only, could be expanded to one letter
#Edit: added protein_SO_term variable to loop
if protein_variant_split[0] == protein_variant_split [2]:
    print("Variant is synonymous")
    protein_SO_term = "synonymous_variant"
elif protein_variant_split[0] != "Ter" and protein_variant_split[2] == "Ter":
    print("Variant is stop gain")
    protein_SO_term = "stop_gain"
elif protein_variant_split[0] == "Ter" and protein_variant_split[2] != "Ter":
    print("Variant is stop loss")
    protein_SO_term = "stop_lost"
elif protein_variant_split[1] == "1" and protein_variant_split[0] == "Met" \
    and protein_variant_split[2] != "Met":
        print("Variant is start lost")
        protein_SO_term = "start_lost"
elif protein_variant_split[0] != "Ter" and protein_variant_split[2] != "Ter" \
    and protein_variant_split[1] != "1" and protein_variant_split[0] != "Met" \
    and protein_variant_split[0] != protein_variant_split[2]:
        print("Variant is missense")
        protein_SO_term = "missense_variant"
else:
    print("Variant type not recognised")
    raise SystemExit(0)

#Open an empty dictionary to store the response
SO_terms_dict = {}

#Add accession + protein SO term to dictionary
SO_terms_dict['Accession'] = variant_accession
SO_terms_dict['SO term'] = protein_SO_term

#Convert the dictionary to a json
SO_terms_output = json.dumps(SO_terms_dict, sort_keys=True, indent=4, separators=(',', ': '))
print(SO_terms_output)
        