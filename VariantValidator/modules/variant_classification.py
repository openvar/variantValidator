# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 17:13:08 2021

@author: naomi
@author: Ali
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


#Object that at the moment simply creates a dictionary and has capacity to
# add further entries. At the moment this is a simple repository which we 
#can use in later editions to provide/populate further variant information.
 
class Ensemble_reference_dict:
#initiator construct that creates a dictionary attribute when an instance of 
#the class is made. 
    def __init__(self):
        self.term_definitions = {}
        
 #Takes values for each SO descriptor and places in a dictionary using 
 # the SO term as a key and displays subsequent information as a list value,
 # for that key.    
    def add_entry(self, term, description, SOnumber, display_term, impact):
        self.term_definitions[term] = [description, SOnumber,display_term, impact]

#Calling an instance of the Ensemble_reference_dict class and populating the
#dictionary.
Ensemble_reference = Ensemble_reference_dict()
Ensemble_reference.add_entry("frameshift_variant", "A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three", "SO:0001589", "Frameshift variant", "High")
Ensemble_reference.add_entry("stop_gained", "A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript", "SO:00015872", "Stop gained", "High")
Ensemble_reference.add_entry("stop_lost", "A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript", "SO:0001578", "Stop lost", "High")
Ensemble_reference.add_entry("start_lost", "A codon variant that changes at least one base of the canonical start codon","SO:0002012", "Start lost", "High")
Ensemble_reference.add_entry("transcript_amplification", "A feature amplification of a region containing a transcript", "SO:0001889", "Transcript amplification", "High")

#couple of test print statements to access all and specific entries.
print(Ensemble_reference.term_definitions)
print(Ensemble_reference.term_definitions['stop_gained'][2])
