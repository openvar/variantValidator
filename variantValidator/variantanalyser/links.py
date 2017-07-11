# -*- coding: utf-8 -*-

# Module containing functions that use hgvs to return sequence files

# IMPORT REQUIRED PYTHON MODULES
import re

# Biotools modules
import biotools.translate

def reference(tx_ac_fasta_title, ac_seq):
	# Open a list to store the fasta file
	tx_ac_fasta = []
	# Append the title line
	tx_ac_fasta.append(tx_ac_fasta_title)
	# Remove any white space from the string (tabs or spaces)
	# ac_seq = textwrap.dedent(ac_seq).strip()
	# Apply a width of 60 characters to the string output
	# ac_seq = textwrap.fill(ac_seq, width=60)
	# Append the sequence
	# tx_ac_fasta.append(seq_text)
	# The sequence string should be 60 bp chunks
	# Split the strings into lists of 60 base chunks
	#seq_list = nsplit(ac_seq, n=60)
	# loop and append as required
	#for lines in seq_list:
		#tx_ac_fasta.append(lines + '\n')
	tx_ac_fasta.append(ac_seq)
	return tx_ac_fasta
	
# Coding Inversions - This comes from VariantValidator, not validation!!!!
def coding_inversion(ref_seq, del_seq, inv_seq, interval_start, interval_end):	
	
	# Open a list to store the fasta file
	sequence = ''
	
	# Use string indexing to check whether the sequences are the same
	test = ref_seq[interval_start-1:interval_end]
	# return test
	if test == del_seq:
		sequence = ref_seq[0:interval_start-1] + inv_seq + ref_seq[interval_end:]
		return sequence
	else:
		sequence = 'error'
		return sequence

# n. Inversions - This comes from VariantValidator, not validation!!!!
def n_inversion(ref_seq, del_seq, inv_seq, interval_start, interval_end):	
	
	# Open a list to store the fasta file
	sequence = ''
	
	# Use string indexing to check whether the sequences are the same
	test = ref_seq[interval_start-1:interval_end]
	# return test
	if test == del_seq:
		sequence = ref_seq[0:interval_start-1] + inv_seq + ref_seq[interval_end:]
		return sequence
	else:
		sequence = 'error'
		return sequence


# Protein variation information
def pro_inv_info(prot_ref_seq, prot_var_seq):
	info = {
			'variant' : 'true',
			'prot_del_seq' : '',
			'prot_ins_seq' : '',
			'edit_start' : 0,
			'edit_end' : 0,
			'terminate' : 'false',
			'ter_pos' : 0,
			'error' : 'false'
			}
	
	# Is there actually any variation?
	if prot_ref_seq == prot_var_seq:
		info['variant'] = 'false'
	else:
		# Deal with terminations
		term = re.compile("\*")
		if term.search(prot_var_seq):
			# Set the termination reporter to true
			info['terminate'] = 'true'
			# The termination position will be equal to the length of the variant sequence because it's a TERMINATOR!!!
			info['ter_pos'] = len(prot_var_seq)
			# cut the ref sequence to == size
			prot_ref_seq = prot_ref_seq[0:info['ter_pos']]
			prot_var_seq = prot_var_seq[0:info['ter_pos']] 
		
			# Whether terminated or not, the sequences should now be the same length
			if len(prot_var_seq) != len(prot_ref_seq):
				info['error'] = 'true'
				return info
			else:
				# Set the counter
				aa_counter = 0

				# Make list copies of the sequences to gather the required info
				ref = list(prot_ref_seq)
				var = list(prot_var_seq)
				
				# Loop through ref list to find the first missmatch position
				for aa in ref:
					if ref[aa_counter] == var[aa_counter]:
						aa_counter = aa_counter + 1
					else:
						break
				
				# Enter the start position
				info['edit_start'] = aa_counter + 1
				# Remove those elements form the list
				del ref[0:aa_counter]
				del var[0:aa_counter]

				# the sequences should now be the same length
				if len(ref) != len(var):
					info['error'] = 'true'
					return info
				else:
					# Reset the aa_counter but to go backwards
					aa_counter = 0
					# reverse the lists
					ref = ref[::-1]
					var = var[::-1]
					# Reverse loop through ref list to find the first missmatch position
					for aa in ref:
						if var[aa_counter] == '\*':
							break
						if aa == var[aa_counter]:
							aa_counter = aa_counter + 1
						else:
							break
					# Remove those elements form the list
					del ref[0:aa_counter]
					del var[0:aa_counter]
					# re-reverse the lists	
					ref = ref[::-1]
					var = var[::-1]
					
					# the sequences should now be the same length
					if len(ref) != len(var):
						info['error'] = 'true'
						return info
					else:
						# Enter the sequences
						info['prot_del_seq'] = ''.join(ref)
						info['prot_ins_seq'] = ''.join(var)
						info['edit_end'] = info['edit_start'] + len(ref) -1 
						return info

	
# Translate sequences passed through
def translate(ed_seq, cds_start):
	# ed_seq = ed_seq.replace('\n', '')
	ed_seq = ed_seq.strip()
	# Ensure the starting codon is in the correct position
	met = ed_seq[cds_start:cds_start+3]
	if (met == 'ATG') or (met == 'atg'):
		# Remove the 5 prime UTR
		sequence = ed_seq[cds_start:]
		# Translate
		trans = biotools.translate.translate(sequence)
		aain = list(trans)
		aaout = []
		count = 0
		while aain:
			if aain[count] != '*':
				aaout.append(aain[count])
				count = count +1
			else:
				aaout.append(aain[count])
				break
		translation = ''.join(aaout)
		# Apply a width of 60 characters to the string output
		# translation = textwrap.fill(translation, width=60)
		return translation
	else:
		translation = 'error'
		return translation  
	
	
# Function to convert tri coded amino acids to unicoded
def one_to_three(seq):
	aacode = {
   'A' : 'Ala', 'C' : 'Cys', 'D' :'Asp', 'E' :'Glu',
    'F' :'Phe', 'G' :'Gly', 'H' :'His' , 'I' :'Ile',
    'K' :'Lys' , 'L' :'Leu', 'M' :'Met', 'N' :'Asn' ,
    'P' :'Pro', 'Q' :'Gln' , 'R' :'Arg' , 'S' :'Ser',
    'T' :'Thr', 'V' :'Val' , 'W' :'Trp', 'Y' :'Tyr',
    '*': '*'}
	
	oned = list(seq)
	out = []
	for aa in oned:
		get_value = aacode.get(aa)
		out.append(get_value)
		
	threed_up = ''.join(out)
	return threed_up
	
	
	
	
	
	
	
	
	
	
	
	
		
