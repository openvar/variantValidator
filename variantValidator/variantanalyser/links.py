# -*- coding: utf-8 -*-
"""
links.py

is an extension module of functions.py. It was ofiginally built to provide necesary
functions for compiling reference sequence alignments
 
The module contains additional VariantValidator sub-functions. The majoirty of these functions require
hgvs Python package top-level functions or sub-functions contained in uta.py and 
seqfetcher.py 

"""

# IMPORT REQUIRED PYTHON MODULES
import re

# BioPython modules
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

"""
Function which predicts the protein effect of c. inversions
"""
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
			# Unless the termination codon has been disrupted
			if len(prot_var_seq) < len(prot_ref_seq):
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
				# Except if the termination codon was removed
				if len(ref) > len(var):
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
					
					# If the var is > ref, the ter has been removed, need to re-add ter to each
					if len(ref) < len(var):
						ref.append('*')
						if prot_var_seq[-1] == '*':
							var.append('*')
					# the sequences should now be the same length
					# Except if the ter was removed
					if len(ref) > len(var):
						info['error'] = 'true'
						return info				
					else:
						# Enter the sequences
						info['prot_del_seq'] = ''.join(ref)
						info['prot_ins_seq'] = ''.join(var)
						info['edit_end'] = info['edit_start'] + len(ref) -1 
						return info
						
"""
Translate c. reference sequences, including those that have been modified 
must have the CDS in the specified position
"""
def translate(ed_seq, cds_start):
	# ed_seq = ed_seq.replace('\n', '')
	ed_seq = ed_seq.strip()
	# Ensure the starting codon is in the correct position
	met = ed_seq[cds_start:cds_start+3]
	if (met == 'ATG') or (met == 'atg'):
		# Remove the 5 prime UTR
		sequence = ed_seq[cds_start:]
		coding_dna = Seq(str(sequence), IUPAC.unambiguous_dna)
		# Translate
		trans = coding_dna.translate()
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

	
	
"""
Convert single letter amino acid code to 3 letter code
"""
def one_to_three(seq):
	aacode = {
   'A' : 'Ala', 'C' : 'Cys', 'D' :'Asp', 'E' :'Glu',
    'F' :'Phe', 'G' :'Gly', 'H' :'His' , 'I' :'Ile',
    'K' :'Lys' , 'L' :'Leu', 'M' :'Met', 'N' :'Asn' ,
    'P' :'Pro', 'Q' :'Gln' , 'R' :'Arg' , 'S' :'Ser',
    'T' :'Thr', 'V' :'Val' , 'W' :'Trp', 'Y' :'Tyr',
    '*': 'Ter'}
	
	oned = list(seq)
	out = []
	for aa in oned:
		get_value = aacode.get(aa)
		out.append(get_value)
		
	threed_up = ''.join(out)
	
	return threed_up
						

"""
marked for removal
"""
# def reference(tx_ac_fasta_title, ac_seq):
# 	# Open a list to store the fasta file
# 	tx_ac_fasta = []
# 	# Append the title line
# 	tx_ac_fasta.append(tx_ac_fasta_title)
# 	# Remove any white space from the string (tabs or spaces)
# 	# ac_seq = textwrap.dedent(ac_seq).strip()
# 	# Apply a width of 60 characters to the string output
# 	# ac_seq = textwrap.fill(ac_seq, width=60)
# 	# Append the sequence
# 	# tx_ac_fasta.append(seq_text)
# 	# The sequence string should be 60 bp chunks
# 	# Split the strings into lists of 60 base chunks
# 	#seq_list = nsplit(ac_seq, n=60)
# 	# loop and append as required
# 	#for lines in seq_list:
# 		#tx_ac_fasta.append(lines + '\n')
# 	tx_ac_fasta.append(ac_seq)
# 	return tx_ac_fasta
	
""" 
Takes a reference sequence and inverts the specified position
"""
# n. Inversions - This comes from VariantValidator, not validation!!!!
def n_inversion(ref_seq, del_seq, inv_seq, interval_start, interval_end):	
	
	# Open a list to store the fasta file
	sequence = ''
	
	# Use string indexing to check whether the sequences are the same
	test = ref_seq[interval_start-1:interval_end]

	if test == del_seq:
		sequence = ref_seq[0:interval_start-1] + inv_seq + ref_seq[interval_end:]
		return sequence
	else:
		sequence = 'error'
		return sequence


# <LICENSE>

# </LICENSE>


"""
The following section contains old functions which may be recycled
Do not delete
"""

		
""" 
mark for removal
"""
# Coding Inversions - This comes from VariantValidator, not validation!!!!
#def coding_inversion(ref_seq, del_seq, inv_seq, interval_start, interval_end):	
# 	
# 	# Open a list to store the fasta file
# 	sequence = ''
# 	
# 	# Use string indexing to check whether the sequences are the same
# 	test = ref_seq[interval_start-1:interval_end]
# 	# return test
# 	if test == del_seq:
# 		sequence = ref_seq[0:interval_start-1] + inv_seq + ref_seq[interval_end:]
# 		return sequence
# 	else:
# 		sequence = 'error'
# 		return sequence
		

"""
Legacy function, May be recycled
"""
# SIMPLE DELETIONS
# def sim_del(tx_ac_fasta_title, ac_seq, interval_start, interval_end, variant):	
# 	
# 	# Open a list to store the fasta file
# 	tx_ac_fasta = []
# 	
# 	# Intronic positions need to be dealt with
# 	pl = re.compile('\+')
# 	mi = re.compile('\-')
# 	if (pl.search(variant)) or (mi.search(variant)):
# 		if interval_start == interval_end:
# 			# Append the title line
# 			tx_ac_fasta.append(tx_ac_fasta_title)
# 			# Append the sequence
# 			tx_ac_fasta.append(ac_seq)
# 			return tx_ac_fasta
# 		else:
# 			pass
# 
# 	# Append the title line
# 	tx_ac_fasta.append(tx_ac_fasta_title)
# 	
# 	# Remove any white space from the string (tabs or spaces)
# 	# ac_text = textwrap.dedent(ac_seq).strip()
# 	ac_text = ac_seq
# 				
# 	# Assign the original sequence to a list
# 	ac_list = list(ac_text)
# 				
# 	# Delete the required bases (elements)
# 	del ac_list[interval_start - 1 : interval_end]
# 				
# 	# Join the sequence back together
# 	seq_text = ''.join(ac_list)
# 					
# 	# Apply a width of 60 characters to the string output
# 	#seq_text = textwrap.fill(seq_text, width=60)
# 
# 	# Make edits DNA
# 	seq_text = seq_text.replace("U", "T")
# 	seq_text = seq_text.replace("u", "t")
# 	
# 	# Append the sequence
# 	tx_ac_fasta.append(seq_text)
# 		
# 	return tx_ac_fasta

"""
legacy function, May be recycled
"""
# SIMPLE INSERTIONS
# def sim_ins(tx_ac_fasta_title, ac_seq, interval_start, interval_end, edit, variant):	
# 	
# 	# Open a dictionary to store the fasta file
# 	tx_ac_fasta_dict = {
# 						'tx_ac_fasta' : [],
# 						'edit_len' : 0
# 						}
# 
# 	# Intronic positions need to be dealt with
# 	pl = re.compile('\+')
# 	mi = re.compile('\-')
# 	if (pl.search(variant)) or (mi.search(variant)):
# 		# Append the title line
# 		tx_ac_fasta_dict['tx_ac_fasta'].append(tx_ac_fasta_title)
# 		# Append the sequence
# 		tx_ac_fasta_dict['tx_ac_fasta'].append(ac_seq)
# 		return tx_ac_fasta_dict
# 	else:
# 		pass
# 	
# 	# Append the title line
# 	tx_ac_fasta_dict['tx_ac_fasta'].append(tx_ac_fasta_title)
# 	
# 	# Remove any white space from the string (tabs or spaces)
# 	# ac_text = textwrap.dedent(ac_seq).strip()
# 	ac_text = ac_seq
# 				
# 	# Extract inserted sequence and change into a list
# 	ins_seq = re.search(r"(([GATCUgatcu]+)$)", edit)
# 	ins_seq = ins_seq.group(1)							
# 	
# 	# Add the edit length to the dictionary
# 	tx_ac_fasta_dict['edit_len'] = len(ins_seq)
# 
# 	# Assign the original sequence to a list
# 	ac_list = list(ac_text)
# 				
# 	# Assign the bases up to the edit to the edit list
# 	edit_list = ac_list[0:interval_start]
# 				
# 	# Delete the completed bases from the sequence list
# 	del ac_list[0:interval_start]
# 	
# 	# Append the insertion and the remaining sequence
# 	for base in ins_list:
# 		edit_list.append(base)
# 
# 	for base in ac_list:
# 		edit_list.append(base)
# 	
# 	# Join the sequence back together
# 	seq_text = ''.join(edit_list)
# 	
# 	# Make edits DNA
# 	seq_text = seq_text.replace("U", "T")
# 	seq_text = seq_text.replace("u", "t")
# 
# 	# Apply a width of 60 characters to the string output
# 	# seq_text = textwrap.fill(seq_text, width=60)
# 		
# 	# Append the sequence
# 	tx_ac_fasta_dict['tx_ac_fasta'].append(seq_text)
# 	
# 	return tx_ac_fasta_dict


"""
legacy function, May be recycled
"""
# DELETION INSERTIONS
# def delins(tx_ac_fasta_title, ac_seq, interval_start, interval_end, edit, variant, type):
# 	# Open a dictionary to store the fasta file
# 	tx_ac_fasta_dict = {
# 						'tx_ac_fasta' : [],
# 						'edit_len' : 0
# 						}
# 
# 	# Intronic positions need to be dealt with
# 	pl = re.compile('\+')
# 	mi = re.compile('\-')
# 	if (pl.search(variant)) or (mi.search(variant)):
# 		if interval_start == interval_end:
# 			# Append the title line
# 			tx_ac_fasta_dict['tx_ac_fasta'].append(tx_ac_fasta_title)
# 			# Append the sequence
# 			tx_ac_fasta_dict['tx_ac_fasta'].append(ac_seq)
# 			# Append the insert length as zero
# 			tx_ac_fasta_dict['edit_len'] = 0
# 			return tx_ac_fasta_dict
# 		else:
# 			pass
# 	
# 	# Append the title line
# 	tx_ac_fasta_dict['tx_ac_fasta'].append(tx_ac_fasta_title)
# 	
# 	# Remove any white space from the string (tabs or spaces)
# 	# ac_text = textwrap.dedent(ac_seq).strip()
# 	ac_text = ac_seq
# 				
# 	# Extract inserted sequence and change into a list
# 	ins_seq = re.search(r"(([GATCUgatcu]+)$)", edit)
# 	ins_seq = ins_seq.group(1)
# 	# Calculate the edit length
# 	edit_len = len(ins_seq)
# 	# Balance out the insertions at exon_intron boundaries
# 	if type != ':g.':
# 		if mi.search(variant):
# 			if edit_len > (interval_end - interval_start +1):
# 				edit_len = (interval_end - interval_start +1)
# 				ins_seq = ins_seq[-edit_len:]
# 		else:
# 			pass
# 		if pl.search(variant):	
# 			if edit_len > (interval_end - interval_start +1):
# 				edit_len = (interval_end - interval_start +1)
# 				ins_seq = ins_seq[0:edit_len]
# 		else:
# 			pass
# 	else:
# 		pass				
# 	
# 	# Add the edit length to the dictionary
# 	tx_ac_fasta_dict['edit_len'] = edit_len
# 
# 	# List the insert
# 	ins_list = list(ins_seq)
# 	
# 	# Assign the original sequence to a list
# 	ac_list = list(ac_text)
# 				
# 	# Delete the required bases (elements) and add the insert list
# 	ac_list[interval_start - 1 : interval_end] = ins_list
# 				
# 	# Join the sequence back together
# 	seq_text = ''.join(ac_list)
# 	
# 	# Make edits DNA
# 	seq_text = seq_text.replace("U", "T")
# 	seq_text = seq_text.replace("u", "t")
# 	
# 	# Apply a width of 60 characters to the string output
# 	# seq_text = textwrap.fill(seq_text, width=60)
# 
# 	# Append the sequence
# 	tx_ac_fasta_dict['tx_ac_fasta'].append(seq_text)
# 	
# 	return tx_ac_fasta_dict
	
"""
legacy function, May be recycled
"""
# POINT MUTATIONS	
# def point(tx_ac_fasta_title, ac_seq, interval_start, edit, variant):
# 	# Open a dictionary to store the fasta file
# 	tx_ac_fasta_dict = {
# 						'tx_ac_fasta' : [],
# 						'flag' : ''
# 						}
# 	
# 	# Intronic positions need to be dealt with
# 	pl = re.compile('\+')
# 	mi = re.compile('\-')
# 	if (pl.search(variant)) or (mi.search(variant)):
# 		# Append the title line
# 		tx_ac_fasta_dict['tx_ac_fasta'].append(tx_ac_fasta_title)
# 		# Append the sequence
# 		tx_ac_fasta_dict['tx_ac_fasta'].append(ac_seq)
# 		return tx_ac_fasta_dict
# 	else:
# 		pass
# 
# 	# Append the title line
# 	tx_ac_fasta_dict['tx_ac_fasta'].append(tx_ac_fasta_title)
# 
# 	start_base = re.search(r"(^([GATCUgatcu]+))", edit)
# 	start_base = start_base.group(1)
# 			
# 	end_base = re.search(r"(([GATCUgatcu]+)$)", edit)
# 	end_base = end_base.group(1)
# 	
# 	if end_base == "U":
# 		end_base = end_base.replace("U", "T")
# 	if end_base == "u":
# 		end_base = end_base.replace("u", "t")
# 	if start_base == "U":
# 		start_base = start_base.replace("U", "T")
# 	if start_base == "u":
# 		start_base = start_base.replace("u", "t")
# 				
# 	# Remove any white space from the string (tabs or spaces)
# 	# ac_text = textwrap.dedent(ac_seq).strip()
# 	ac_text = ac_seq
# 				
# 	# Assign the original sequence to a list
# 	ac_list = list(ac_text)
# 				
# 	# Search the list at the correct location for the edit base
# 	if ac_list[interval_start -1] == start_base:			
# 		# Make the edit
# 		ac_list[interval_start -1] = end_base								
# 		# Join the sequence back together
# 		seq_text = ''.join(ac_list)			
# 		# Make edits DNA
# 		seq_text = seq_text.replace("U", "T")
# 		seq_text = seq_text.replace("u", "t")
# 		# Apply a width of 60 characters to the string output
# 		# seq_text = textwrap.fill(seq_text, width=60)
# 		# Append the sequence
# 		tx_ac_fasta_dict['tx_ac_fasta'].append(seq_text)
# 		return tx_ac_fasta_dict
# 	
# 	# Search the list at the correct location for the edit base
# 	if ac_list[interval_start -1] == start_base.lower():			
# 		# Make the edit
# 		ac_list[interval_start -1] = end_base.lower()								
# 		# Join the sequence back together
# 		seq_text = ''.join(ac_list)			
# 		# Make edits DNA
# 		seq_text = seq_text.replace("U", "T")
# 		seq_text = seq_text.replace("u", "t")
# 		# Apply a width of 60 characters to the string output
# 		#seq_text = textwrap.fill(seq_text, width=60)
# 		# Append the sequence
# 		tx_ac_fasta_dict['tx_ac_fasta'].append(seq_text)
# 		return tx_ac_fasta_dict
# 
# 	# Search the list at the correct location for the edit base
# 	else:			
# 		# Make the edit
# 		ac_list[interval_start -1] = end_base								
# 		# Join the sequence back together
# 		seq_text = ''.join(ac_list)			
# 		# Make edits DNA
# 		seq_text = seq_text.replace("U", "T")
# 		seq_text = seq_text.replace("u", "t")
# 		# Apply a width of 60 characters to the string output
# 		# seq_text = textwrap.fill(seq_text, width=60)
# 		# Append the sequence
# 		tx_ac_fasta_dict['tx_ac_fasta'].append(seq_text)
# 		if (pl.search(variant)) or (mi.search(variant)):
# 			tx_ac_fasta_dict['flag'] = ''
# 		else:	
# 			tx_ac_fasta_dict['flag'] = 'Warning: Variant sequence does not agree with reference sequence'
# 			#tx_ac_fasta_dict['flag'] = ac_seq
# 		return tx_ac_fasta_dict

"""
legacy function, May be recycled
"""
# DUPLICATIONS		
# def dupn(tx_ac_fasta_title, ac_seq, interval_start, interval_end, edit, variant):
# 	# Open a dictionary to store the fasta file
# 	tx_ac_fasta_dict = {
# 						'tx_ac_fasta' : [],
# 						'edit_len' : 0
# 						}
# 						
# 	# Intronic positions need to be dealt with
# 	pl = re.compile('\+')
# 	mi = re.compile('\-')
# 	if (pl.search(variant)) or (mi.search(variant)):
# 		# Append the title line
# 		tx_ac_fasta_dict['tx_ac_fasta'].append(tx_ac_fasta_title)
# 		# Append the sequence
# 		tx_ac_fasta_dict['tx_ac_fasta'].append(ac_seq)
# 		return tx_ac_fasta_dict
# 	else:
# 		pass
# 	
# 	# Append the title line
# 	tx_ac_fasta_dict['tx_ac_fasta'].append(tx_ac_fasta_title)			
# 						
# 	# Extract duplicated bases if sequence available
# 	bases = re.compile('[GATCUgatcu]+$')
# 	if bases.search(edit):
# 		dup_base = re.search(r"(([GATCUgatcu]+)$)", edit)
# 		dup_base = dup_base.group(1)
# 	else:
# 		# The duplicated sequence needs to be extracted based on coordinates,
# 		dup_base = ac_seq[interval_start -1:interval_end]
# 	
# 	# Append the duplicated sequence length
# 	tx_ac_fasta_dict['edit_len'] = len(dup_base)
# 	
# 	# Remove any white space from the string (tabs or spaces)
# 	# ac_text = textwrap.dedent(ac_seq).strip()
# 	ac_text = ac_seq
# 				
# 	# Assign the original sequence to a list
# 	ac_list = list(ac_text)
# 
# 	# Search the list at the correct location for the edit base
# 	if ac_list[interval_start -1] == dup_base[0]:
# 		#return 'Giraffes'			
# 		# Make the edit
# 		edit_list = ac_list[0:interval_start -1]
# 		# Remove the processed bases
# 		del ac_list[0:interval_start -1]
# 		# List the insert and append the bases
# 		ins_list = list(dup_base)
# 		for base in ins_list:
# 			edit_list.append(base)
# 		# Append the remaining sequence
# 		for base in ac_list:
# 			edit_list.append(base)
# 		
# 	# Join the sequence back together
# 	seq_text = ''.join(edit_list)
# 	#return seq_text
# 	
# 	# Make edits DNA
# 	seq_text = seq_text.replace("U", "T")
# 	seq_text = seq_text.replace("u", "t")
# 						
# 	# Apply a width of 60 characters to the string output
# 	# seq_text = textwrap.fill(seq_text, width=60)
# 
# 	# Append the sequence
# 	tx_ac_fasta_dict['tx_ac_fasta'].append(seq_text)
# 		
# 	return tx_ac_fasta_dict


"""
legacy function, May be recycled
"""
# CHUNKING
# Split a string into chunks of a given length and return a list
# def nsplit(s, n):
#     return [s[k:k+n] for k in xrange(0, len(s), n)]
	
	
"""
legacy function, May be recycled
"""
# PERFORM ALIGNING OF TWO KNOWN SEQUENCES
# Takes sequence strings and creates an alignment list
# def format_alignment(align1, align2, edit_type, interval_start, interval_end, begin, end, ins_len, type, frame):
# 	"""format_alignment(align1, align2, score, begin, end) -> string
# 	Format the alignment prettily into a string.
# 	"""
# 	s = []
# 	l1 = []
# 	l2 = []
# 	l3 = []
# 	
# 	al_l1 = list(align1)
# 	al_l2 = list(align2)
# 	element = 0
# 	
# 	
# 	# Simple alignments
# 	###################
# 	if edit_type == 'sim_align':
# 		del_l1 = []
# 		del_l2 = []
# 		# Handle the start of the sequence where the bases should align
# 		# Set the counter
# 		count = 0
# 		
# 		for base in al_l2:
# 			if base != al_l1[count]:
# 				break
# 			else:
# 				del_l2.append(al_l2[count])
# 				del_l1.append(al_l1[count]) 
# 				count = count + 1
# 			
# 		# ASSEMBLE THE ALIGNMENT
# 		# Sequences should be equal length at this stage
# 		if len(del_l1) == len(del_l2):
# 			for elements in del_l2:
# 				l1.append(del_l1[element])
# 				l3.append(del_l2[element])
# 				
# 				if del_l1[element] == del_l2[element]:
# 					l2.append("|")
# 				else:
# 					l2.append(" ")
# 				
# 				# Add 1 to the element counter	
# 				element = element +1
# 
# 	
# 	# Select for out of frame protein variants
# 	##########################################
# 	if edit_type == 'sim_prot':
# 		# Select for out of frame protein variants
# 		if frame == 'false':
# 			
# 			#return 'girffe'
# 			# Basic alignment
# 			del_l1 = []
# 			del_l2 = []
# 			
# 			# Handle the start of the sequence where the bases should align
# 			# Set the counter
# 			count = 0
# 			while al_l2:
# 				if al_l2[count] != al_l1[count]:
# 					break
# 				else:
# 					del_l2.append(al_l2[count])
# 					del_l1.append(al_l1[count]) 
# 					count = count + 1
# 			
# 			add1 = count
# 			add2 = count
# 			
# 			while add1 < len(al_l1):
# 				del_l1.append(al_l1[add1])
# 				add1 = add1 + 1
# 			
# 			while add2 < len(al_l2):
# 				del_l2.append(al_l2[add2])
# 				add2 = add2 + 1
# 			
# 			#return del_l1
# 			#return del_l2
# 			
# 			# ASSEMBLE THE ALIGNMENT
# 			# Sequences should be equal length at this stage
# 			if len(del_l1) >= len(del_l2):
# 				element = 0
# 				for elements in del_l2:
# 					l1.append(del_l1[element])
# 					l3.append(del_l2[element])
# 				
# 					if del_l1[element] == del_l2[element]:
# 						l2.append("|")
# 					else:
# 						l2.append(" ")
# 				
# 					# Add 1 to the element counter	
# 					element = element +1
# 			else:
# 				element = 0
# 				# return 'look, another giraffe!'
# 				for elements in del_l1:
# 					l1.append(del_l1[element])
# 					l3.append(del_l2[element])
# 				
# 					if del_l1[element] == del_l2[element]:
# 						l2.append("|")
# 					else:
# 						l2.append(" ")
# 				
# 					# Add 1 to the element counter	
# 					element = element +1
# 	
# 	# Align point mutations
# 	#######################
# 	# Simple alignment which need no alterations accounting for protein changes
# 	if edit_type == 'point':
# 	
# 		# Sequences should be equal length at this stage
# 		if len(al_l1) == len(al_l2):
# 			for elements in al_l1:
# 				l1.append(al_l1[element])
# 				l3.append(al_l2[element])
# 			
# 				if al_l1[element] == al_l2[element]:
# 					l2.append("|")
# 				else:
# 					l2.append(" ")
# 			
# 				# Add 1 to the element coounter
# 				element = element +1
# 		
# 	# Join together the elements into alignment strings
# 	s1 = ''.join(l1)
# 	s2 = ''.join(l2)
# 	s3 = ''.join(l3)
# 
# 
# 	# Align simple deletions
# 	########################
# 	if edit_type == 'sim_del':
# 	
# 		# Select for in frame protein or non protein variants
# 		if (frame == 'true') or (type == ':g.')	or (type == ':c.') or (type == ':n.'):
# 			
# 			# List of the sequence with the deletion
# 			del_l1 = []
# 			del_l2 = []
# 		
# 			# The gap has been set to end - start then + 1 (so 23 - 22 = 1, but gap actually - 2 so plus 1)
# 			gaps = interval_end - interval_start + 1
# 			# Set the counter
# 			count = 0
# 			while count < interval_start - begin -1:
# 				if al_l2[count] != al_l1[count]:
# 					break
# 				else:
# 					del_l2.append(al_l2[count])
# 					del_l1.append(al_l1[count]) 
# 					count = count + 1
# 			
# 			#return align2
# 			
# 			# reset the counter
# 			count = 0
# 			# Fill the gap with  - and append
# 			while count < gaps:
# 				base = "-"
# 				del_l1.append(base)
# 				count = count +1 			
# 			
# 			# To keep the next loop simple, we need to remove the extraneous 5 prime sequence
# 			# We cut off the length of the current edit list
# 			cut = len(del_l1)
# 		
# 			# Remove the processed section from each list
# 			del al_l1[0:cut - gaps]
# 			del_l2 = al_l2[:]
# 			del del_l2[0:cut]
# 		
# 			# Now repeat the loop on the truncated sequences to complete the alignments
# 			count = 0
# 			for base in al_l1:
# 				if  al_l1[count] != del_l2[count]:
# 					break
# 				else:
# 					del_l1.append(al_l1[count])
# 					count = count +1 
# 		
# 			# Sequences should be equal length at this stage
# 			if len(del_l1) == len(al_l2):
# 				for elements in del_l1:
# 					l1.append(del_l1[element])
# 					l3.append(al_l2[element])
# 				
# 					if del_l1[element] == al_l2[element]:
# 						l2.append("|")
# 					else:
# 						l2.append(" ")
# 				
# 					# Add 1 to the element counter	
# 					element = element +1
# 			#else:
# 				#return del_l2
# 				#l1.append('Giraffes')		
# 				#l1.append(del_l1)
# 	
# 
# 	# Align simple insertions and duplications
# 	##########################################
# 	if edit_type == 'sim_ins' or edit_type == 'dupn':
# 	
# 		# Select for in frame protein or non protein variants
# 		if (frame == 'true') or (type == ':g.')	or (type == ':c.') or (type == ':r.'):
# 			
# 			# List of the sequence with the deletion
# 			ins_l2 = []
# 			ins_l1 = al_l1[:]
# 			
# 			# Set the counter and gap length
# 			count = 0
# 			gaps = ins_len
# 			
# 			#return str(ins_len)
# 			
# 			if edit_type == 'dupn':
# 				# While the bases align, append the element to the edited list
# 				while count < interval_end - begin:
# 					if  al_l1[count] != al_l2[count]:
# 						break
# 					else:
# 						ins_l2.append(al_l2[count])
# 						count = count +1
# 			else:
# 				# While the bases align, append the element to the edited list
# 				while count < interval_start - begin:
# 					if  al_l1[count] != al_l2[count]:
# 						break
# 					else:
# 						ins_l2.append(al_l2[count])
# 						count = count +1
# 			
# 			# reset the counter
# 			count = 0
# 			# Fill the gap with  - and append
# 			while count < gaps:
# 				base = "-"
# 				ins_l2.append(base)
# 				count = count +1 			
# 			
# 			# To keep the next loop simple, we need to remove the extraneous 5 prime sequence
# 			# We cut off the length of the current edit list
# 			cut = len(ins_l2)
# 		
# 			# Remove the processed section from each list
# 			del al_l2[0:cut - gaps]
# 			del al_l1[0:cut]
# 			
# 			# Now repeat the loop on the truncated sequences to complete the alignments
# 			count = 0
# 			for base in al_l2:
# 				if  al_l1[count] != al_l2[count]:
# 					break
# 				else:
# 					ins_l2.append(al_l2[count])
# 					count = count +1 
# 		
# 			#return ins_l2
# 			
# 			# Sequences should be equal length at this stage
# 			if len(ins_l1) == len(ins_l2):
# 				for elements in ins_l1:
# 					l1.append(ins_l1[element])
# 					l3.append(ins_l2[element])
# 				
# 					if ins_l1[element] == ins_l2[element]:
# 						l2.append("|")
# 					else:
# 						l2.append(" ")
# 				
# 					# Add 1 to the element counter
# 					element = element +1
# 			#else:
# 				#return del_l2
# 				#l1.append('Giraffes')		
# 				#l1.append(del_l1)
# 	
# 
# 	# Align delins mutations
# 	########################
# 	if edit_type == 'delins':
# 		# List of the sequence with the deletion
# 		del_l1 = []
# 		del_l2 = []
# 		
# 		# return al_l2
# 		
# 		# Set the total gap length
# 		gaps = interval_end - interval_start + 1
# 		insertion = ins_len 
# 
# 		# return str(ins_len)
# 		
# 		# Select for in frame protein or non protein variants
# 		if (frame == 'true') or (type == ':g.')	or (type == ':c.') or (type == ':r.'):
# 			# Handle the start of the sequence where the bases should align
# 			# Set the counter
# 			count = 0
# 			while count < interval_start - begin -1:
# 				if al_l2[count] != al_l1[count]:
# 					break
# 				else:
# 					del_l2.append(al_l2[count])
# 					del_l1.append(al_l1[count]) 
# 					count = count + 1
# 			
# 			# Handle the deletion
# 			# append the deleted section to the reference list
# 			# blank out the variant base
# 			count = 0
# 			done = len(del_l2)
# 			while count < gaps:
# 				base = '-'
# 				del_l1.append(base)
# 				del_l2.append(al_l2[done])
# 				count = count + 1
# 				done = done + 1
# 			
# 			# Handle the insertion
# 			# reset the counter
# 			count = 0
# 			done = len(del_l2)
# 			while count < ins_len:
# 				base = '-'
# 				del_l2.append(base)
# 				del_l1.append(al_l1[done - gaps])
# 				count = count + 1
# 				done = done + 1
# 			
# 			# Remove the processed section from each list
# 			# Lists are of equal length
# 			done = len(del_l1)
# 			del al_l1[0:done - gaps]
# 			del al_l2[0:done - ins_len]
# 			
# 			# Now repeat the loop on the truncated sequences to complete the alignments
# 			count = 0
# 			for base in al_l1:
# 				if  del_l1[count] != del_l2[count]:
# 					break
# 				else:
# 					del_l1.append(al_l1[count])
# 					del_l2.append(al_l2[count])
# 					count = count +1 
# 		
# 			# ASSEMBLE THE ALIGNMENT
# 			# Sequences should be equal length at this stage
# 			if len(del_l2) == len(del_l1):
# 				for elements in del_l2:
# 					l1.append(del_l1[element])
# 					l3.append(del_l2[element])
# 				
# 					if del_l1[element] == del_l2[element]:
# 						l2.append("|")
# 					else:
# 						l2.append(" ")
# 				
# 					# Add 1 to the element counter	
# 					element = element +1
# 	
# 		# Select for out of frame protein variants
# 		if frame == 'false':
# 			# Basic alignment
# 
# 			# Handle the start of the sequence where the bases should align
# 			# Set the counter
# 			count = 0
# 			while al_l2:
# 				if al_l2[count] != al_l1[count]:
# 					break
# 				else:
# 					del_l2.append(al_l2[count])
# 					del_l1.append(al_l1[count]) 
# 					count = count + 1
# 			
# 			add1 = count
# 			add2 = count
# 			
# 			while add1 < len(al_l1):
# 				del_l1.append(al_l1[add1])
# 				add1 = add1 + 1
# 			
# 			while add2 < len(al_l2):
# 				del_l2.append(al_l2[add2])
# 				add2 = add2 + 1
# 			
# 			# ASSEMBLE THE ALIGNMENT
# 			if len(del_l1) >= len(del_l2):
# 				for elements in del_l2:
# 					l1.append(del_l1[element])
# 					l3.append(del_l2[element])
# 				
# 					if del_l1[element] == del_l2[element]:
# 						l2.append("|")
# 					else:
# 						l2.append(" ")
# 				
# 					# Add 1 to the element counter	
# 					element = element +1
# 			else:
# 				for elements in del_l1:
# 					l1.append(del_l1[element])
# 					l3.append(del_l2[element])
# 				
# 					if del_l1[element] == del_l2[element]:
# 						l2.append("|")
# 					else:
# 						l2.append(" ")
# 				
# 					# Add 1 to the element counter	
# 					element = element +1
# 	
# 	
# 	# String the alignment lists
# 	############################
# 	
# 	# Join together the elements into alignment strings
# 	s1 = ''.join(l1)
# 	s2 = ''.join(l2)
# 	s3 = ''.join(l3)
# 		
# 		
# 	# Append the alignment strings into a list (s)
# 	s.append(s1)
# 	s.append(s2)
# 	s.append(s3)
# 		
# 	return s

	
	
	
	
	
	
	
	
	
	
	
	
	
		
