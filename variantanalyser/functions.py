# -*- coding: utf-8 -*-
"""
functions.py
 
Module containing VariantValidator sub-functions. The majoirty of these functions require
hgvs Python package top-level functions or sub-functions contained in uta.py and 
seqfetcher.py 
"""

# Config Section Mapping function
def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

# IMPORT REQUIRED PYTHON MODULES
import re
import os
import sys
import copy
import requests

# Set up paths
# FUNCTIONS_ROOT = os.path.dirname(os.path.abspath(__file__))
ENTREZ_ID = os.environ.get('ENTREZ_ID')
if ENTREZ_ID is None:	
	from configparser import ConfigParser
	CONF_ROOT = os.environ.get('CONF_ROOT')
	Config = ConfigParser()
	Config.read(os.path.join(CONF_ROOT, 'config.ini'))
	ENTREZ_ID = ConfigSectionMap("EntrezID")['entrezid']
	
# IMPORT HGVS MODULES and create instances
import hgvs
import hgvs.exceptions
from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSUnsupportedOperationError
from hgvs.dataproviders import uta, seqfetcher
import hgvs.normalizer
import hgvs.validator
import hgvs.parser
import hgvs.variantmapper

# Connect to UTA
hdp = hgvs.dataproviders.uta.connect(pooling=True)
# Create normalizer
hn = hgvs.normalizer.Normalizer(hdp,
		cross_boundaries=False,
		shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
		alt_aln_method='splign'
		)	

# Validator
vr = hgvs.validator.Validator(hdp)	
# parser
hp = hgvs.parser.Parser()
# Variantmapper
vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=True) #, normalize=False)
# SeqFetcher
sf = hgvs.dataproviders.seqfetcher.SeqFetcher()

# variantanalyser modules
import dbControls
import supported_chromosome_builds

# BioPython
from Bio import Entrez
from Bio import SeqIO

# HGNC rest variables
import httplib2 as http
import json
try:
 	from urlparse import urlparse
except ImportError:
 	from urllib.parse import urlparse

"""
usr_input
collect the input from the form and convert to a hgvs readable string
	Removes brackets and contained information -if given
	Identifies variant type (p. c. etc)
	Returns a dictionary containing a formated input string which is optimal for hgvs 
	parsing and the variant type
	Accepts c, g, n, r currently. And now P also 15.07.15
"""
def user_input(input):
	raw_variant = input.strip()
	
	# Set regular expressions for if statements
	pat_g = re.compile("\:g\.")			# Pattern looks for :g.
	pat_gene = re.compile('\(.+?\)')	# Pattern looks for (....)
	pat_c = re.compile("\:c\.") 		# Pattern looks for :c.
	pat_r = re.compile("\:r\.") 		# Pattern looks for :r.
	pat_n = re.compile("\:n\.") 		# Pattern looks for :n.
	pat_p = re.compile("\:p\.") 		# Pattern looks for :p. 
	pat_m = re.compile("\:m\.")			# Pattern looks for :m.
	pat_est = re.compile("\d\:\d")		# Pattern looks for number:number
	
	# If statements
	if pat_g.search(raw_variant): # If the :g. pattern is present in the raw_variant, g_in is linked to the raw_variant
		if pat_gene.search(raw_variant):				# If pat gene is present in the raw_variant
			variant = pat_gene.sub('', raw_variant) 	# variant is set to the raw_variant string with the pattern (...) substituted out
			formated = {'variant' : variant, 'type' : ':g.'}
			return formated
		else:
			variant = raw_variant						# Otherwise it is set to raw_variant
			formated = {'variant' : variant, 'type' : ':g.'}
			return formated
	
	elif pat_r.search(raw_variant): 
		if pat_gene.search(raw_variant):				
			variant = pat_gene.sub('', raw_variant) 	
			formated = {'variant' : variant, 'type' : ':r.'}
			return formated
		else:
			variant = raw_variant						
			formated = {'variant' : variant, 'type' : ':r.'}
			return formated

	elif pat_n.search(raw_variant): 
		if pat_gene.search(raw_variant):				
			variant = pat_gene.sub('', raw_variant) 	
			formated = {'variant' : variant, 'type' : ':n.'}
			return formated
		else:
			variant = raw_variant						
			formated = {'variant' : variant, 'type' : ':n.'}
			return formated
	
	elif pat_c.search(raw_variant): 
		if pat_gene.search(raw_variant):				
			variant = pat_gene.sub('', raw_variant) 	
			formated = {'variant' : variant, 'type' : ':c.'}
			return formated
		else:
			variant = raw_variant						
			formated = {'variant' : variant, 'type' : ':c.'}
			return formated

	elif pat_p.search(raw_variant): 
		variant = raw_variant						
		formated = {'variant' : variant, 'type' : ':p.'}
		return formated
	
	elif pat_m.search(raw_variant): 
		variant = raw_variant						
		formated = {'variant' : variant, 'type' : ':m.'}
		return formated
	elif pat_est.search(raw_variant): 
		variant = raw_variant						
		formated = {'variant' : variant, 'type' : 'est'}
		return formated
	else:
		formatted = 'invalid'
		return formatted
		
"""
r_to_c
parses r. variant strings into hgvs object and maps to the c. equivalent. 
""" 
def r_to_c(variant, evm, hp):
	# convert the input string into a hgvs object by parsing
	var_r = hp.parse_hgvs_variant(variant)	
	# map to the coding sequence
	var_c = evm.r_to_c(var_r)	#  coding level variant
	variant = str(var_c)
	c_from_r = {'variant' : variant, 'type' : ':c.'}
	return c_from_r
	
"""	
Maps transcript variant descriptions onto specified RefSeqGene reference sequences
Return an hgvs object containing the genomic sequence variant relative to the RefSeqGene 
acession
refseq_ac = RefSeqGene ac
"""
def refseq(variant, vm, refseq_ac, hp, evm, hdp, primary_assembly):
	vr = hgvs.validator.Validator(hdp)
	# parse the variant into hgvs object
	var_c = hp.parse_hgvs_variant(variant)
	# map to the genomic co-ordinates using the easy variant mapper set to alt_aln_method = alt_aln_method
	var_g = myevm_t_to_g(var_c, evm, hdp, primary_assembly)
	# Get overlapping transcripts - forcing a splign alignment
	start_i = var_g.posedit.pos.start.base
	end_i = var_g.posedit.pos.end.base
	alt_ac = var_g.ac
	alt_aln_method = 'splign'
	transcripts = hdp.get_tx_for_region(alt_ac,alt_aln_method,start_i-1,end_i)
	# Take the first transcript
	for trans in transcripts:
		tx_ac = trans[0]
		try:
			ref_c = vm.g_to_t(var_g, tx_ac, alt_aln_method='splign')
		except:
			continue
		else:
			# map the variant co-ordinates to the refseq Gene accession using vm
			ref_g_dict = {
						'ref_g' : '',
						'error' : 'false'
						}
			try:	
				ref_g_dict['ref_g'] = vm.t_to_g(ref_c, alt_ac=refseq_ac, alt_aln_method='splign')
			except:
				e = sys.exc_info()[0]
				ref_g_dict['error'] = e
			try: 
				vr.validate(ref_g_dict['ref_g'])
			except:	
				e = sys.exc_info()[0]
				ref_g_dict['error'] = e			
			if ref_g_dict['error'] == 'false':
				return ref_g_dict
			else:
				continue	
	# Return as an error if all fail
	return ref_g_dict
	

"""
Parses genomic variant strings into hgvs objects
Maps genomic hgvs object into a coding hgvs object if the c accession string is provided
returns a c. variant description string
"""
def g_to_c(var_g, tx_ac, hp, evm):
	pat_g = re.compile("\:g\.") 		# Pattern looks for :g.
	# If the :g. pattern is present in the input variant
	if pat_g.search(var_g): 
		# convert the input string into a hgvs object by parsing
		var_g = hp.parse_hgvs_variant(var_g)
		# Map to coding variant
		var_c = str(evm.g_to_c(var_g, tx_ac))
		return var_c
		

"""
Parses genomic variant strings into hgvs objects
Maps genomic hgvs object into a non-coding hgvs object if the n accession string is provided
returns a n. variant description string
"""
def g_to_n(var_g, tx_ac, hp, evm):
	pat_g = re.compile("\:g\.") 		# Pattern looks for :g.
	# If the :g. pattern is present in the input variant
	if pat_g.search(var_g): 
		# convert the input string into a hgvs object by parsing
		var_g = hp.parse_hgvs_variant(var_g)
		# Map to coding variant
		var_n = str(evm.g_to_n(var_g, tx_ac))
		return var_n


"""
Ensures variant strings are transcript c. or n.
returns parsed hgvs c. or n. object
"""
def coding(variant, hp):
	# If the :c. pattern is present in the input variant
	if re.search(':c.', variant) or re.search(':n.', variant): 
		# convert the input string into a hgvs object
		var_c = hp.parse_hgvs_variant(variant)
		return var_c
		

"""
Mapping transcript to genomic position
Ensures variant strings are transcript c. or n.
returns parsed hgvs g. object
"""
def genomic(variant, evm, hp, hdp, primary_assembly):
	# Set regular expressions for if statements
	pat_g = re.compile("\:g\.")			# Pattern looks for :g.
	pat_n = re.compile("\:n\.")
	pat_c = re.compile("\:c\.") 		# Pattern looks for :c.
	
	# If the :c. pattern is present in the input variant
	if  pat_c.search(variant) or pat_n.search(variant): 
		error = 'false'
		hgvs_var = hp.parse_hgvs_variant(variant)
		try:
			var_g = myevm_t_to_g(hgvs_var, evm, hdp, primary_assembly)	# genomic level variant
		except hgvs.exceptions.HGVSError as e:
			error = e
		if error != 'false':
			var_g = 'error ' + str(e)
		return var_g

	# If the :g. pattern is present in the input variant
	elif  (pat_g.search(variant)): # or (pat_n.search(variant)): 
		# convert the input string into a hgvs object
		var_g = hp.parse_hgvs_variant(variant)
		return var_g

"""
Mapping transcript to protein prediction
Ensures variant strings are transcript c.
returns parsed hgvs p. object
"""
def protein(variant, evm, hp):
	variant = str(variant)
	# Set regular expressions for if statements
	pat_c = re.compile("\:c\.") 		# Pattern looks for :c. Note (gene) has been removed
	
	# If the :c. pattern is present in the input variant
	if  pat_c.search(variant): 
		# convert the input string into a hgvs object
		var_c = hp.parse_hgvs_variant(variant)		
		# map to the genomic sequence
		var_p = evm.c_to_p(var_c)	# genomic level variant
		return var_p
	if re.search(':n.', variant):
		var_p = hp.parse_hgvs_variant(variant)
		var_p.ac = 'Non-coding transcript '
		var_p.posedit = ''
		return var_p
	

"""
Marked for removal
"""
# Return an hgvs object containing the rna sequence variant
# def rna(variant, evm, hp):
# 	Set regular expressions for if statements
# 	pat_c = re.compile("\:c\.") 		# Pattern looks for :c. Note (gene) has been removed
# 	If the :c. pattern is present in the input variant
# 	if  pat_c.search(variant): 
# 		convert the input string into a hgvs object
# 		var_c = hp.parse_hgvs_variant(variant)
# 		map to the genomic sequence
# 		var_r = evm.c_to_n(var_c)	# rna level variant
# 		return var_r

"""
Marked for removal
"""
# def hgvs_rna(variant, hp):
# 	# Set regular expressions for if statements
# 	pat_r = re.compile("\:n\.") 		# Pattern looks for :n. Note (gene) has been removed
# 	# If the :r. pattern is present in the input variant
# 	if  pat_r.search(variant): 
# 		# convert the input string into a hgvs object
# 		var_r = hp.parse_hgvs_variant(variant)
# 		return var_r


"""
Ensures variant strings are g.
returns parsed hgvs g. object
"""
def hgvs_genomic(variant, hp):
	# Set regular expressions for if statements
	pat_g = re.compile("\:g\.") 		# Pattern looks for :g. Note (gene) has been removed
	# If the :g. pattern is present in the input variant
	if  pat_g.search(variant): 
		# convert the input string into a hgvs object
		var_g = hp.parse_hgvs_variant(variant)
		return var_g


"""
Enhanced transcript to genome position mapping function using evm
Deals with mapping from transcript positions that do not exist in the genomic sequence
i.e. the stated position aligns to a genomic gap!
Trys to ensure that a genomic position is always returned even if the c. or n. transcript
will not map to the specified genome build primary assembly.
Deals with transcript mapping to several genomic assemblies
Order 
Map to a single NC_ for the specified genome build primary assembly
Map to a single NC_ for an alternate genome build primary assembly
Map to an NT_ from the specified genome build
Map to an NT_ from an alternative genome build
Map to an NW_ from the specified genome build
Map to an NW_ from an alternative genome buildRequires parsed c. or n. object
returns parsed hgvs g. object
"""
def myevm_t_to_g(hgvs_c, evm, hdp, primary_assembly):

# 	create no_norm_evm
 	no_norm_evm = hgvs.assemblymapper.AssemblyMapper(hdp,
 			assembly_name=primary_assembly, 
 			alt_aln_method='splign', 
 			normalize=False, 
 			replace_reference=True
 			)
		
	# store the input
	stored_hgvs_c = copy.deepcopy(hgvs_c)
	
	# Set expansion variable
	expand_out = 'false'
	if hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del' or hgvs_c.posedit.edit.type =='delins' or hgvs_c.posedit.edit.type == 'dup' or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins':
				
		# if NM_ need the n. position
		if re.match('NM_', str(hgvs_c.ac)):
			hgvs_c = no_norm_evm.c_to_n(hgvs_c)
		# Check for intronic
		try: 
			hn.normalize(hgvs_c)
		except hgvs.exceptions.HGVSError as e:
			error = str(e)
			if re.search('intronic variant', error):
				pass
			elif re.search('Length implied by coordinates must equal sequence deletion length', error) and hgvs_c.type == 'n':
				hgvs_c.posedit.pos.end.base = hgvs_c.posedit.pos.start.base + len(hgvs_c.posedit.edit.ref) - 1
		
		# Check again before continuing
		if re.search('\d+\+', str(hgvs_c.posedit.pos)) or re.search('\d+\-', str(hgvs_c.posedit.pos)) or re.search('\*\d+\+', str(hgvs_c.posedit.pos)) or re.search('\*\d+\-', str(hgvs_c.posedit.pos)):
			pass		

		else:		
			try:	
				# For non-intronic sequence
				hgvs_t = copy.deepcopy(hgvs_c)
				if hgvs_c.posedit.edit.type == 'dup':
					# hgvs_t = reverse_normalize.normalize(hgvs_t)
					pre_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
					post_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
					alt = pre_base + hgvs_t.posedit.edit.ref + hgvs_t.posedit.edit.ref + post_base				
					ref = pre_base + hgvs_t.posedit.edit.ref + post_base
					dup_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base - 1) + '_' + str((hgvs_t.posedit.pos.start.base + len(ref)) -2) + 'del' + ref + 'ins' + alt
					hgvs_t = hp.parse_hgvs_variant(dup_to_delins)
				elif hgvs_c.posedit.edit.type == 'ins':	
					pass
				else:	
					if str(hgvs_t.posedit.edit.alt) == 'None':
						hgvs_t.posedit.edit.alt = ''	
					pre_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
					post_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
					hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
					hgvs_t.posedit.edit.alt = pre_base + hgvs_t.posedit.edit.alt + post_base
					hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
					start = hgvs_t.posedit.pos.start.base
					hgvs_t.posedit.pos.start.base = start + 1
					hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
					end = hgvs_t.posedit.pos.end.base
					hgvs_t.posedit.pos.start.base = start
					hgvs_t.posedit.pos.end.base = end
					hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(end) + str(hgvs_t.posedit.edit)
					hgvs_t = hp.parse_hgvs_variant(hgvs_str)
				hgvs_c = copy.deepcopy(hgvs_t)
			
				# Set expanded out test to true
				expand_out = 'true'

			except Exception:
				hgvs_c = hgvs_c			
			
		if re.match('NM_', str(hgvs_c.ac)):
			try:
				hgvs_c = no_norm_evm.n_to_c(hgvs_c)
			except hgvs.exceptions.HGVSError as e:
				hgvs_c = copy.deepcopy(stored_hgvs_c)				
				
							
		# Ensure the altered c. variant has not crossed intro exon boundaries
		hgvs_check_boundaries = copy.deepcopy(hgvs_c)
		try:
			h_variant = hn.normalize(hgvs_check_boundaries)
		except hgvs.exceptions.HGVSError as e:
			error = str(e) 
			if re.search('spanning the exon-intron boundary', error):
				hgvs_c = copy.deepcopy(stored_hgvs_c)
			else:
				hgvs_c = copy.deepcopy(stored_hgvs_c)
				
	try:
		hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
		hn.normalize(hgvs_genomic) # Check the validity of the mapping
		# This will fail on multiple refs for NC_
	except hgvs.exceptions.HGVSError as e:

		# Recover all available mapping options from UTA
		mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)	
		if mapping_options == []:
			raise HGVSDataNotAvailableError("no g. mapping options available")
		for option in mapping_options:
			if re.match('NC_', option[1]):
				chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
				if chr_num != 'false':
					try:
						hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
						break
					except:
						continue

		try:
			hn.normalize(hgvs_genomic)
		except:	
			for option in mapping_options:
				if re.match('NC_', option[1]):
					try:
						hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
						break
					except:
						continue
			try:
				hn.normalize(hgvs_genomic)
			except:
				for option in mapping_options:
					if re.match('NT_', option[1]):
						try:
							hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
							break
						except:
							continue
				try:
					hn.normalize(hgvs_genomic)
				except:							
					for option in mapping_options:
						if re.match('NW_', option[1]):
							try:
								hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
								break
							except:
								continue
		 
	if hgvs_genomic.posedit.edit.type == 'ins':
		try:
			hgvs_genomic = hn.normalize(hgvs_genomic)
		except hgvs.exceptions.HGVSError as e:
			error = str(e)
			if error == 'insertion length must be 1':
				ref = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.end.base)
				hgvs_genomic.posedit.edit.ref = ref
				hgvs_genomic.posedit.edit.alt = ref[0:1] + hgvs_genomic.posedit.edit.alt + ref[-1:]
				hgvs_genomic = hn.normalize(hgvs_genomic)
				
	# Remove identity bases
	if hgvs_c == stored_hgvs_c:
		expanded_out = 'false'
	elif expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) >= 3:
		hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
		hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
		hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
		if hgvs_genomic.posedit.edit.alt is not None:
			hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1] 
	else:
		pass		
		
	return hgvs_genomic

"""
USE WITH MAPPER THAT DOES NOT REPLACE THE REFERENCE GENOMIC BASES AND DOED NOT NORMALIZE

Enhanced transcript to genome position mapping function using evm
Trys to ensure that a genomic position is always returned even if the c. or n. transcript
will not map to the specified genome build primary assembly.
Deals with transcript mapping to several genomic assemblies
Order 
Map to a single NC_ (or ALT) for the specified genome build
returns parsed hgvs g. object
"""
def noreplace_myevm_t_to_g(hgvs_c, evm, hdp, primary_assembly):
	try:
		hgvs_genomic = evm.t_to_g(hgvs_c)
		# This will fail on multiple refs for NC_
	except hgvs.exceptions.HGVSError:
		# Recover all available mapping options from UTA
		mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)	
		if mapping_options == []:
			raise HGVSDataNotAvailableError("no g. mapping options available")
		for option in mapping_options:
			if re.match('NC_', option[1]):
				chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
				if chr_num != 'false':
					try:
						hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
						break
					except Exception as e:
						continue
		try:
			hgvs_genomic
		except:
			for option in mapping_options:
				if re.match('NC_', option[1]):
					chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
					try:
						hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
						break
					except Exception as e:
						continue

			try:
				hgvs_genomic
			except:
				for option in mapping_options:
					if re.match('NT_', option[1]):
						chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
						if chr_num != 'false':
							try:
								hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
								break
							except Exception as e:
								continue
				try:
					hgvs_genomic
				except:							
					for option in mapping_options:
						if re.match('NW_', option[1]):
							chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
							if chr_num != 'false':
								try:
									hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
									break
								except Exception as e:
									continue
					# All have failed, likely a variant description error
					try:
						hgvs_genomic
					except:
						for option in mapping_options:
							if re.match('NC_', option[1]):
								chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
								# Do not trap the error
								hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
								break
	return hgvs_genomic

"""
Enhanced transcript to genome position on a specified genomic reference using vm
Deals with mapping from transcript positions that do not exist in the genomic sequence
i.e. the stated position aligns to a genomic gap!
returns parsed hgvs g. object
"""
def myvm_t_to_g(hgvs_c, alt_chr, vm, hn, hdp, primary_assembly):

	# create no_norm_evm
	no_norm_evm = hgvs.assemblymapper.AssemblyMapper(hdp,
			assembly_name=primary_assembly, 
			alt_aln_method='splign', 
			normalize=False, 
			replace_reference=True
			)		
	
	# store the input
	stored_hgvs_c = copy.deepcopy(hgvs_c)
	
	# Set expansion variable
	expand_out = 'false'
	if hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del' or hgvs_c.posedit.edit.type =='delins' or hgvs_c.posedit.edit.type == 'dup' or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins':
		
		# if NM_ need the n. position
		if re.match('NM_', str(hgvs_c.ac)):
			hgvs_c = no_norm_evm.c_to_n(hgvs_c)
		
		# Check for intronic
		try: 
			hn.normalize(hgvs_c)
		except hgvs.exceptions.HGVSError as e:
			error = str(e)
			if re.search('intronic variant', error):
				pass
			elif re.search('Length implied by coordinates must equal sequence deletion length', error) and hgvs_c.type == 'n':
				hgvs_c.posedit.pos.end.base = hgvs_c.posedit.pos.start.base + len(hgvs_c.posedit.edit.ref) - 1

		# Check again before continuing
		if re.search('\d+\+', str(hgvs_c.posedit.pos)) or re.search('\d+\-', str(hgvs_c.posedit.pos)) or re.search('\*\d+\+', str(hgvs_c.posedit.pos)) or re.search('\*\d+\-', str(hgvs_c.posedit.pos)):
			pass

		else:		
			try:
				# For non-intronic sequence
				hgvs_t = copy.deepcopy(hgvs_c)
				if hgvs_c.posedit.edit.type == 'dup':
					# hgvs_t = reverse_normalize.normalize(hgvs_t)
					pre_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
					post_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
					alt = pre_base + hgvs_t.posedit.edit.ref + hgvs_t.posedit.edit.ref + post_base				
					ref = pre_base + hgvs_t.posedit.edit.ref + post_base
					dup_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base - 1) + '_' + str((hgvs_t.posedit.pos.start.base + len(ref)) -2) + 'del' + ref + 'ins' + alt
					hgvs_t = hp.parse_hgvs_variant(dup_to_delins)
				elif hgvs_c.posedit.edit.type == 'ins':	
					pass
				else:	
					if str(hgvs_t.posedit.edit.alt) == 'None':
						hgvs_t.posedit.edit.alt = ''	
					pre_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
					post_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
					hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
					hgvs_t.posedit.edit.alt = pre_base + hgvs_t.posedit.edit.alt + post_base
					hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
					start = hgvs_t.posedit.pos.start.base
					hgvs_t.posedit.pos.start.base = start + 1
					hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
					end = hgvs_t.posedit.pos.end.base
					hgvs_t.posedit.pos.start.base = start
					hgvs_t.posedit.pos.end.base = end
					hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(end) + str(hgvs_t.posedit.edit)
					hgvs_t = hp.parse_hgvs_variant(hgvs_str)
				hgvs_c = copy.deepcopy(hgvs_t)
			
				# Set expanded out test to true
				expand_out = 'true'
			
			except Exception:
				hgvs_c = hgvs_c		
		
		if re.match('NM_', str(hgvs_c.ac)):
			try:
				hgvs_c = no_norm_evm.n_to_c(hgvs_c)
			except hgvs.exceptions.HGVSError as e:
				hgvs_c = copy.deepcopy(stored_hgvs_c)
		
		# Ensure the altered c. variant has not crossed intro exon boundaries
		hgvs_check_boundaries = copy.deepcopy(hgvs_c)
		try:
			h_variant = hn.normalize(hgvs_check_boundaries)
		except hgvs.exceptions.HGVSError as e:
			error = str(e) 
			if re.search('spanning the exon-intron boundary', error):
				hgvs_c = copy.deepcopy(stored_hgvs_c)

	hgvs_genomic = vm.t_to_g(hgvs_c, alt_chr)	
	if hgvs_genomic.posedit.edit.type == 'ins':
		try:
			hgvs_genomic = hn.normalize(hgvs_genomic)
		except hgvs.exceptions.HGVSError as e:
			error = str(e)
			if error == 'insertion length must be 1':
				ref = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.end.base)
				hgvs_genomic.posedit.edit.ref = ref
				hgvs_genomic.posedit.edit.alt = ref[0:1] + hgvs_genomic.posedit.edit.alt + ref[-1:]
				hgvs_genomic = hn.normalize(hgvs_genomic)

	# Remove identity bases
	if hgvs_c == stored_hgvs_c:
		expanded_out = 'false'
	elif expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) >= 3:
		hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
		hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
		hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
		if hgvs_genomic.posedit.edit.alt is not None:
			hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1]  
	else:
		pass

	return hgvs_genomic
	

"""
Simple hgvs g. to c. or n. mapping
returns parsed hgvs c. or n. object
"""
def myevm_g_to_t(hdp, evm, hgvs_genomic, alt_ac):
	hgvs_t = evm.g_to_t(hgvs_genomic, alt_ac)					
	return hgvs_t
	
"""
parse p. strings into hgvs p. objects
"""
def hgvs_protein(variant, hp):
	# Set regular expressions for if statements
	pat_p = re.compile("\:p\.") 		# Pattern looks for :g. Note (gene) has been removed
	# If the :p. pattern is present in the input variant
	if  pat_p.search(variant): 
		# convert the input string into a hgvs object
		var_p = hp.parse_hgvs_variant(variant)
		return var_p

"""
Convert r. into c.
"""
def hgvs_r_to_c(hgvs_object):
	hgvs_object.type = 'c'
	edit = str(hgvs_object.posedit.edit)
	edit = edit.upper()
	# lowercase the supported variant types
	edit = edit.replace('DEL', 'del')
	edit = edit.replace('INS', 'ins')
	edit = edit.replace('INV', 'inv')
	edit = edit.replace('DUP', 'dup')
	# edit = edit.replace('CON', 'con')
	# edit = edit.replace('TRA', 'tra')
	edit = edit.replace('U', 'T')
	hgvs_object.posedit.edit = edit
	return hgvs_object

"""
Convert c. into r.
"""
def hgvs_c_to_r(hgvs_object):
	hgvs_object.type = 'r'
	edit = str(hgvs_object.posedit.edit)
	edit = edit.lower()
	edit = edit.replace('t', 'u')
	hgvs_object.posedit.edit = edit
	return hgvs_object	

"""
Input c. r. n. variant string
Use uta.py (hdp) to return the identity information for the transcript variant 
see hgvs.dataproviders.uta.py for details
"""
def tx_identity_info(variant, hdp):
	# Set regular expressions for if statements
	pat_c = re.compile("\:c\.") 		# Pattern looks for :c. Note (gene) has been removed
	pat_n = re.compile("\:n\.") 		# Pattern looks for :c. Note (gene) has been removed
	pat_r = re.compile("\:r\.") 		# Pattern looks for :c. Note (gene) has been removed
	
	# If the :c. pattern is present in the input variant
	if  pat_c.search(variant):
		# Remove all text to the right and including pat_c
		tx_ac = variant[:variant.index(':c.') + len(':c.')]
		tx_ac = pat_c.sub('', tx_ac)
		# Interface with the UTA database via get_tx_identity in uta.py
		tx_id_info = hdp.get_tx_identity_info(tx_ac)
		# NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
		return tx_id_info
	
	# If the :n. pattern is present in the input variant
	if pat_n.search(variant):
		# Remove all text to the right and including pat_c
		tx_ac = variant[:variant.index(':n.') + len(':n.')]
		tx_ac = pat_n.sub('', tx_ac)
		# Interface with the UTA database via get_tx_identity in uta.py
		tx_id_info = hdp.get_tx_identity_info(tx_ac)
		# NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
		return tx_id_info
		
	# If the :r. pattern is present in the input variant
	if pat_r.search(variant):
		# Remove all text to the right and including pat_c
		tx_ac = variant[:variant.index(':r.') + len(':r.')]
		tx_ac = pat_r.sub('', tx_ac)
		# Interface with the UTA database via get_tx_identity in uta.py
		tx_id_info = hdp.get_tx_identity_info(tx_ac)
		# NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
		return tx_id_info
		
"""
Input c. r. nd accession string
Use uta.py (hdp) to return the identity information for the transcript variant 
see hgvs.dataproviders.uta.py for details
"""
def tx_id_info(alt_ac, hdp):
	tx_id_info = hdp.get_tx_identity_info(alt_ac)
	# NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
	return tx_id_info

	
"""
Use uta.py (hdp) to return the transcript information for a specified gene (HGNC SYMBOL)
see hgvs.dataproviders.uta.py for details
"""
def tx_for_gene(hgnc, hdp):
	# Interface with the UTA database via get_tx_for_gene in uta.py
	tx_for_gene = hdp.get_tx_for_gene(hgnc)
	return tx_for_gene 
	

"""
Extract RefSeqGene Accession from transcript information
see hgvs.dataproviders.uta.py for details
"""
def ng_extract(tx_for_gene):
	# Set regular expressions for if statements
	pat_NG = re.compile("^NG_")			# Pattern looks for NG_ at beginning of a string
	# For each list in the list of lists tx_for_gene
	for list in tx_for_gene:
		# If the pattern NG_ is found in element 4
		if pat_NG.search(list[4]):
			# The gene accession is set to list element 4 
			gene_ac = list[4]
			return gene_ac
	
"""
marked for removal
"""
# def int_start(var_g):
# 	start = var_g.posedit.pos.start
# 	# Stringify to get start co-ords
# 	start = str(start)
# 	# Make into an integer
# 	int_start = int(start)
# 	return int_start

"""
marked for removal
"""
# def int_end(var_g):
# 	end = var_g.posedit.pos.end
# 	# Stringify to get start co-ords
# 	end = str(end)
# 	# Make into an integer
# 	int_end = int(end)
# 	return int_end
	
"""
Returns exon information for a given transcript
e.g. how the exons align to the genomic reference
see hgvs.dataproviders.uta.py for details
"""
def tx_exons(tx_ac, alt_ac, alt_aln_method, hdp):
	
	# Interface with the UTA database via get_tx_exons in uta.py
	try:
		tx_exons = hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
	except hgvs.exceptions.HGVSError as e:
		e
		tx_exons = 'hgvs Exception: ' + str(e)
		return tx_exons
	try:
		completion = tx_exons[0]['alt_strand']
	except TypeError:
		tx_exons = 'error'
		return tx_exons
	# If on the reverse strand, reverse the order of elements
	if tx_exons[0]['alt_strand'] == -1:
		tx_exons = tx_exons[::-1]
		return tx_exons
	else:
		return tx_exons

"""
Automatically maps genomic positions onto all overlapping transcripts
"""
def relevant_transcripts(hgvs_genomic, evm, hdp, alt_aln_method):
	# Pass relevant transcripts for the input variant to rts
	rts = evm.relevant_transcripts(hgvs_genomic)
	
	# Project genomic variants to new transcripts 
	# and  populate a code_var list
	#############################################
	# Open a list to store relevant transcripts
	code_var = []
	# Populate transcripts - The keys become the list elements from rel_trs
	for x in rts:
		y = x.rstrip()	# Chomp any whitespace from the right of x ($_) - Assign to y
		# Easy variant mapper used to map the input variant to the relevant transcripts
		# Check for coding transcripts
		try:
			variant = evm.g_to_t(hgvs_genomic, y)
		except hgvs.exceptions.HGVSError as e:
			# Check for non-coding transcripts
			try:
				variant = evm.g_to_t(hgvs_genomic, y)
			except hgvs.exceptions.HGVSError as e:
				continue
		except:
			continue		
			
		# Corrective Normalisation of intronic descriptions in the antisense oriemtation
		pl = re.compile('\+')
		mi = re.compile('\-')
		ast = re.compile('\*')
		if pl.search(str(variant)) or mi.search(str(variant)) or ast.search(str(variant)):
			tx_ac = variant.ac
			alt_ac = hgvs_genomic.ac

			# Interface with the UTA database via get_tx_exons in uta.py
			try:
				tx_exons = hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
			except hgvs.exceptions.HGVSError as e:
				e
				tx_exons = 'hgvs Exception: ' + str(e)
				return tx_exons
			try:
				completion = tx_exons[0]['alt_strand']
			except TypeError:
				tx_exons = 'error'
				return tx_exons
			# If on the reverse strand, reverse the order of elements
			if tx_exons[0]['alt_strand'] == -1:
				tx_exons = tx_exons[::-1]
			else:
				pass
			
			# Gene orientation
			if tx_exons[0]['alt_strand'] == -1: 
				antisense = 'true'
			else:
				antisense = 'false'
				
			# Pass if antisense = 'false'
			if antisense == 'false':
				pass
			else:
				reverse_hn = hgvs.normalizer.Normalizer(hdp,
				cross_boundaries=False,
				shuffle_direction=5,
				alt_aln_method=alt_aln_method
				)
				# Reverse normalize hgvs_genomic
				rev_hgvs_genomic = reverse_hn.normalize(hgvs_genomic)
				# map back to coding
				variant = evm.g_to_t(rev_hgvs_genomic, tx_ac)
		code_var.append(str(variant))
	return code_var


"""
Take HGVS string, parse into hgvs object and validate
"""					
def validate(input, hp, vr):
	hgvs_input = hp.parse_hgvs_variant(input)
	g = re.compile(":g.")
	p = re.compile(":p.")
	if p.search(input):
		if hasattr(hgvs_input.posedit.pos.start, 'offset'):
			pass
		else:
			hgvs_input.posedit.pos.start.offset = 0
		if hasattr(hgvs_input.posedit.pos.end, 'offset'):
			pass
		else:
			hgvs_input.posedit.pos.end.offset = 0
		if hasattr(hgvs_input.posedit.pos.start, 'datum'):
			pass
		else:
			hgvs_input.posedit.pos.start.datum = 0
		if hasattr(hgvs_input.posedit.pos.end, 'datum'):
			pass
		else:
			hgvs_input.posedit.pos.end.datum = 0
		if hasattr(hgvs_input.posedit.edit, 'ref_n'):
			pass
		else:
			hgvs_input.posedit.edit.ref_n = hgvs_input.posedit.pos.end.base - hgvs_input.posedit.pos.start.base +1

	try:
		vr.validate( hgvs_input )
	except hgvs.exceptions.HGVSError as e:
	
		error = e
		return error
	
	else:
		error = 'false'
		return error 

"""
marked for removal
"""
# def sequence_extractor(ac, hdp):
# 	ac_seq = hdp.get_tx_seq(ac)
# 	return ac_seq

"""
marked for removal
""" 
# def ref_replace(e, hgvs_variant):
# 	error = str(e)
# 	match = re.findall('\(([GATC]+)\)', error)
# 	new_ref = match[1] 
# 	hgvs_variant.posedit.edit.ref = new_ref
# 	return hgvs_variant	
	
"""
Search HGNC rest
"""	
def hgnc_rest(path):	
	data = {
		'record' : '',
		'error' : 'false'
		}
	# HGNC server
	headers = {
 		'Accept': 'application/json',
		}
	uri = 'http://rest.genenames.org'
	target = urlparse(uri+path)
	method = 'GET'
	body = ''
	h = http.Http()
	# collect the response
	response, content = h.request(
 		target.geturl(),
	 	method,
 		body,
 		headers)
	if response['status'] == '200':
		# assume that content is a json reply
 		# parse content with the json module 
 		data['record'] = json.loads(content)	
	else:
		data['error'] = "Unable to contact the HGNC database: Please try again later"
	return data
	
	
"""
Search Entrez databases with efetch and SeqIO
"""
def entrez_efetch(db, id, rettype, retmode):			
	# IMPORT Bio modules
	#from Bio import Entrez
	Entrez.email = ENTREZ_ID
	#from Bio import SeqIO
	handle = Entrez.efetch(db=db, id=id, rettype=rettype, retmode=retmode)
	# Get record
	record = SeqIO.read(handle, "gb")
	# Place into text
	# text = handle.read()
	handle.close()
	return record
	
	
"""
search Entrez databases with efetch and read
"""
def entrez_read(db, id,retmode):	
	# IMPORT Bio modules
	#from Bio import Entrez
	Entrez.email = ENTREZ_ID
	#from Bio import SeqIO
	handle = Entrez.efetch(db=db, id=id, retmode=retmode)
	# Get record
	record = Entrez.read(handle)
	# Place into text
	# text = handle.read()
	handle.close()
	return record
	
"""
Simple reverse complement function for nucleotide sequences
"""
def revcomp(bases):
	l2 = []
	l = list(bases)
	element = 0
	for base in l:
		element = element+1
		if base == 'G':
			l2.append('C')
		if base == 'C':
			l2.append('G')
		if base == 'A':
			l2.append('T')
		if base == 'T':
			l2.append('A')
	revcomp = ''.join(l2)
	revcomp = revcomp[::-1]
	return revcomp
	
# <LICENSE>

# </LICENSE>	
