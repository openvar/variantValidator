# -*- coding: utf-8 -*-

# Module containing functions that use hgvs to return variant data

# IMPORT REQUIRED PYTHON MODULES
import re
import os
import copy

# Set up paths
FUNCTIONS_ROOT = os.path.dirname(os.path.abspath(__file__))

# IMPORT HGVS MODULES
import hgvs
import hgvs.exceptions
# Create SeqFetcher
import hgvs.dataproviders
sf = hgvs.dataproviders.seqfetcher.SeqFetcher()							
import hgvs.normalizer
import hgvs.validator
import hgvs.parser
hp = hgvs.parser.Parser()
		
# Enables Ensembl Rest variables
import requests, sys
# HGNC rest variables
import httplib2 as http
import json
try:
 	from urlparse import urlparse
except ImportError:
 	from urllib.parse import urlparse

# XML variables
# import xml.etree.ElementTree as ET

# usr_input function: collect the input from the form and convert to a hgvs readable string
	# Removes gene name (if given)
	# Identifies variant type
	# Returns a dictionary containing the formated input string and the variant type
	# Accepts c, g, n, r currently. And now P also 15.07.15
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
		#if pat_gene.search(raw_variant):				
		#	variant = pat_gene.sub('', raw_variant) 	
		#	formated = {'variant' : variant, 'type' : ':p.'}
		#	return formated
		#else:
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

# Maps the r variant to the c variant
def r_to_c(variant, evm, hp):
	# convert the input string into a hgvs object by parsing
	var_r = hp.parse_hgvs_variant(variant)	
	# map to the coding sequence
	var_c = evm.r_to_c(var_r)	#  coding level variant
	variant = str(var_c)
	c_from_r = {'variant' : variant, 'type' : ':c.'}
	return c_from_r
	

# Return an hgvs object containing the genomic sequence variant relative to the refseq acession
# This approach is required to handle alt_aln_method other than splign
# Return an hgvs object containing the genomic sequence variant relative to the refseq acession
# This approach is required to handle alt_aln_method other than splign
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
	transcripts = hdp.get_tx_for_region(alt_ac,alt_aln_method,start_i,end_i)
	# Take the first transcript
#	gbtrs = re.compile('^NM_')
	for trans in transcripts:
#		if gbtrs.search(trans[0]):
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

# Maps genomic coordinates to coding if the c accession is known
def g_to_c(var_g, tx_ac, hp, evm):
	pat_g = re.compile("\:g\.") 		# Pattern looks for :g.
	# If the :g. pattern is present in the input variant
	if pat_g.search(var_g): 
		# convert the input string into a hgvs object by parsing
		var_g = hp.parse_hgvs_variant(var_g)
		# Map to coding variant
		var_c = str(evm.g_to_c(var_g, tx_ac))
		return var_c

# Maps genomic coordinates to coding if the c accession is known
def g_to_n(var_g, tx_ac, hp, evm):
	pat_g = re.compile("\:g\.") 		# Pattern looks for :g.
	# If the :g. pattern is present in the input variant
	if pat_g.search(var_g): 
		# convert the input string into a hgvs object by parsing
		var_g = hp.parse_hgvs_variant(var_g)
		# Map to coding variant
		var_n = str(evm.g_to_n(var_g, tx_ac))
		return var_n


# Return an hgvs object containing the coding sequence variant
def coding(variant, hp):
	# Set regular expressions for if statements
# 	pat_c = re.compile("\:c\.") 		# Pattern looks for :c. Note (gene) has been removed
	# If the :c. pattern is present in the input variant
	if re.search(':c.', variant) or re.search(':n.', variant): 
		# convert the input string into a hgvs object
		var_c = hp.parse_hgvs_variant(variant)
		return var_c
		

# Return an hgvs object containing the genomic sequence variant
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

# Return an hgvs object containing the genomic sequence variant
def protein(variant, evm, hp):
	variant = str(variant)
	# Set regular expressions for if statements
	pat_c = re.compile("\:c\.") 		# Pattern looks for :c. Note (gene) has been removed
	
	# If the :c. pattern is present in the input variant
	if  pat_c.search(variant): 
		# # # # # print variant		
		# convert the input string into a hgvs object
		var_c = hp.parse_hgvs_variant(variant)
		# # # # # print var_c
		# # # # # print var_c.ac		
		# map to the genomic sequence
		var_p = evm.c_to_p(var_c)	# genomic level variant
		# # # # print var_p		
		return var_p
	if re.search(':n.', variant):
		var_p = hp.parse_hgvs_variant(variant)
		var_p.ac = 'Non-coding transcript '
		var_p.posedit = ''
		# # # # # print var_p
		return var_p
	

# Return an hgvs object containing the rna sequence variant
def rna(variant, evm, hp):
	# Set regular expressions for if statements
	pat_c = re.compile("\:c\.") 		# Pattern looks for :c. Note (gene) has been removed
	# If the :c. pattern is present in the input variant
	if  pat_c.search(variant): 
		# convert the input string into a hgvs object
		var_c = hp.parse_hgvs_variant(variant)
		# map to the genomic sequence
		var_r = evm.c_to_n(var_c)	# rna level variant
		return var_r


def hgvs_rna(variant, hp):
	# Set regular expressions for if statements
	pat_r = re.compile("\:n\.") 		# Pattern looks for :n. Note (gene) has been removed
	# If the :r. pattern is present in the input variant
	if  pat_r.search(variant): 
		# convert the input string into a hgvs object
		var_r = hp.parse_hgvs_variant(variant)
		return var_r

def hgvs_genomic(variant, hp):
	# Set regular expressions for if statements
	pat_g = re.compile("\:g\.") 		# Pattern looks for :g. Note (gene) has been removed
	# If the :g. pattern is present in the input variant
	if  pat_g.search(variant): 
		# convert the input string into a hgvs object
		var_g = hp.parse_hgvs_variant(variant)
		return var_g

# Replacement for straightforward c_to_g takes into account mappings to multiple chromosomes
# def myevm_t_to_g(hgvs_c, evm, hdp, primary_assembly):
# 	Create normalizer
# 	hn = hgvs.normalizer.Normalizer(hdp,
# 			cross_boundaries=False,
# 			shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
# 			alt_aln_method='splign'
# 			)	
# 
# 	Create reverse normalizer
# 	reverse_normalize = hgvs.normalizer.Normalizer(hdp, 
# 			cross_boundaries=False, 
# 			shuffle_direction=5, 
# 			alt_aln_method='splign'
# 			)
# 
# 	create no_norm_evm
# 	no_norm_evm = hgvs.assemblymapper.AssemblyMapper(hdp,
# 			assembly_name=primary_assembly, 
# 			alt_aln_method='splign', 
# 			normalize=False, 
# 			replace_reference=True
# 			)
# 	Validator
# 	vr = hgvs.validator.Validator(hdp)			
# 	
# 	store the input
# 	stored_hgvs_c = copy.deepcopy(hgvs_c)
# 	
# 	Test mapping options
# 	mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)
# 	for line in mapping_options:
# 	# # print line
# 	
# 	# # print '\n\n\nTesting'
# 	# # print hgvs_c
# 	identity sub del delins ins dup inv
# 	Working
# 	Set expansion variable
# 	expand_out = 'false'
# 	if hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del' or hgvs_c.posedit.edit.type =='delins' or hgvs_c.posedit.edit.type == 'dup' or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins':
# 		
# 		if NM_ need the n. position
# 		if re.match('NM_', str(hgvs_c.ac)):
# 			hgvs_c = no_norm_evm.c_to_n(hgvs_c)
# 		Check for intronic
# 		try: 
# 			hn.normalize(hgvs_c)
# 		except hgvs.exceptions.HGVSError as e:
# 			error = str(e)
# 			if re.search('intronic variant', error):
# 				pass
# 		else:		
# 			For non-intronic sequence
# 			hgvs_t = copy.deepcopy(hgvs_c)
# 			if hgvs_c.posedit.edit.type == 'dup':
# 				hgvs_t = reverse_normalize.normalize(hgvs_t)
# 				pre_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
# 				post_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
# 				alt = pre_base + hgvs_t.posedit.edit.ref + hgvs_t.posedit.edit.ref + post_base				
# 				ref = pre_base + hgvs_t.posedit.edit.ref + post_base
# 				dup_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base - 1) + '_' + str((hgvs_t.posedit.pos.start.base + len(ref)) -2) + 'del' + ref + 'ins' + alt
# 				hgvs_t = hp.parse_hgvs_variant(dup_to_delins)
# 			elif hgvs_c.posedit.edit.type == 'ins':	
# 				pass
# 			else:	
# 				if str(hgvs_t.posedit.edit.alt) == 'None':
# 					hgvs_t.posedit.edit.alt = ''	
# 				pre_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
# 				post_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
# 				hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
# 				hgvs_t.posedit.edit.alt = pre_base + hgvs_t.posedit.edit.alt + post_base
# 				hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
# 				start = hgvs_t.posedit.pos.start.base
# 				hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
# 			
# 			if hgvs_t.posedit.edit.type =='delins':
# 				t = hgvs_t.ac + ':' + hgvs_c.type + '.' + str(start) + '_' + str(start + (len(hgvs_t.posedit.edit.ref) - 1)) + str(hgvs_t.posedit.edit)
# 				hgvs_c = hp.parse_hgvs_variant(t)
# 				# # print 'new ' + str(hgvs_c)
# 			else:
# 			hgvs_c = copy.deepcopy(hgvs_t)
# 			
# 			Set expanded out test to true
# 			expand_out = 'true'
# 			
# 		if re.match('NM_', str(hgvs_c.ac)):
# 			hgvs_c = no_norm_evm.n_to_c(hgvs_c)
# 		
# 		Ensure the altered c. variant has not crossed intro exon boundaries
# 		hgvs_check_boundaries = copy.deepcopy(hgvs_c)
# 		try:
# 			h_variant = hn.normalize(hgvs_check_boundaries)
# 		except hgvs.exceptions.HGVSError as e:
# 			error = str(e) 
# 			if re.search('spanning the exon-intron boundary', error):
# 				hgvs_c = copy.deepcopy(stored_hgvs_c)		 	
# 			
# 	# # print 'End test\n\n\n'
# 	try:
# 		hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
# 		This will fail on multiple refs for NC_
# 	except hgvs.exceptions.HGVSError as e:
#  		# # print '\n\n\n'
#  		# # print str(e)
#  		# # print '\n\n\n'
# 		Some transcripts map to ALT LOCI
# 		import dbControls
# 		import supported_chromosome_builds
# 		Make a variantmapper instance
# 		from hgvs import variantmapper
# 		from hgvs import validator
# 		vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=True) #, normalize=False)
# 		vr = hgvs.validator.Validator(hdp)
# 		get hgnc gene symbol assigned to transcript
# 		hgnc_symbol = hdp.get_tx_identity_info(hgvs_c.ac)[6]
# 		Search for current symbol
# 		current = hgnc_rest(path = "/search/prev_symbol/" + hgnc_symbol)
# 		if int(current['record']['response']['numFound']) == 0:
# 			pass
# 		else:
# 			hgnc_symbol = current['record']['response']['docs'][0]['symbol']
# 		Get chromosome location from database
# 		if primary_assembly == 'GRCh37':
# 			table = 'genePos37'
# 		if primary_assembly == 'GRCh38':
# 			table = 'genePos38'
# 		symbol = 'sym:' + hgnc_symbol
# 		gene_data = dbControls.data.in_entries(symbol, table)
# 		specified_chr = gene_data['chr']
# 		Get the alternative accession
# 		chr_accession = supported_chromosome_builds.to_accession(str(specified_chr), primary_assembly)
# 		if str(chr_accession) == 'None':
# 			Recover all available mapping options from UTA
# 			mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)	
# 			for option in mapping_options:
# 				if re.match('NT_', option[1]):
# 					try:
# 						hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
# 						break
# 					except:
# 						continue
# 			for option in mapping_options:
# 				if re.match('NW_', option[1]):
# 					try:
# 						hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
# 						break
# 					except:
# 						continue
# 		else:
# 			hgvs_genomic = vm.t_to_g(hgvs_c, chr_accession)		
# 	 
# 	if hgvs_genomic.posedit.edit.type == 'ins':
# 		try:
# 			hgvs_genomic = hn.normalize(hgvs_genomic)
# 		except hgvs.exceptions.HGVSError as e:
# 			error = str(e)
# 			if error == 'insertion length must be 1':
# 				Get orientation of the gene wrt genome and a list of exons mapped to the genome
# 				ori = tx_exons(tx_ac=hgvs_c.ac, alt_ac=hgvs_genomic.ac, alt_aln_method='splign', hdp=hdp)
# 				orientation = int(ori[0]['alt_strand'])
# 				if orientation	== 1:
# 					specified_base = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.start.base)
# 					post_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
# 					end_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base-1,hgvs_t.posedit.pos.end.base)
# 					start_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-1,hgvs_t.posedit.pos.start.base)
# 					pre_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
# 					if end_c_base == specified_base:	
# 						hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base - 1
# 						hgvs_genomic.posedit.edit.alt = start_c_base + hgvs_genomic.posedit.edit.alt
# 					if start_c_base == specified_base:
# 						hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base + 1
# 						hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt + end_c_base
# 				if orientation	== -1:
# 					specified_base = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.start.base)
# 					post_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
# 					end_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base-1,hgvs_t.posedit.pos.end.base)
# 					start_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-1,hgvs_t.posedit.pos.start.base)
# 					pre_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
# 					specified_base = revcomp(specified_base)
# 					post_c_base = revcomp(post_c_base)
# 					end_c_base = revcomp(end_c_base)
# 					start_c_base = revcomp(start_c_base)
# 					pre_c_base = revcomp(pre_c_base)				
# 					if end_c_base == specified_base:	
# 						hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base + 1
# 						hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt + start_c_base
# 					if start_c_base == specified_base:	
# 						hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base -1 
# 						hgvs_genomic.posedit.edit.alt = end_c_base + hgvs_genomic.posedit.edit.alt
# 				
# 	Remove identity bases
# 	elif expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) >= 3:
# 		hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
# 		hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
# 		hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
# 		hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1] 
# 	else:
# 		pass		
# 		
# 	# # print hgvs_genomic
# 	return hgvs_genomic

# Replacement for straightforward c_to_g takes into account mappings to multiple chromosomes
def myevm_t_to_g(hgvs_c, evm, hdp, primary_assembly):
	# Create normalizer
	hn = hgvs.normalizer.Normalizer(hdp,
			cross_boundaries=False,
			shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
			alt_aln_method='splign'
			)	

	# Create reverse normalizer
# 	reverse_normalize = hgvs.normalizer.Normalizer(hdp, 
# 			cross_boundaries=False, 
# 			shuffle_direction=5, 
# 			alt_aln_method='splign'
# 			)

	# create no_norm_evm
	no_norm_evm = hgvs.assemblymapper.AssemblyMapper(hdp,
			assembly_name=primary_assembly, 
			alt_aln_method='splign', 
			normalize=False, 
			replace_reference=True
			)
	# Validator
	vr = hgvs.validator.Validator(hdp)			
	
	# store the input
	stored_hgvs_c = copy.deepcopy(hgvs_c)
	
	# Test mapping options
	#mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)
	#for line in mapping_options:
	#	# # print line
	
	## # print '\n\n\nTesting'
	## # print hgvs_c
	# identity sub del delins ins dup inv
	# Working
	# Set expansion variable
	expand_out = 'false'
	if hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del' or hgvs_c.posedit.edit.type =='delins' or hgvs_c.posedit.edit.type == 'dup' or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins':
		
		# print 'in to function'
		# print hgvs_c
		
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
		try: 
			hn.normalize(hgvs_c)
		except hgvs.exceptions.HGVSError as e:
			error = str(e)
			if re.search('intronic variant', error):
				pass		

		else:		
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
			
		if re.match('NM_', str(hgvs_c.ac)):
			hgvs_c = no_norm_evm.n_to_c(hgvs_c)
		
		# print 'after expansion'
		# print hgvs_c
		
		# Ensure the altered c. variant has not crossed intro exon boundaries
		hgvs_check_boundaries = copy.deepcopy(hgvs_c)
		try:
			h_variant = hn.normalize(hgvs_check_boundaries)
		except hgvs.exceptions.HGVSError as e:
			error = str(e) 
			if re.search('spanning the exon-intron boundary', error):
				hgvs_c = copy.deepcopy(stored_hgvs_c)
			else:
				print 'boundary scan other error'			 	
				print error
				hgvs_c = copy.deepcopy(stored_hgvs_c)
				
	# # # print 'End test\n\n\n'
	try:
		hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
		# This will fail on multiple refs for NC_
	except hgvs.exceptions.HGVSError as e:
 		# print '\n\n\n'
 		# print hgvs_c
 		# print str(e)
 		# print '\n\n\n'
		# Some transcripts map to ALT LOCI
		import dbControls
		import supported_chromosome_builds
		# Make a variantmapper instance
		from hgvs import variantmapper
		from hgvs import validator
		vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=True) #, normalize=False)
		vr = hgvs.validator.Validator(hdp)
# 		get hgnc gene symbol assigned to transcript
# 		hgnc_symbol = hdp.get_tx_identity_info(hgvs_c.ac)[6]
# 		Search for current symbol
# 		current = hgnc_rest(path = "/search/prev_symbol/" + hgnc_symbol)
# 		if int(current['record']['response']['numFound']) == 0:
# 			pass
# 		else:
# 			hgnc_symbol = current['record']['response']['docs'][0]['symbol']
# 		Get chromosome location from database
# 		if primary_assembly == 'GRCh37':
# 			table = 'genePos37'
# 		if primary_assembly == 'GRCh38':
# 			table = 'genePos38'
# 		symbol = 'sym:' + hgnc_symbol
# 		gene_data = dbControls.data.in_entries(symbol, table)
# 		specified_chr = gene_data['chr']
# 		Get the alternative accession
# 		chr_accession = supported_chromosome_builds.to_accession(str(specified_chr), primary_assembly)
# 		if str(chr_accession) == 'None':
			# Recover all available mapping options from UTA
		mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)	
		# print mapping_options
		for option in mapping_options:
			# # print mapping_options
			if re.match('NC_', option[1]):
				chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
				if chr_num != 'false':
					# # # print '\n\nfound you mr ' + primary_assembly + '!!!\n\n' 
					try:
						hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
						break
					except:
						continue

		try:
			hgvs_genomic
		except:	
			for option in mapping_options:
				if re.match('NC_', option[1]):
					try:
						hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
						break
					except:
						continue
			try:
				hgvs_genomic
			except:
				for option in mapping_options:
					if re.match('NT_', option[1]):
						try:
							hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
							break
						except:
							continue
				try:
					hgvs_genomic
				except:							
					for option in mapping_options:
						if re.match('NW_', option[1]):
							try:
								hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
								break
							except:
								continue
# 		else:
# 			hgvs_genomic = vm.t_to_g(hgvs_c, chr_accession)		
	 
	if hgvs_genomic.posedit.edit.type == 'ins':
		# # print 'INS'
		# # print hgvs_genomic
		try:
			hgvs_genomic = hn.normalize(hgvs_genomic)
		except hgvs.exceptions.HGVSError as e:
			error = str(e)
			if error == 'insertion length must be 1':
				ref = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.end.base)
				hgvs_genomic.posedit.edit.ref = ref
				hgvs_genomic.posedit.edit.alt = ref[0:1] + hgvs_genomic.posedit.edit.alt + ref[-1:]
				hgvs_genomic = hn.normalize(hgvs_genomic)
				# Get orientation of the gene wrt genome and a list of exons mapped to the genome
# 				ori = tx_exons(tx_ac=hgvs_c.ac, alt_ac=hgvs_genomic.ac, alt_aln_method='splign', hdp=hdp)
# 				orientation = int(ori[0]['alt_strand'])
# 				if orientation	== 1:
# 					# # print '1'
# 					specified_base = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.start.base)
# 					post_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
# 					end_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base-1,hgvs_t.posedit.pos.end.base)
# 					start_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-1,hgvs_t.posedit.pos.start.base)
# 					pre_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
# 					if end_c_base == specified_base:	
# 						# # print 'a'
# 						hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base - 1
# 						hgvs_genomic.posedit.edit.alt = start_c_base + hgvs_genomic.posedit.edit.alt
# 						# # print hgvs_genomic
# 					if start_c_base == specified_base:
# 						# # print 'b'
# 						hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base + 1
# 						hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt + end_c_base
# 						# # print hgvs_genomic
# 					else:
# 						# # print 'c'
# 				if orientation	== -1:
# 					# # print '-1'
# 					specified_base = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.start.base)
# 					post_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
# 					end_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base-1,hgvs_t.posedit.pos.end.base)
# 					start_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-1,hgvs_t.posedit.pos.start.base)
# 					pre_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
# 					# specified_base = revcomp(specified_base)
# 					post_c_base = revcomp(post_c_base)
# 					end_c_base = revcomp(end_c_base)
# 					start_c_base = revcomp(start_c_base)
# 					pre_c_base = revcomp(pre_c_base)				
# 					if end_c_base == specified_base:	
# 						hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base + 1
# 						hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt + start_c_base
# 					if start_c_base == specified_base:	
# 						hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base -1 
# 						hgvs_genomic.posedit.edit.alt = end_c_base + hgvs_genomic.posedit.edit.alt
				
	# Remove identity bases
	elif expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) >= 3:
		hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
		hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
		hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
		hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1] 
	else:
		pass		
		
	# print 'Mapped to genomic'
	# print hgvs_genomic
	return hgvs_genomic

# Replacement for straightforward c_to_g takes into account mappings to multiple chromosomes
def noreplace_myevm_t_to_g(hgvs_c, evm, hdp, primary_assembly):
	try:
		hgvs_genomic = evm.t_to_g(hgvs_c)
		# This will fail on multiple refs for NC_
	except hgvs.exceptions.HGVSError:
		import dbControls
		import supported_chromosome_builds
		# Make a variantmapper instance
		from hgvs import variantmapper
		from hgvs import validator
		vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=False)
		vr = hgvs.validator.Validator(hdp)
		# get hgnc gene symbol assigned to transcript
		hgnc_symbol = hdp.get_tx_identity_info(hgvs_c.ac)[6]
		# Search for current symbol
		current = hgnc_rest(path = "/search/prev_symbol/" + hgnc_symbol)
		if int(current['record']['response']['numFound']) == 0:
			pass
		else:
			hgnc_symbol = current['record']['response']['docs'][0]['symbol']
		# Get chromosome location from database
		if primary_assembly == 'GRCh37':
			table = 'genePos37'
		if primary_assembly == 'GRCh38':
			table = 'genePos38'
		symbol = 'sym:' + hgnc_symbol
		gene_data = dbControls.data.in_entries(symbol, table)
		specified_chr = gene_data['chr']
		# Get the alternative accession
		chr_accession = supported_chromosome_builds.to_accession(str(specified_chr), primary_assembly)
		hgvs_genomic = vm.t_to_g(hgvs_c, chr_accession)		
	return hgvs_genomic


# VM method
def myvm_t_to_g(hgvs_c, alt_chr, vm, hn, hdp, primary_assembly):

	# create no_norm_evm
	no_norm_evm = hgvs.assemblymapper.AssemblyMapper(hdp,
			assembly_name=primary_assembly, 
			alt_aln_method='splign', 
			normalize=False, 
			replace_reference=True
			)
	# Validator
	vr = hgvs.validator.Validator(hdp)			
	
	# store the input
	stored_hgvs_c = copy.deepcopy(hgvs_c)
	
	# Test mapping options
	#mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)
	#for line in mapping_options:
	#	# # print line
	
	## # print '\n\n\nTesting'
	## # print hgvs_c
	# identity sub del delins ins dup inv
	# Working
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
		try: 
			hn.normalize(hgvs_c)
		except hgvs.exceptions.HGVSError as e:
			error = str(e)
			if re.search('intronic variant', error):
				pass

		else:		
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
			
		if re.match('NM_', str(hgvs_c.ac)):
			hgvs_c = no_norm_evm.n_to_c(hgvs_c)
		
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

# 				Get orientation of the gene wrt genome and a list of exons mapped to the genome
# 				ori = tx_exons(tx_ac=hgvs_c.ac, alt_ac=hgvs_genomic.ac, alt_aln_method='splign', hdp=hdp)
# 				orientation = int(ori[0]['alt_strand'])
# 				if orientation	== 1:
# 					specified_base = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.start.base)
# 					post_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
# 					end_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base-1,hgvs_t.posedit.pos.end.base)
# 					start_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-1,hgvs_t.posedit.pos.start.base)
# 					pre_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
# 					if end_c_base == specified_base:	
# 						hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base - 1
# 						hgvs_genomic.posedit.edit.alt = start_c_base + hgvs_genomic.posedit.edit.alt
# 					if start_c_base == specified_base:
# 						hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base + 1
# 						hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt + end_c_base
# 				if orientation	== -1:
# 					specified_base = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.start.base)
# 					post_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
# 					end_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base-1,hgvs_t.posedit.pos.end.base)
# 					start_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-1,hgvs_t.posedit.pos.start.base)
# 					pre_c_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
# 					specified_base = revcomp(specified_base)
# 					post_c_base = revcomp(post_c_base)
# 					end_c_base = revcomp(end_c_base)
# 					start_c_base = revcomp(start_c_base)
# 					pre_c_base = revcomp(pre_c_base)				
# 					if end_c_base == specified_base:	
# 						hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base + 1
# 						hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt + start_c_base
# 					if start_c_base == specified_base:	
# 						hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base -1 
# 						hgvs_genomic.posedit.edit.alt = end_c_base + hgvs_genomic.posedit.edit.alt

	# Remove identity bases
	elif expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) >= 3:
		hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
		hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
		hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
		hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1] 
	else:
		pass

	return hgvs_genomic

# Replacement for straightforward g_to_c/n takes into account mappings to multiple chromosomes
def myevm_g_to_t(hdp, evm, hgvs_genomic, alt_ac):
	from hgvs import validator
	vr = hgvs.validator.Validator(hdp)
	hgvs_t = evm.g_to_t(hgvs_genomic, alt_ac)
	#try:
	#	vr.validate(hgvs_t)
	#except (hgvs.exceptions.HGVSInvalidVariantError, hgvs.exceptions.HGVSUnsupportedOperationError) as e:
	#	error = str(e)		
	#	if re.search('does not agree with reference sequence', error):
	#		hgvs_t = ref_replace(e, hgvs_t)
	#	else:
	#		pass					
	return hgvs_t
	
def hgvs_protein(variant, hp):
	# Set regular expressions for if statements
	pat_p = re.compile("\:p\.") 		# Pattern looks for :g. Note (gene) has been removed
	# If the :p. pattern is present in the input variant
	if  pat_p.search(variant): 
		# convert the input string into a hgvs object
		var_p = hp.parse_hgvs_variant(variant)
		return var_p

# hgvs_r_to_c
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

# hgvs_c_to_r
def hgvs_c_to_r(hgvs_object):
	hgvs_object.type = 'r'
	edit = str(hgvs_object.posedit.edit)
	edit = edit.lower()
	edit = edit.replace('t', 'u')
	hgvs_object.posedit.edit = edit
	return hgvs_object	


# Return the identity information for the transcript variant (see uta.py)
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
		
# Alternative to the above but accepts the accession directly. Try to incorporate
def tx_id_info(alt_ac, hdp):
	tx_id_info = hdp.get_tx_identity_info(alt_ac)
	# NOTE The hgnc id is the 6th element in this list tx_ac is the 0th element in the list
	return tx_id_info

	

# Return tx information for a named hgnc gene (see uta.py)
def tx_for_gene(hgnc, hdp):
	# Interface with the UTA database via get_tx_for_gene in uta.py
	tx_for_gene = hdp.get_tx_for_gene(hgnc)
	return tx_for_gene 
	

# Extract Genomic refseq ID from tx_for_gene dictionary 
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
	
# Get genomic co-ordinates for variant start position
def int_start(var_g):
	start = var_g.posedit.pos.start
	# Stringify to get start co-ords
	start = str(start)
	# Make into an integer
	int_start = int(start)
	return int_start

# Get genomic co-ordinates for variant end position
def int_end(var_g):
	end = var_g.posedit.pos.end
	# Stringify to get start co-ords
	end = str(end)
	# Make into an integer
	int_end = int(end)
	return int_end
	

# Returns exon table of tx_acession exons aligned to the genomic refseq ID
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
	


# Compile a table containing exon data coding accession aligned to genomic refseq ID
def exon_table(tx_exons, gen_int_start, gen_int_end, cd_offset):
	
	# Set exon and intron counter
	ex_num = 0
	intron_num = -1
	
	# Set regular expressions for if statements
	pat_eq = re.compile("\=")			# Pattern looks for =
	
	# Set the table line variables
	start_exon_info = {}
	end_exon_info = {}

	# Set empty dictionary to compile table
	exon_table = {
		'start_info' : {},  
		'end_info' : {}
	}
	
	
	# LOOP THROUGH TX EXONS 
	for exon in tx_exons:
		# Add 1 to the exon number
		ex_num = ex_num + 1
		# Add 1 to intron_num
		intron_num = intron_num + 1 # Forst loop = 0

		# Deal with zero based offsets regarding exon 1 which includes UTR in the transcript
		if ex_num == 1:
			# INTERVALS IN EXONS - all controlled by the tx_exons table
			# Interval start
			# Extract the information about the variant start position wrt exon affected
			if exon['alt_strand'] == -1:
				g_ori = 'antisense'
		
				if (gen_int_start > exon['alt_start_i']) and (gen_int_start <= exon['alt_end_i']):
					# Remove the eq sign from the end of 'cigar'
					cigar = (pat_eq.sub('', exon['cigar']))
			
					exon_table['start_info'] = {
						'id' : 'exon ' + str(ex_num), 
						'c_start' : str(exon['tx_start_i'] + 1), 
						'c_end' : str(exon['tx_end_i']),
						'tr_start' : str(exon['tx_start_i'] + cd_offset), 
						'tr_end' : str(exon['tx_end_i'] + cd_offset), 
						'g_end' : str(exon['alt_start_i']), 
						'g_start' : str(exon['alt_end_i']), 
						'cigar' : cigar,
						'g_ori' : g_ori
					}
				
				# Interval end
				# Extract the information about the variant end position wrt exon affected
				if (gen_int_end > exon['alt_start_i']) and (gen_int_end <= exon['alt_end_i']):
					# Remove the eq sign from the end of 'cigar'
					cigar = (pat_eq.sub('', exon['cigar']))
			
					exon_table['end_info'] = {
						'id' : 'exon ' + str(ex_num), 
						'c_start' : str(exon['tx_start_i'] + 1), 
						'c_end' : str(exon['tx_end_i']),
						'tr_start' : str(exon['tx_start_i'] + cd_offset), 
						'tr_end' : str(exon['tx_end_i'] + cd_offset), 
						'g_end' : str(exon['alt_start_i']), 
						'g_start' : str(exon['alt_end_i']), 
						'cigar' : cigar,
						'g_ori' : g_ori
					}
			
			else:
				g_ori = 'sense'
				if (gen_int_start > exon['alt_start_i']) and (gen_int_start <= exon['alt_end_i']):
					# Remove the eq sign from the end of 'cigar'
					cigar = (pat_eq.sub('', exon['cigar']))
			
					exon_table['start_info'] = {
						'id' : 'exon ' + str(ex_num), 
						'c_start' : str(exon['tx_start_i'] + 1), 
						'c_end' : str(exon['tx_end_i']),
						'tr_start' : str(exon['tx_start_i'] + cd_offset), 
						'tr_end' : str(exon['tx_end_i'] + cd_offset), 
						'g_start' : str(exon['alt_start_i']), 
						'g_end' : str(exon['alt_end_i']), 
						'cigar' : cigar,
						'g_ori' : g_ori	
					}
				
				# Interval end
				# Extract the information about the variant end position wrt exon affected
				if (gen_int_end > exon['alt_start_i']) and (gen_int_end <= exon['alt_end_i']):
					# Remove the eq sign from the end of 'cigar'
					cigar = (pat_eq.sub('', exon['cigar']))
			
					exon_table['end_info'] = {
						'id' : 'exon ' + str(ex_num), 
						'c_start' : str(exon['tx_start_i'] + 1), 
						'c_end' : str(exon['tx_end_i']),
						'tr_start' : str(exon['tx_start_i'] + cd_offset), 
						'tr_end' : str(exon['tx_end_i'] + cd_offset),
						'g_start' : str(exon['alt_start_i']), 
						'g_end' : str(exon['alt_end_i']), 
						'cigar' : cigar,
						'g_ori' : g_ori
					}
		
		else:
			# INTERVALS IN EXONS - all controlled by the tx_exons table
			# Interval start
			# Extract the information about the variant start position wrt exon affected
			if exon['alt_strand'] == -1:
				g_ori = 'antisense'
		
				if (gen_int_start > exon['alt_start_i']) and (gen_int_start <= exon['alt_end_i']):
					# Remove the eq sign from the end of 'cigar'
					cigar = (pat_eq.sub('', exon['cigar']))
			
					exon_table['start_info'] = {
						'id' : 'exon ' + str(ex_num), 
						'c_start' : str(exon['tx_start_i'] + 1), 
						'c_end' : str(exon['tx_end_i']),
						# Here
						'tr_start' : str(exon['tx_start_i'] + cd_offset + 1), 
						'tr_end' : str(exon['tx_end_i'] + cd_offset), 
						'g_end' : str(exon['alt_start_i']), 
						'g_start' : str(exon['alt_end_i']), 
						'cigar' : cigar,
						'g_ori' : g_ori
					}
				
				# Interval end
				# Extract the information about the variant end position wrt exon affected
				if (gen_int_end > exon['alt_start_i']) and (gen_int_end <= exon['alt_end_i']):
					# Remove the eq sign from the end of 'cigar'
					cigar = (pat_eq.sub('', exon['cigar']))
			
					exon_table['end_info'] = {
						'id' : 'exon ' + str(ex_num), 
						'c_start' : str(exon['tx_start_i'] + 1), 
						'c_end' : str(exon['tx_end_i']),
						'tr_start' : str(exon['tx_start_i'] + cd_offset + 1), 
						'tr_end' : str(exon['tx_end_i'] + cd_offset), 
						'g_end' : str(exon['alt_start_i']), 
						'g_start' : str(exon['alt_end_i']), 
						'cigar' : cigar,
						'g_ori' : g_ori
					}
			
			else:
				g_ori = 'sense'
				if (gen_int_start > exon['alt_start_i']) and (gen_int_start <= exon['alt_end_i']):
					# Remove the eq sign from the end of 'cigar'
					cigar = (pat_eq.sub('', exon['cigar']))
			
					exon_table['start_info'] = {
						'id' : 'exon ' + str(ex_num), 
						'c_start' : str(exon['tx_start_i'] + 1), 
						'c_end' : str(exon['tx_end_i']),
						'tr_start' : str(exon['tx_start_i'] + cd_offset + 1), 
						'tr_end' : str(exon['tx_end_i'] + cd_offset), 
						'g_start' : str(exon['alt_start_i']), 
						'g_end' : str(exon['alt_end_i']), 
						'cigar' : cigar,
						'g_ori' : g_ori	
					}
				
				# Interval end
				# Extract the information about the variant end position wrt exon affected
				if (gen_int_end > exon['alt_start_i']) and (gen_int_end <= exon['alt_end_i']):
					# Remove the eq sign from the end of 'cigar'
					cigar = (pat_eq.sub('', exon['cigar']))
			
					exon_table['end_info'] = {
						'id' : 'exon ' + str(ex_num), 
						'c_start' : str(exon['tx_start_i'] + 1), 
						'c_end' : str(exon['tx_end_i']),
						'tr_start' : str(exon['tx_start_i'] + cd_offset + 1), 
						'tr_end' : str(exon['tx_end_i'] + cd_offset),
						'g_start' : str(exon['alt_start_i']), 
						'g_end' : str(exon['alt_end_i']), 
						'cigar' : cigar,
						'g_ori' : g_ori
					}
				
		# INTERVALS IN INTRONS - Take care of if statements controlling strandedness
		if exon['alt_strand'] == -1:
			g_ori = 'antisense'
			
			# If the interval start is < the start position of the current exon but > than the end position 
			# of the next exon in the table tx_exons
			if (gen_int_start <= exon['alt_start_i']) and (gen_int_start > tx_exons[intron_num + 1]['alt_end_i']):
				cigar = exon['alt_start_i'] - tx_exons[intron_num + 1]['alt_end_i']
			
				exon_table['start_info'] = {
					'id' : 'intron ' + str(intron_num + 1), 
					'c_start' : 'na', 
					'c_end' : 'na',
					'tr_start' : 'na',
					'tr_end' : 'na',
					'g_start' : str(exon['alt_start_i']),
					'g_end' : str(tx_exons[intron_num + 1]['alt_end_i']),
					'cigar' : cigar,
					'g_ori' : g_ori	
				}
			
			# If the interval end is < the start position of the current exon but > than the end position 
			# of the next exon in the table tx_exons
			if (gen_int_end <= exon['alt_start_i']) and (gen_int_end > tx_exons[intron_num + 1]['alt_end_i']):
				cigar = exon['alt_start_i'] - tx_exons[intron_num + 1]['alt_end_i']
			
				exon_table['end_info'] = {
					'id' : 'intron ' + str(intron_num + 1), 
					'c_start' : 'na', 
					'c_end' : 'na',
					'tr_start' : 'na',
					'tr_end' : 'na',
					'g_start' : str(exon['alt_start_i']),
					'g_end' : str(tx_exons[intron_num + 1]['alt_end_i']),
					'cigar' : cigar,
					'g_ori' : g_ori	
				}
		
		else:
			g_ori = 'sense'
			
			# If the interval start is < the start position of the current exon but > than the end position 
			# of the next exon in the table tx_exons
			if (gen_int_start > exon['alt_end_i']) and (gen_int_start < tx_exons[intron_num + 1]['alt_end_i']):
				cigar = tx_exons[intron_num + 1]['alt_start_i'] - exon['alt_end_i']
			
				exon_table['start_info'] = {
					'id' : 'intron ' + str(intron_num + 1), 
					'c_start' : 'na', 
					'c_end' : 'na',
					'tr_end' : 'na', 
					'tr_start' : 'na', 
					'g_end' : str(exon['alt_start_i']),
					'g_start' : str(tx_exons[intron_num - 1]['alt_end_i']),
					'cigar' : cigar,
					'g_ori' : g_ori	
				}
			
			# If the interval end is < the start position of the current exon but > than the end position 
			# of the next exon in the table tx_exons
			if (gen_int_end > exon['alt_end_i']) and (gen_int_end < tx_exons[intron_num + 1]['alt_end_i']):
				cigar = tx_exons[intron_num + 1]['alt_start_i'] - exon['alt_end_i']
			
				exon_table['end_info'] = {
					'id' : 'intron ' + str(intron_num + 1), 
					'c_start' : 'na', 
					'c_end' : 'na',
					'tr_end' : 'na', 
					'tr_start' : 'na', 
					'g_end' : str(exon['alt_start_i']),
					'g_start' : str(tx_exons[intron_num - 1]['alt_end_i']),
					'cigar' : cigar,
					'g_ori' : g_ori	
				}
	
	# END FOR LOOP	
	
	return exon_table

# Return relevant transcripts 
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
			
			# # # # # print coding	
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
				# # # # # print coding
		code_var.append(str(variant))
	return code_var


# Validate the						
def validate(input, hp, vr):
	hgvs_input = hp.parse_hgvs_variant(input)
	g = re.compile(":g.")
	p = re.compile(":p.")
	#if g.search(input):
	#	if hasattr(hgvs_input.posedit.pos.start, 'offset'):
	#		pass
	#	else:
	#		hgvs_input.posedit.pos.start.offset = 0
	#	if hasattr(hgvs_input.posedit.pos.end, 'offset'):
	#		pass
	#	else:
	#		hgvs_input.posedit.pos.end.offset = 0
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

# Extract accession sequences
def sequence_extractor(ac, hdp):
	ac_seq = hdp.get_tx_seq(ac)
	
	return ac_seq

# Update reference. Should only be used when evm has supplied the variant 
def ref_replace(e, hgvs_variant):
	error = str(e)
	match = re.findall('\(([GATC]+)\)', error)
	new_ref = match[1]
	#if re.search('=', str(hgvs_variant.posedit.edit)):
	#	hgvs_variant.posedit.edit.ref = new_ref
	#	hgvs_variant.posedit.edit.alt = new_ref
	#else: 
	hgvs_variant.posedit.edit.ref = new_ref
	return hgvs_variant	
	
# Search the ensembl database with rest
def ensembl_rest(ext, primary_assembly):
	# hash for decoded
	decoded = {
		'record' : '',
		'error' : 'false'
		}
	# Current fix in place
	primary_assembly = 'GRCh37'
	# ENSEMBLREST SERVER
	enr = "http://" + primary_assembly + ".rest.ensembl.org"	
	record = requests.get(enr+ext, headers={ "Content-Type" : "application/json"})
	# Check that the response is ok
	if not record.ok:
  		record.raise_for_status()
  		#sys.exit()
		decoded['error'] = "Unable to contact the Ensemble database: Please try again later" 
	else:
		# Decode
		decoded['record'] = record.json()
	return decoded
	
# Search the hgnc database with rest	
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
	
	
# Search Entrez databases with efetch and SeqIO
def entrez_efetch(db, id, rettype, retmode):			
	# IMPORT Bio modules
	from Bio import Entrez
	Entrez.email = 'pjf9@le.ac.uk'
	from Bio import SeqIO
	handle = Entrez.efetch(db=db, id=id, rettype=rettype, retmode=retmode)
	# Get record
	record = SeqIO.read(handle, "gb")
	# Place into text
	# text = handle.read()
	handle.close()
	return record
	
	
# search Entrez databases with efetch and read
def entrez_read(db, id,retmode):	
	# IMPORT Bio modules
	from Bio import Entrez
	Entrez.email = 'pjf9@le.ac.uk'
	from Bio import SeqIO
	handle = Entrez.efetch(db=db, id=id, retmode=retmode)
	# Get record
	record = Entrez.read(handle)
	# Place into text
	# text = handle.read()
	handle.close()
	return record
	
				
# NG to NC mapping
def ng_to_nc(rsg_id, primary_assembly, gen_acc):				
	ngnc = {
		'rsg_start' : '',
		'rsg_end' : ''
		}
	# Set the alignment file to be used
	if primary_assembly == 'GRCh37':
		ng = re.compile(rsg_id)
		# with open(os.path.join(APP_DATA, 'gene_RefSeqGene'), 'r') as inf
		with open(os.path.join(FUNCTIONS_ROOT, 'GCF_000001405.25_refseqgene_alignments-1.gff3'), 'r') as alnmt:
			for line in alnmt:
				line = line.rstrip()
				cells = line.split('\t')
				if cells[0] == gen_acc and ng.search(cells[8]):
					ngnc['rsg_start'] = cells[3]
					ngnc['rsg_end'] = cells[4]						
					break
				else:
					pass
	if primary_assembly == 'GRCh38':
		with open(os.path.join(FUNCTIONS_ROOT, 'GCF_000001405.28_refseqgene_alignments.gff3'), 'r') as alnmt:
			for line in alnmt:
				line = line.rstrip()
				cells = line.split('\t')
				if cells[0] == gen_acc and ng.search(cells[8]):
					ngnc['rsg_start'] = cells[3]
					ngnc['rsg_end'] = cells[4]						
					break
				else:
					pass	
	alnmt.close()
	return ngnc


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


def vcf2hgvs_submit(APP_ROOT, args):
	# Submit query to vcf2hgvs.py
	import subprocess
	# Select the correct python environment
	dev_ = re.compile("dev_batchValidator")
	if dev_.search(APP_ROOT):
		# nohup python bg_process.py and arguements
		script = os.path.join(APP_ROOT, 'vcf2hgvs.py')
		subprocess.Popen(["/local/python/2.7.12/bin/python", str(script)] + args)
	else:
		# nohup python bg_process.py and arguements
		script = os.path.join(APP_ROOT, 'vcf2hgvs.py')
		subprocess.Popen(["/local/python/2.7.12/bin/python", str(script)] + args)
	return

def vcf2hgvs_v2_submit(APP_ROOT, args):
	# Submit query to vcf2hgvs.py
	import subprocess
	# Select the correct python environment
	dev_ = re.compile("dev_batchValidator")
	if dev_.search(APP_ROOT):
		# nohup python bg_process.py and arguements
		script = os.path.join(APP_ROOT, 'dev_vcf2hgvs.py')
		subprocess.Popen(["/local/python/2.7.12/bin/python", str(script)] + args)
	else:
		# nohup python bg_process.py and arguements
		script = os.path.join(APP_ROOT, 'dev_vcf2hgvs.py')
		subprocess.Popen(["/local/python/2.7.12/bin/python", str(script)] + args)
	return	

def batch_validator_submit(APP_ROOT, args):
	# Submit query to Validator.py
	import subprocess
	# Select the correct python environment
	dev_ = re.compile("dev_batchValidator")
	if dev_.search(APP_ROOT):
		# nohup python bg_process.py and arguements
		script = os.path.join(APP_ROOT, 'validator.py')
		subprocess.Popen(["/local/python/2.7.12/bin/python", str(script)] + args)
	else:
		# nohup python bg_process.py and arguements
		script = os.path.join(APP_ROOT, 'validator.py')
		subprocess.Popen(["/local/python/2.7.12/bin/python", str(script)] + args)
	return
	
	
	
	
	
