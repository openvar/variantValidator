# curl examples
# http://127.0.0.1:5000/GRCh37/'Chr16:2099572TC>T' -d "transcript_ids=The Giraffe Was Ere"
# http://127.0.0.1:5000/GRCh37/'Chr16:2099572TC>T'/'NM_000548.3'
# http://127.0.0.1:5000/GRCh37/'Chr16:2099572TC>T'/'NM_001318829.1'
# http://127.0.0.1:5000/GRCh37/'Chr16:2099572TC>T'/'all'

# IMPORT HGVS MODULES
import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.dataproviders.seqfetcher		
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.sequencevariant
import hgvs.validator
import hgvs.exceptions
import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.normalizer

# IMPORT PYTHON MODULES
import re
import time
from contextlib import closing
import copy
import os
import sys
from operator import itemgetter

# Import Biopython
from Bio.Seq import Seq

# Import variantanalyser and peripheral VV modules
import ref_seq_type
import variantanalyser
from variantanalyser import functions as va_func
from variantanalyser import dbControls as va_dbCrl
from variantanalyser import hgvs2vcf as va_H2V
from variantanalyser import batch as va_btch
from variantanalyser import liftover as va_lo
from variantanalyser import g_to_g as va_g2g
from variantanalyser import supported_chromosome_builds as va_scb 

# Set debug mode
VALIDATOR_DEBUG = os.environ.get('VALIDATOR_DEBUG')
if VALIDATOR_DEBUG is not None:
	# Logging
	import logging
	logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

# PRE COMPILE VARIABLES
hdp = hgvs.dataproviders.uta.connect(pooling=True)
# From the hgvs parser import, create an instance of hgvs.parser.Parser
hp = hgvs.parser.Parser() 			
# Validator
vr = hgvs.validator.Validator(hdp)
# Variant mapper
vm = hgvs.variantmapper.VariantMapper(hdp)
# Create seqfetcher object
sf = hgvs.dataproviders.seqfetcher.SeqFetcher()
# create normalizers
hn = hgvs.normalizer.Normalizer(hdp,
		cross_boundaries=False,
		shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
		alt_aln_method='splign'
		)
reverse_normalize = hgvs.normalizer.Normalizer(hdp, 
		cross_boundaries=False, 
		shuffle_direction=5, 
		alt_aln_method='splign'
		)


# Error types
class variantValidatorException(Exception):
	pass

# method for final validation and stringifying parsed hgvs variants prior to # # # printing/passing to html
def valstr(hgvs_variant):
	# import re
	if re.search('del', str(hgvs_variant.posedit.edit)) or re.search('dup', str(hgvs_variant.posedit.edit)):
		if len(hgvs_variant.posedit.edit.ref) <= 4:
			pass
		else:	
			hgvs_variant.posedit.edit.ref = ''
		hgvs_variant = str(hgvs_variant)
	else:
		hgvs_variant = str(hgvs_variant)
	
	return hgvs_variant

# Defining reference sequence type from accession
def ref_type_assign(accession):
	if re.match('NC_', accession) or re.match('NG_', accession) or re.match('NT_', accession) or re.match('NW_', accession):
		ref_type = ':g.'		
	elif re.match('NM_', accession):
		ref_type = ':c.'
	elif re.match('NR_', accession):
		ref_type = ':n.'
	elif re.match('NP_', accession):
		ref_type = ':p.'	
	return ref_type		

# Check environ variables
def locate_dbs():
	HGVS_SEQREPO_DIR = os.environ.get('HGVS_SEQREPO_DIR')
	UTA_DB_URL = os.environ.get('UTA_DB_URL')
	VALIDATOR_DB_URL = os.environ.get('VALIDATOR_DB_URL')
	
	locate = {
			'HGVS_SEQREPO_DIR' : HGVS_SEQREPO_DIR,
			'UTA_DB_URL' : UTA_DB_URL,
			'VALIDATOR_DB_URL' : VALIDATOR_DB_URL
			}
	
	return locate

# Validator code
def validator(batch_variant, selected_assembly, select_transcripts):	

	# Set pre defined variables
	gene_symbols = ''
	alt_aln_method = 'splign'
	# SeqFetcher
	sf = hgvs.dataproviders.seqfetcher.SeqFetcher()
	
	try:	
		# Validation
		############

		# Create a dictionary of transcript ID : ''
		if select_transcripts != 'all':
			select_transcripts_list = select_transcripts.split('|')
			select_transcripts_dict = {}
			select_transcripts_dict_plus_version = {}
			for id in select_transcripts_list:
				id = id.strip()
				if re.match('LRG', id):
					id = va_dbCrl.data.get_RefSeqTranscriptID_from_lrgTranscriptID(id)
					if id == 'none':
						continue
				select_transcripts_dict_plus_version[id] = ''
				id = id.split('.')[0]
				select_transcripts_dict[id] = ''	
		# Set up gene list dictionary
		input_genes = {}

		# Remove genes if transcripts selected
		if select_transcripts != 'all':
			gene_symbols = ''
		else:
			pass	

		if gene_symbols != '':
			# List the genes
			gene_symbols = gene_symbols.strip()
			genes_list = gene_symbols.split()
			# Extract the required data
			for symbol in genes_list:
				symbol = symbol.strip()
				if symbol == '':
					continue
				else:
					query = va_dbCrl.data.get_transcribed_span_for_gene(symbol, selected_assembly)
					input_genes[symbol] = query	
		else:
			pass

		# split the batch queries into a list 
		batch_queries = batch_variant.split()

		# Turn each variant into a dictionary. The dictionary will be compiled during validation
		batch_list = []
		for queries in batch_queries:
			queries = queries.strip()
			query = {'quibble' : queries, 'id' : queries, 'warnings' : '', 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : 'false', 'order' : 'false'}
			batch_list.append(query) 

		# Create List to carry batch data
		batch_out = []

		# Ensure batch_list is pulled into the function so that it can be appended to
		batch_list = batch_list

		# Enter the validation loop
		###########################
		# Allow order by input
		ordering = 0
		for validation in batch_list:
			# Re-set cautions and automaps
			caution = ''
			automap = ''
			
			# This will be used to order the final output
			if str(validation['order']) == 'false':
				ordering = ordering + 1
				validation['order'] = ordering
			else:
				pass
			# Bug catcher 
			try:
				# Note, ID is not touched. It is always the input variant description. Quibble will be altered but id will not if type = g.
				input = validation['quibble']
				# Set the primary_assembly
				if validation['primary_assembly'] == 'false':
					primary_assembly = selected_assembly
					validation['primary_assembly'] = primary_assembly
				else:
					primary_assembly = validation['primary_assembly']
		
				# Set variables that batch will not use but are required
				crossing = 'false'
				boundary = 'false'
		
				# VCF type 1
				if re.search('-\d+-[GATC]+-[GATC]+', input):
					pre_input = copy.deepcopy(input)
					vcf_elements = pre_input.split('-')
					input = '%s:%s%s>%s' %(vcf_elements[0],vcf_elements[1],vcf_elements[2],vcf_elements[3])
				elif re.search('-\d+-[GATC]+-', input):
					pre_input = copy.deepcopy(input)
					vcf_elements = pre_input.split('-')
					input = '%s:%s%s>%s' %(vcf_elements[0],vcf_elements[1],vcf_elements[2],'del')		
				elif re.search('-\d+--[GATC]+', input):	
					pre_input = copy.deepcopy(input)
					vcf_elements = pre_input.split('-')
					input = '%s:%s%s>%s' %(vcf_elements[0],vcf_elements[1], 'ins', vcf_elements[2])		
			
				# API type non-HGVS
				# Chr16:2099572TC>T
				if re.search('\w+\:', input) and not re.search('\w+\:[gcnmrp].', input):
					sf = hgvs.dataproviders.seqfetcher.SeqFetcher()
					try:
						pre_input = copy.deepcopy(input)
						input_list = input.split(':')
						pos_ref_alt = str(input_list[1]) 
						positionAndEdit = input_list[1]
						if not re.match('N[CGTWMRP]_', input) and not re.match('LRG_', input):
							chr_num = str(input_list[0])
							chr_num = chr_num.upper()
							chr_num = chr_num.strip()
							chr_num = chr_num.replace('CHR', '')
							accession = va_scb.to_accession(chr_num, primary_assembly)	
						else:
							accession = input_list[0]	
						if re.search('>', pre_input):
							if re.search('del', pre_input):
								pos = re.match('\d+', pos_ref_alt)
								position = pos.group(0)
								old_ref,old_alt = pos_ref_alt.split('>')
								old_ref = old_ref.replace(position, '')
								position = int(position) - 1
								required_base = sf.fetch_seq(accession, start_i=position-1, end_i=position)
								ref = required_base + old_ref
								alt = required_base
								positionAndEdit = str(position) + ref + '>' + alt
							elif re.search('ins', pre_input):
								pos = re.match('\d+', pos_ref_alt)
								position = pos.group(0)
								old_ref,old_alt = pos_ref_alt.split('>')
								old_ref = old_ref.replace(position, '')
								position = int(position) - 1
								required_base = sf.fetch_seq(accession, start_i=position-1, end_i=position)
								ref = required_base
								alt = required_base + old_alt
								positionAndEdit = str(position) + ref + '>' + alt
						# Assign reference sequence type
						ref_type = ref_seq_type.ref_type_assign(accession)
						if re.match('LRG_', accession):
							if ref_type == ':g.':
								accession = va_dbCrl.data.get_RefSeqGeneID_from_lrgID(accesion)
							else:
								accession = va_dbCrl.data.get_RefSeqTranscriptID_from_lrgTranscriptID(accession)
						else:
							accession = accession
						input = str(accession) + ref_type + str(positionAndEdit)
						automap = 'Non-HGVS compliant variant description ' + input + ' automapped to ' + input + ': '
						validation['warnings'] = validation['warnings'] + ': ' + automap
					except:
						pass	
				

				# Ambiguous chr reference 
				if re.search('\w+\:[gcnmrp].', input) and not re.match('N[CGTWMR]_', input):
					if not re.match('LRG_', input) and not re.match('ENS', input):
						try:
							pre_input = copy.deepcopy(input)
							input_list = input.split(':')
							query_a_symbol = input_list[0]
							is_it_a_gene = va_dbCrl.data.get_current_hgnc_symbol(query_a_symbol, primary_assembly)
							if is_it_a_gene == 'none':
								pos_ref_alt = str(input_list[1]) 
								positionAndEdit = input_list[1]
								chr_num = str(input_list[0])
								chr_num = chr_num.upper()
								chr_num = chr_num.strip()
								chr_num = chr_num.replace('CHR', '')
								accession = va_scb.to_accession(chr_num, primary_assembly)	
								input = str(accession) + ':' + str(positionAndEdit)
								caution = 'Non-HGVS compliant variant description ' + input + ' automapped to ' + input + ':'
								validation['warnings'] = validation['warnings'] + ': ' + caution							
							else:
								pass
						except:
							pass

				# GENE_SYMBOL:c. n. types
				if re.search('\w+\:[cn].', input):
					try:
						pre_input = copy.deepcopy(input)
						query_a_symbol = pre_input.split(':')[0]
						tx_edit = pre_input.split(':')[1]
						is_it_a_gene = va_dbCrl.data.get_current_hgnc_symbol(query_a_symbol, primary_assembly)
						if is_it_a_gene != 'none':
							uta_symbol = va_dbCrl.data.get_uta_symbol(is_it_a_gene)
							available_transcripts = hdp.get_tx_for_gene(uta_symbol)
							select_from_these_transcripts = {}
							for tx in available_transcripts:
								if re.match('NM_', tx[3]) or re.match('NR_', tx[3]):
									if tx[3] not in select_from_these_transcripts.keys():
										select_from_these_transcripts[tx[3]] = ''
									else:
										continue											
								else:
									continue
							select_from_these_transcripts = '|'.join(select_from_these_transcripts.keys())
							if select_transcripts != 'all':
								validation['write'] = 'false'
								for transcript in select_transcripts_dict_plus_version.keys():
									validation['warnings'] = 'HGVS variant nomenclature does not allow the use of a gene symbol (' + query_a_symbol + ') in place of a valid reference sequence'
									refreshed_description = transcript + ':' + tx_edit
									query = {'quibble' : refreshed_description, 'id' : validation['id'], 'warnings' : validation['warnings'], 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
									batch_list.append(query)
							else:
								validation['warnings'] = validation['warnings'] + ': ' + 'HGVS variant nomenclature does not allow the use of a gene symbol (' + query_a_symbol + ') in place of a valid reference sequence: Re-submit ' + input + ' and specify transcripts from the following - ' + select_from_these_transcripts
							continue
						else:
							pass	
					except:
						pass

				# NG_:c. or NC_:c.
				if re.search('\w+\:[cn]', input):
					try:
						if re.match('^NG_', input):
							refSeqGeneID = input.split(':')[0]
							tx_edit = input.split(':')[1]
							gene_symbol = va_dbCrl.data.get_gene_symbol_from_refSeqGeneID(refSeqGeneID)
							if gene_symbol != 'none':
								uta_symbol = va_dbCrl.data.get_uta_symbol(gene_symbol)
								available_transcripts = hdp.get_tx_for_gene(uta_symbol)
								select_from_these_transcripts = {}
								for tx in available_transcripts:
									if re.match('NM_', tx[3]) or re.match('NR_', tx[3]):
										if tx[3] not in select_from_these_transcripts.keys():
											select_from_these_transcripts[tx[3]] = ''
										else:
											continue											
									else:
										continue
								select_from_these_transcripts = '|'.join(select_from_these_transcripts.keys())
								if select_transcripts != 'all':
									validation['write'] = 'false'
									for transcript in select_transcripts_dict_plus_version.keys():
										validation['warnings'] = 'NG_:c.PositionVariation descriptions should not be used unless a transcript reference sequence has also been provided e.g. NG_(NM_):c.PositionVariation'										
										refreshed_description = refSeqGeneID + '(' + transcript + ')' + ':' +  tx_edit
										query = {'quibble' : refreshed_description, 'id' : validation['id'], 'warnings' : validation['warnings'], 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
										batch_list.append(query)
								else:
									validation['warnings'] = validation['warnings'] + ': ' + 'A transcript reference sequence has not been provided e.g. NG_(NM_):c.PositionVariation. Re-submit ' + input + ' but also specify transcripts from the following - ' + str(select_from_these_transcripts)
								continue									
							else:
								validation['warnings'] = validation['warnings'] + ': ' + 'A transcript reference sequence has not been provided e.g. NG_(NM_):c.PositionVariation'
							continue
						elif re.match('^NC_', input):
							validation['warnings'] = validation['warnings'] + ': ' + 'A transcript reference sequence has not been provided e.g. NC_(NM_):c.PositionVariation. Unable to predict available transripts because chromosomal position is not specified'
							continue
						else:
							pass
					except:
						pass						

				# Find not_sub type in input e.g. GGGG>G
				not_sub = copy.deepcopy(input)
				not_sub_find = re.compile("([GATCgatc]+)>([GATCgatc]+)")
				if not_sub_find.search(not_sub):
					try:
						# If the length of either side of the substitution delimer (>) is >1
						matches = not_sub_find.search(not_sub)
						if len(matches.group(1)) > 1 or len(matches.group(2)) > 1:
							# Search for and remove range
							range = re.compile("([0-9]+)_([0-9]+)")
							if range.search(not_sub):
								m = not_sub_find.search(not_sub)
								start = m.group(1)
								delete = m.group(2)
								beginning_string,middle_string = not_sub.split(':')
								middle_string = middle_string.split('_')[0]
								end_string = start + '>' + delete
								not_sub = beginning_string + ':' +  middle_string + end_string
							# Split description
							split_colon = not_sub.split(':')
							ref_ac = split_colon[0]
							remainder = split_colon[1]
							split_dot = remainder.split('.')
							ref_type = split_dot[0]
							remainder = split_dot[1]
							posedit = remainder
							split_greater = remainder.split('>')
							insert = split_greater[1]
							remainder = split_greater[0]
							# Split remainder using matches
							r = re.compile("([0-9]+)([GATCgatc]+)")
							try:
								m = r.search(remainder)
								start = m.group(1)
								delete = m.group(2)
								starts = posedit.split(delete)[0]
								re_try = ref_ac + ':' + ref_type + '.' +  starts + 'del' + delete[0] + 'ins' +  insert
								hgvs_re_try = hp.parse_hgvs_variant(re_try)
								hgvs_re_try.posedit.edit.ref = delete
								start_pos = str(hgvs_re_try.posedit.pos.start)
								if re.search('\-', start_pos):
									base,offset = start_pos.split('-')
									new_offset = 0-int(offset) + (len(delete))
									end_pos = int(base)
									hgvs_re_try.posedit.pos.end.base = int(end_pos)
									hgvs_re_try.posedit.pos.end.offset = int(new_offset)-1
									not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
								elif re.search('\+', start_pos):
									base,offset = start_pos.split('+')
									end_pos = int(base) + (len(delete) - int(offset) - 1)
									new_offset = 0+int(offset) + (len(delete) - 1)
									hgvs_re_try.posedit.pos.end.base = int(end_pos)
									hgvs_re_try.posedit.pos.end.offset = int(new_offset)
									not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(hgvs_re_try.posedit.pos.end) + 'del' + delete + 'ins' + insert
								else:
									end_pos = int(start_pos) + (len(delete) - 1)
									not_delins = ref_ac + ':' + ref_type + '.' + start_pos + '_' + str(end_pos) + 'del' + delete + 'ins' + insert	
							except:
								not_delins = not_sub
							# Parse into hgvs object
							try:
								hgvs_not_delins = hp.parse_hgvs_variant(not_delins)
							except hgvs.exceptions.HGVSError as e:
								error = str(e)
								issue_link = ''
								validation['warnings'] = validation['warnings'] + ': ' + error
								continue

							try:
								not_delins = str(hn.normalize(hgvs_not_delins))
							except hgvs.exceptions.HGVSError as e:
								error = str(e)
								if re.search('Normalization of intronic variants is not supported', error):
									not_delins = not_delins
								else:
									issue_link = ''
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									continue	
							# Create warning
							caution = 'Variant description ' + input + ' is not HGVS compliant'
							automap = input + ' automapped to ' + not_delins				
							validation['warnings'] = validation['warnings'] + ': ' + automap
							# Change input to normalized variant
							input = not_delins
						else:
							pass
					except:
						pass			
				else:
					pass

				# Tackle edit1234 type
				edit_pass = re.compile('_\d+$')
				edit_fail = re.compile('\d+$')
				if edit_fail.search(input):
					if edit_pass.search(input):
						pass
					else:
						error = 'false'
						issue_link = 'false'
						failed = copy.deepcopy(input)
						# Catch the trailing digits
						digits = re.search(r"(\d+$)", failed)
						digits = digits.group(1)
						# Remove them so that the string SHOULD parse
						try:
							hgvs_failed = hp.parse_hgvs_variant(failed)
						except hgvs.exceptions.HGVSError as e:
							error = str(e)
							error = 'The syntax of the input variant description is invalid '
							if re.search('ins\d+', failed): 
								issue_link = 'http://varnomen.hgvs.org/recommendations/DNA/variant/insertion/'
								error = error + ' please refer to ' + issue_link 
							validation['warnings'] = validation['warnings'] + error
							continue
						hgvs_failed = hp.parse_hgvs_variant(failed)
						hgvs_failed.posedit.edit = str(hgvs_failed.posedit.edit).replace(digits, '')
						failed = str(hgvs_failed)
						hgvs_failed = hp.parse_hgvs_variant(failed)
						automap = 'Non HGVS compliant variant description '+ input + ' automapped to ' + failed
						validation['warnings'] = validation['warnings'] + ': ' + automap
						input = failed 

		
				# Tackle compound variant descriptions NG or NC (NM_) i.e. correctly input NG/NC_(NM_):c. 
				caution = ''
				compounder = re.compile('\(NM_')
				compounder_b = re.compile('\(ENST')
				if compounder.search(input):
					# Find pattern e.g. +0000 and assign to a variable
					transy = re.search(r"(\NM_.+)", input)
					transy = transy.group(1)
					transy = transy.replace(')', '')
					input = transy
			
				# INITIAL USER INPUT FORMATTING
				# Removes whitespace from the ends of the string
				# Removes anything in brackets
				# Identifies variant type
				# Returns a dictionary containing the formated input string and the variant type
				# Accepts c, g, n, r currently
				formatted = va_func.user_input(input)

				# Validator specific variables, note, not all will be necessary for batch, but keep to ensure that batch works
				vars = []
				hgnc = ''
				refseq_gene = ''
				exon_table = ''
				relevant = ''
				warning = ''
				automap = 'false'
				vmapped = 'false'
				coords = 'false'
				ensembl_gene = 'false'
				hgnc_gene_info = 'false'
				issue_link = 'false'
				cr_available = 'false'
				rcmds_tab = 'false'
				hgnc = ''

				# Check the initial validity of the input
				if formatted == 'invalid':
					validation['warnings'] = validation['warnings'] + ': ' + 'Input variant descriptions is not HGVS compliant'
					continue
				else:
					variant = formatted['variant']
					input = formatted['variant']
					type = formatted['type']

				# Conversions
				conversion = re.compile('con')
				if conversion.search(variant):
					validation['warnings'] = validation['warnings'] + ': ' + 'Gene conversions currently unsupported'
					continue

				# Primary check that hgvs will accept the variant 
				error = 'false'
				# Change RNA bases to upper case but nothing else
				if type == ":r.":
					variant = variant.upper()
					variant = variant.replace(':R.', ':r.')
					# lowercase the supported variant types
					variant = variant.replace('DEL', 'del')
					variant = variant.replace('INS', 'ins')
					variant = variant.replace('INV', 'inv')
					variant = variant.replace('DUP', 'dup')

				try:
					input_parses = hp.parse_hgvs_variant(variant)
				except hgvs.exceptions.HGVSError as e:
					error = str(e)
				if error == 'false':
					pass
				else:
					validation['warnings'] = validation['warnings'] + ': ' + str(error)
					continue	

				if re.match('^ENST', str(input_parses)):
					trap_ens_in = str(input_parses)
					sim_tx = hdp.get_similar_transcripts(input_parses.ac)
					for line in sim_tx:
						if str(line[2]) == 'True' and str(line[3]) == 'True' and str(line[4]) == 'True' and str(line[5]) == 'True' and str(line[6]) == 'True':
							input_parses.ac = (line[1])
							input = str(input_parses)
							variant = input
							break
					if re.match('^ENST', str(input_parses)):
						reason = 'Unable to continue validation'
						error = 'Unable to map ' + str(input_parses.ac) + ' to an equivalent RefSeq transcript'  
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue			
					else:
						validation['warnings'] = validation['warnings'] + ': ' + str(trap_ens_in) + ' automapped to equivalent RefSeq transcript ' + variant 

				# Catch interval end > interval start
				astr = re.compile('\*')
				if astr.search(str(input_parses.posedit)):
					input_parses_copy = copy.deepcopy(input_parses)
					input_parses_copy.type = "c"
					# Map to n. position
					# Create easy variant mapper (over variant mapper) and splign locked evm
					evm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=primary_assembly, alt_aln_method='splign', normalize=True, replace_reference=True)
					try:
						to_n = evm.c_to_n(input_parses_copy)
					except hgvs.exceptions.HGVSError as e:
						error = str(e)
						reason = 'Invalid variant description'
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue
					if to_n.posedit.pos.end.base < to_n.posedit.pos.start.base:
						error = 'Interval end position < interval start position '
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue
				elif input_parses.posedit.pos.end.base < input_parses.posedit.pos.start.base:
					error = 'Interval end position '  + str(input_parses.posedit.pos.end.base) + ' < interval start position ' + str(input_parses.posedit.pos.start.base)
					validation['warnings'] = validation['warnings'] + ': ' + str(error)
					continue
				else:
					pass

				# Catch missing version number in refseq
				ref_type = re.compile("^N\w\w\d")
				is_version = re.compile("\d\.\d")
				en_type = re.compile('^ENS')
				lrg_type = re.compile('LRG')
				if (ref_type.search(str(input_parses)) and is_version.search(str(input_parses))) or (en_type.search(str(input_parses))):
					pass
				else:	
					if lrg_type.search(str(input_parses)):
						#return render_template('bootstrap/variantError.html', title=title, user=input, error='LRG reference sequences are not currently supported', reason = 'Invalid Variant description', issue_link=issue_link)
						pass
					if ref_type.search(str(input_parses)):
						error='RefSeq variant accession numbers MUST include a version number'
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue

				# Ensure ENS sequences are always genebuild
				ens = re.compile('^ENS')
				if ens.search(variant):
					if alt_aln_method != 'genebuild':
						caution = '"'+alt_aln_method+'"' + " alignments are not available for Ensembl sequences"
						alt_aln_method = 'genebuild'
						automap = 'Automap has set the alignment method to "genebuild"'
						validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(automap)
						# No additional information required
						vmapped = 'false'

				# Allows bridging between alignment types - ENS to refseq
				elif alt_aln_method == 'genebuild' and type != ':g.':
					error = 'genebuild alignments available for Ensembl (ENS) sequences only'
					validation['warnings'] = validation['warnings'] + ': ' + str(error)
					continue
				else:
					pass

				# Check whether supported genome build is requested for non g. descriptions
				historic_assembly = 'false'
				mapable_assemblies = {
					'GRCh37' : 'true',
					'GRCh38' : 'true',
					'NCBI36' : 'false'
					}	
				is_mapable = mapable_assemblies.get(primary_assembly)
				if is_mapable == 'true':
					# Create easy variant mapper (over variant mapper) and splign locked evm
					evm = hgvs.assemblymapper.AssemblyMapper(hdp, 
							assembly_name=primary_assembly, 
							alt_aln_method=alt_aln_method, 
							normalize=True, 
							replace_reference=True
							)
					# create normalizer
# 					hn = hgvs.normalizer.Normalizer(hdp,
# 							cross_boundaries=False,
# 							shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
# 							alt_aln_method=alt_aln_method
# 							)
					# Setup a reverse normalize instance and non-normalize evm
					no_norm_evm = hgvs.assemblymapper.AssemblyMapper(hdp,
							assembly_name=primary_assembly, 
							alt_aln_method=alt_aln_method, 
							normalize=False, 
							replace_reference=True
							)
# 					reverse_normalize = hgvs.normalizer.Normalizer(hdp, 
# 							cross_boundaries=False, 
# 							shuffle_direction=5, 
# 							alt_aln_method=alt_aln_method
# 							)							

				# Will extract any genomic variants and append to list		
				elif is_mapable == 'false' and type == ':g.':
					historic_assembly = copy.deepcopy(primary_assembly)
					primary_assembly = 'GRCh37'
					# Create easy variant mapper (over variant mapper) and splign locked evm
					evm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=primary_assembly, alt_aln_method=alt_aln_method, normalize=True, replace_reference=True)

					# Create normalizer instance for primary/which will be target_assembly
# 					hn = hgvs.normalizer.Normalizer(hdp,
# 							cross_boundaries=False,
# 							shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
# 							alt_aln_method=alt_aln_method
# 							)		
					# Liftover
					hgvs_variant = hp.parse_hgvs_variant(variant)
					target_hgvs_variant = va_lo.hgvs_liftover(hgvs_variant, primary_assembly=historic_assembly, target_assembly=primary_assembly, hn=hn)
					# Check liftover has completed
					if str(target_hgvs_variant) != str(hgvs_variant):
						try:
							vr.validate(target_hgvs_variant)
						except hgvs.exceptions.HGVSError as e:
							error = valstr(hgvs_variant) + ' (' + historic_assembly + ') lifted over to ' + valstr(target_hgvs_variant) + ' (' + primary_assembly + ')' 
							error = error + ': ' + str(e)
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
							continue
						else:
							# Set variant
							variant = str(target_hgvs_variant)
							# Tag relevant sections of output line
							validation['genomic_g'] = str(target_hgvs_variant)
							validation['primary_assembly'] = primary_assembly
							error = valstr(hgvs_variant) + ' (' + historic_assembly + ') lifted over to ' + valstr(target_hgvs_variant) + ' (' + primary_assembly + ')'
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
					else:
						error='Liftover from ' + historic_assembly + ' to ' + primary_assembly + ' failed.'
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue
				else:
					error='Mapping of ' + variant + ' to genome assembly ' + primary_assembly + ' is not supported'
					validation['warnings'] = validation['warnings'] + ': ' + str(error)
					continue			
		
				# handle LRG inputs
				if re.match('^LRG', str(input_parses)):
					if re.match('^LRG\d+', str(input_parses.ac)):
						string = str(input_parses.ac)
						reference = string.replace('LRG', 'LRG_')	
						input_parses.ac = reference
						caution = string + ' updated to ' + reference
					if not re.match('^LRG_\d+', str(input_parses)):
						pass
					elif re.match('^LRG_\d+:g.', str(input_parses)) or re.match('^LRG_\d+:p.', str(input_parses)) or re.match('^LRG_\d+:c.', str(input_parses)) or re.match('^LRG_\d+:n.', str(input_parses)):
						lrg_reference, variation = str(input_parses).split(':')
						refseqgene_reference = va_dbCrl.data.get_RefSeqGeneID_from_lrgID(lrg_reference)
						if refseqgene_reference != 'none':
							input_parses.ac = refseqgene_reference
							variant = str(input_parses)
							input = str(input_parses)
							if caution == '':
								caution = lrg_reference + ':' + variation + ' automapped to ' + refseqgene_reference + ':' + variation
							else:
								caution = caution + ': ' +	lrg_reference + ':' + variation + ' automapped to ' + refseqgene_reference + ':' + variation
							validation['warnings'] = validation['warnings'] + ': ' + str(caution)
					elif re.match('^LRG_\d+t\d+:c.', str(input_parses)) or re.match('^LRG_\d+t\d+:n.', str(input_parses)) or re.match('^LRG_\d+t\d+:p.', str(input_parses)) or  re.match('^LRG_\d+t\d+:g.', str(input_parses)):		
						lrg_reference, variation = str(input_parses).split(':')
						refseqtranscript_reference = va_dbCrl.data.get_RefSeqTranscriptID_from_lrgTranscriptID(lrg_reference)
						if refseqtranscript_reference != 'none':
							input_parses.ac = refseqtranscript_reference
							variant = str(input_parses)
							input = str(input_parses)
							if caution == '':
								caution = lrg_reference + ':' + variation + ' automapped to ' + refseqtranscript_reference + ':' + variation
							else:
								caution = caution + ': ' + lrg_reference + ':' + variation + ' automapped to ' + refseqtranscript_reference + ':' + variation
							validation['warnings'] = validation['warnings'] + ': ' + str(caution)
					else:
						pass

				# Additional Incorrectly input variant capture training
				# NM_ .g
				if (re.search('^NM_', variant) or re.search('^NR_', variant)) and re.search(':g.', variant):
					suggestion = input.replace(':g.', ':c.')
					error = 'Transcript reference sequence input as genomic (g.) reference sequence. Did you mean ' + suggestion + '?'
					validation['warnings'] = validation['warnings'] + ': ' + error
					continue
				# NR_ c.
				if re.search('^NR_', input) and re.search(':c.', input):
					suggestion = input.replace(':c.', ':n.')
					error = 'Non-coding transcript reference sequence input as coding (c.) reference sequence. Did you mean ' + suggestion + '?'
					validation['warnings'] = validation['warnings'] + ': ' + error
					continue
				# NM_ n.
				if re.search('^NM_', input) and re.search(':n.', input):
					suggestion = input.replace(':n.', ':c.')
					error = 'Coding transcript reference sequence input as non-coding transcript (n.) reference sequence. Did you mean ' + suggestion + '?'
					validation['warnings'] = validation['warnings'] + ': ' + error
					continue

				# NM_ NC_ NG_ NR_ p.
				if (re.search('^NM_', variant) or re.search('^NR_', variant) or re.search('^NC_', variant) or re.search('^NG_', variant)) and re.search(':p.', variant):
					issue_link = 'http://varnomen.hgvs.org/recommendations/protein/'
					error = 'Coding transcript reference sequence input as non-coding transcript (n.) reference sequence. Did you mean ' + suggestion + '?'
					validation['warnings'] = validation['warnings'] + ': ' + error
					continue
					
				# NG_ c or NC_c..
				if (re.search('^NG_', variant) or re.search('^NC_', variant)) and re.search(':c.', variant):
					suggestion = ': For additional assistance, submit ' + str(variant) + ' to VariantValidator'
					error = 'NG_:c.PositionVariation descriptions should not be used unless a transcript reference sequence has also been provided e.g. NG_(NM_):c.PositionVariation' + suggestion
					validation['warnings'] = validation['warnings'] + ': ' + error
					continue
					
				# Primary validation of the input
				# re-make input_parses
				input_parses = hp.parse_hgvs_variant(input)
				if input_parses.type == 'g':
					if re.match('^NC_', input_parses.ac) or re.match('^NG_', input_parses.ac) or re.match('^NT_', input_parses.ac) or re.match('^NW_', input_parses.ac):
						pass
					else:
						error = 'Invalid reference sequence identifier (' + input_parses.ac + ')'
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue
					try: 
						vr.validate(input_parses)
					except hgvs.exceptions.HGVSError as e:
						error = str(e)
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue
					else:
						pass
						
				elif input_parses.type == 'c':  	
					if re.search('\*', str(input_parses)) or re.search('c.\-', str(input_parses)):
						# Catch variation in UTRs
						# These should be in the sequence so can be directly validated. Need to pass to n.
						try: 
							vr.validate(input_parses)
						except hgvs.exceptions.HGVSError as e:
							error = str(e)
							if re.search('intronic variant', error):
								pass
							elif re.search('datums is ill-defined', error):
								called_ref = input_parses.posedit.edit.ref
								to_n = evm.c_to_n(input_parses)
								actual_ref = to_n.posedit.edit.ref
								if called_ref != actual_ref:
									error = 'Variant reference (' + called_ref + ') does not agree with reference sequence (' + actual_ref +')'
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									continue
								else:					
									input_parses.posedit.edit.ref = ''
									variant = str(input_parses)	
							else:
								reason = 'Invalid variant description'
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue								
						try:
							input_parses = evm.c_to_n(input_parses)
						except hgvs.exceptions.HGVSError as e:
							error = str(e)
							reason = 'Invalid variant description'
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
							continue 

						if re.search('n.1-', str(input_parses)):
							input_parses = evm.n_to_c(input_parses)
							error = 'Using a transcript reference sequence to specify an intergenic variant position that lies 5' + "'" + ' to the transcript reference sequence is not HGVS compliant. Instead use '
							genomic_position = va_func.myevm_t_to_g(input_parses, evm, hdp, primary_assembly)
							error = error + str(genomic_position)
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
							continue
						else:
							pass

						# Re-map input_parses back to c. variant					
						input_parses = evm.n_to_c(input_parses)			

					elif re.search('\d\-', str(input_parses)) or re.search('\d\+', str(input_parses)):  
						# Create a specific minimal evm with no normalizer and no replace_reference
						min_evm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=primary_assembly, alt_aln_method=alt_aln_method, normalize=False, replace_reference=False)
						# Have to use this method due to potential multi chromosome error, note, normalizes but does not replace sequence
						try:
							output = va_func.noreplace_myevm_t_to_g(input_parses, min_evm, hdp, primary_assembly)
						except hgvs.exceptions.HGVSDataNotAvailableError as e:
							tx_ac = input_parses.ac
							try:
								gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(tx_ac)
							except:
								gene_symbol = 'none'	
							if gene_symbol == 'none':
								error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting '  + tx_ac + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
							else:
								error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
							continue						
						except ValueError as e:
							error = str(e)
							if re.search('> end', error):
								error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end position ' + str(input_parses.posedit.pos.end)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
						except hgvs.exceptions.HGVSInvalidVariantError as e:
							error = str(e)
							if re.search('base start position must be <= end position', error):
								correction = copy.deepcopy(input_parses)
								st = input_parses.posedit.pos.start
								ed = input_parses.posedit.pos.end
								correction.posedit.pos.start = ed
								correction.posedit.pos.end = st
								error = error + ': Did you mean ' + str(correction) + '?' 
								error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end position ' + str(input_parses.posedit.pos.end)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
						try: 
							vr.validate(output)
						except hgvs.exceptions.HGVSError as e:
							error = str(e)
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
							continue				
				
					else:
						# All other variation 
						try:
							vr.validate(input_parses)
						except hgvs.exceptions.HGVSUnsupportedOperationError:
							pass
						except hgvs.exceptions.HGVSInvalidVariantError as e:
							error = str(e)
							if re.search('Length implied by coordinates', error):
								# Applies to del and inv
								# NOTE, there has been no normalization at all so this error is valid here
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue				 
							# Will apply to > del and inv
							if re.search('does not agree with reference sequence', error):
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue						
								# ensures x_y for insertions
							if re.search('insertion length must be 1', error):
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue					
								# This catches errors in introns
							if re.search('base start position must be <= end position', error):
								correction = copy.deepcopy(input_parses)
								st = input_parses.posedit.pos.start
								ed = input_parses.posedit.pos.end
								correction.posedit.pos.start = ed
								correction.posedit.pos.end = st
								error = error + ': Did you mean ' + str(correction) + '?' 
								error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end position ' + str(input_parses.posedit.pos.end)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue

						except hgvs.exceptions.HGVSDataNotAvailableError as e:
							error = e
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
							continue
						except hgvs.exceptions.HGVSError as e:
							error = str(e)			
							if re.search('outside the bounds', error):
								error = error + ' (' + input_parses.ac + ')'
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
								

				elif input_parses.type == 'n':
					if re.search('\+', str(input_parses)) or re.search('\-', str(input_parses)):
						# Catch variation in UTRs
						# These should be in the sequence so can be directly validated. Need to pass to n.
						try: 
							vr.validate(input_parses)
						except hgvs.exceptions.HGVSError as e:
							error = str(e)
							if re.search('intronic variant', error):
								pass
							elif re.search('datums is ill-defined', error):
								called_ref = input_parses.posedit.edit.ref
								to_n = evm.c_to_n(input_parses)
								actual_ref = to_n.posedit.edit.ref
								if called_ref != actual_ref:
									error = 'Variant reference (' + called_ref + ') does not agree with reference sequence (' + actual_ref +')'
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									continue
								else:					
									input_parses.posedit.edit.ref = ''
									variant = str(input_parses)

							elif re.search('base must be >=1 for datum = SEQ_START or CDS_END', error):
									error = 'The given coordinate is outside the bounds of the reference sequence.'
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									continue									
							else:
								reason = 'Invalid variant description'
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue							

					if re.search('n.1-', str(input_parses)):
						error = 'Using a transcript reference sequence to specify an intergenic variant position that lies 5' + "'" + ' to the transcript reference sequence is not HGVS compliant. Instead use '
						genomic_position = va_func.myevm_t_to_g(input_parses, evm, hdp, primary_assembly)
						error = error + str(genomic_position)
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue
					else:
						pass

					if re.search('\d\-', str(input_parses)) or re.search('\d\+', str(input_parses)):  
						# Create a specific minimal evm with no normalizer and no replace_reference
						min_evm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=primary_assembly, alt_aln_method=alt_aln_method, normalize=False, replace_reference=False)
						# Have to use this method due to potential multi chromosome error, note, normalizes but does not replace sequence
						try:
							output = va_func.noreplace_myevm_t_to_g(input_parses, min_evm, hdp, primary_assembly)
						except hgvs.exceptions.HGVSDataNotAvailableError as e:
							tx_ac = input_parses.ac
							try:
								gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(tx_ac)
							except:
								gene_symbol = 'none'	
							if gene_symbol == 'none':
								error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting '  + tx_ac + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
							else:
								error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
							continue						
						except ValueError as e:
							error = str(e)
							if re.search('> end', error):
								error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end position ' + str(input_parses.posedit.pos.end)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
						except hgvs.exceptions.HGVSInvalidVariantError as e:
							error = str(e)
							if re.search('base start position must be <= end position', error):
								correction = copy.deepcopy(input_parses)
								st = input_parses.posedit.pos.start
								ed = input_parses.posedit.pos.end
								correction.posedit.pos.start = ed
								correction.posedit.pos.end = st
								error = error + ': Did you mean ' + str(correction) + '?' 
								error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end position ' + str(input_parses.posedit.pos.end)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
						try: 
							vr.validate(output)
						except hgvs.exceptions.HGVSError as e:
							error = str(e)
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
							continue

					else:
						# All other variation 
						try:
							vr.validate(input_parses)
						except hgvs.exceptions.HGVSUnsupportedOperationError:
							pass
						except hgvs.exceptions.HGVSInvalidVariantError as e:
							error = str(e)
							if re.search('Length implied by coordinates', error):
								# Applies to del and inv
								# NOTE, there has been no normalization at all so this error is valid here
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue					
								# Will apply to > del and inv
							if re.search('does not agree with reference sequence', error):
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue					
								# ensures x_y for insertions
							if re.search('insertion length must be 1', error):
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue					
								# This catches errors in introns
							if re.search('base start position must be <= end position', error):
								correction = copy.deepcopy(input_parses)
								st = input_parses.posedit.pos.start
								ed = input_parses.posedit.pos.end
								correction.posedit.pos.start = ed
								correction.posedit.pos.end = st
								error = error + ': Did you mean ' + str(correction) + '?' 
								error = 'Interval start position ' + str(input_parses.posedit.pos.start) + ' > interval end position ' + str(input_parses.posedit.pos.end)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
						except hgvs.exceptions.HGVSDataNotAvailableError as e:
								error = e
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue		
						except hgvs.exceptions.HGVSError as e:
							error = str(e)			
							if re.search('outside the bounds', error):
								error = error + ' (' + input_parses.ac + ')'
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue					
				else:
					pass
		

				# Mitochondrial variants
				if type == ':m.':
					error='Mitochondrial variation is not currently supported'
					validation['warnings'] = validation['warnings'] + ': ' + str(error)
					continue

				# handle :p.
				if type == ':p.':		
					error = 'false'
					# Try to validate the variant
					error = va_func.validate(variant, hp, vr)
					if error != 'false':	
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue
					else:
						# Get accurate descriptions from the relevant databases
						# RefSeq databases
						if alt_aln_method != 'genebuild':		
							# Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
							# accession number
							hgvs_object = hp.parse_hgvs_variant(variant)	
							accession = hgvs_object.ac
							# Look for the accession in our database
							# Connect to database and send request
							record = va_func.entrez_efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
							try:
								description = record.description
							except:
								description = 'Unable to recover the description of ' + accession + ' from Entrez'
							try:
								vr.validate(hgvs_object)
							except hgvs.exceptions.HGVSError as e:
								error = str(e)
							else:
								error = str(hgvs_object) + ' is HGVS compliant and contains a valid reference amino acid description' 
							reason = 'Protein level variant descriptions are not fully supported due to redundancy in the genetic code'
							validation['warnings'] = validation['warnings'] + ': ' + str(reason) + ': ' + str(error)
							continue


				# handle :r.
				trapped_input = input
				if type == ':r.':
					input = hp.parse_hgvs_variant(input) # Traps the hgvs variant of r. for further use
					# Change to coding variant
					type =':c.'
					# Change input to reflect!
					hgvs_c = va_func.hgvs_r_to_c(input)
					input = str(hgvs_c)
					variant = str(hgvs_c)


				# Convert and handle :n.
				# Note, map to :c. currently. May need to look into non-coding transcripts
				if type == ':n.':
					# Convert to hgvs object
					hgvs_n = hp.parse_hgvs_variant(variant)
					# try to map to :c.
					error = 'false'
					try:
						hgvs_c = evm.n_to_c(hgvs_n)
					except hgvs.exceptions.HGVSError as e:
						error = str(e)
					if error == 'false':
						# Re-set input and variant
						old_in = input
						old_var = variant
						input = str(hgvs_c)
						variant = valstr(hgvs_c)
						type = ':c.'
						caution = 'The use of n. to describe a variant with respect to a protein-coding transcript does not comply with HGVS variant nomenclature' 
						automap = old_in + ' automapped to ' + variant
						validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(automap)

					else:
						# Keep everything the same
						pass

				# Select for requested Genes if any
				if len(input_genes) == 0:
					# If no genes selected, simply add the VCF call
					pass
				else:
					# Map to genomic position
					if type != ':g.':
						genomic_position = va_func.myevm_t_to_g(input_parses, evm, hdp, primary_assembly)
					else:
						genomic_position = input_parses	
					# Select the genomic positions
					start_pos = genomic_position.posedit.pos.start.base
					end_pos = genomic_position.posedit.pos.end.base
					# Do we keep it?
					keeper = 'false'
					for gene in input_genes:
						if int(start_pos) >= int(input_genes[gene]['start']) and int(start_pos) <= int(input_genes[gene]['end']):
							keeper = 'true'
							break
						elif int(end_pos) >= int(input_genes[gene]['start']) and int(end_pos) <= int(input_genes[gene]['end']):
							keeper = 'true'
							break
						else:
							keeper = 'false'
							continue
					# If not get rid of it!
					if keeper == 'false':
						# By marking it as Do Not Write and continuing through the validation loop
						validation['write'] = 'false'
						continue
		
				# COLLECT gene symbol, name and ACCESSION INFORMATION
				# Gene symbol
				if (type != ':g.'):
					error = 'false'
					hgvs_vt= hp.parse_hgvs_variant(variant)
					try:
						tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
					except hgvs.exceptions.HGVSError as e:
						error = str(e)
					if error != 'false':
						error = 'Please inform UTA admin of the following error: ' + str(error)
						issue_link = "https://bitbucket.org/biocommons/uta/issues?status=new&status=open"
						reason = "VariantValidator can not recover information for transcript " + str(hgvs_vt.ac) + ' beacuse it is not available in the Universal Transcript Archive'
						validation['warnings'] = validation['warnings'] + ': ' + str(reason)
						continue
					else:
						# Get hgnc Gene name from command
						hgnc = tx_id_info[6]
						issue_link = 'false'

					# ACCESS THE GENE INFORMATION RECORDS ON THE UTA DATABASE 
					# Refseq accession
					tx_for_gene = va_func.tx_for_gene(hgnc, hdp)
					refseq_ac = va_func.ng_extract(tx_for_gene)

					# Additional gene info
					gene_info = hdp.get_gene_info(hgnc)
					# Chromosomal location
					try:
						maploc = gene_info[1]
					except:
						maploc = ''	
					chr_loc = ("Chromosome location: " + maploc)

					# Get accurate transcript descriptions from the relevant databases
					# RefSeq databases
					if alt_aln_method != 'genebuild':		
						# Gene description  - requires GenBank search to get all the required info, i.e. transcript variant ID
						# accession number
						hgvs_object = hp.parse_hgvs_variant(variant)	
						accession = hgvs_object.ac
						# Look for the accession in our database
						# Connect to database and send request
						entry = va_dbCrl.data.in_entries(accession.split(',')[0], 'transcript_info')

						# Analyse the returned data and take the necessary actions
						# If the error key exists
						if 'error' in entry:
							# Open a hgvs exception log file in append mode
							error=entry['description']
							validation['warnings'] = validation['warnings'] + ': ' + str(error) + ': A Database error occurred, please contact admin'
							continue

						# If the accession key is found
						elif 'accession' in entry:
							description = entry['description']
							# If the current entry is too old
							if entry['expiry'] == 'true':
								dbaction = 'update'
								entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method, accession=accession, dbaction=dbaction, hp=hp, evm=evm, hdp=hdp)
								hgnc_gene_info = entry['description']
							else:
								hgnc_gene_info = entry['description']
						# If the none key is found add the description to the database
						elif 'none' in entry:
							dbaction = 'insert'
							try:
								entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method, accession=accession, dbaction=dbaction, hp=hp, evm=evm, hdp=hdp)
							except:
								error = 'Unable to assign transcript identity records to ' + accession + ', potentially an obsolete record :'
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
							hgnc_gene_info = entry['description']
				
						# If no correct keys are found
						else:
							# Open a hgvs exception log file in append mode
							error='Unknown error type'
							validation['warnings'] = validation['warnings'] + ': ' + str(error) + ': A Database error occurred, please contact admin'
							continue
		
					# Ensembl databases
					else:
						# accession number
						hgvs_object = hp.parse_hgvs_variant(variant)	
						accession = hgvs_object.ac
						# Look for the accession in our database
						# Connect to database and send request
						entry = va_dbCrl.data.in_entries(accession.split(',')[0], 'transcript_info')
		
						# Analyse the returned data and take the necessary actions
						# If the error key exists
						if 'error' in entry:
							# Open a hgvs exception log file in append mode
							error=entry['description']
							validation['warnings'] = validation['warnings'] + ': ' + str(error) + ': A Database error occurred, please contact admin'
							continue

						# If the accession key is found
						elif 'accession' in entry:
							description = entry['description']
							# If the current entry is too old
							if entry['expiry'] == 'true':
								dbaction = 'update'
								entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method, accession=accession, dbaction=dbaction, hp=hp, evm=evm, hdp=hdp)
								hgnc_gene_info = entry['description']
							else:
								hgnc_gene_info = entry['description']
						# If the none key is found add the description to the database
						elif 'none' in entry:
							dbaction = 'insert'
							try:
								entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method, accession=accession, dbaction=dbaction, hp=hp, evm=evm, hdp=hdp)
							except:
								error = 'Unable to assign transcript identity records to ' + accession + ', potentially an obsolete record :'
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
							hgnc_gene_info = entry['description']
				
						# If no correct keys are found
						else:
							# Open a hgvs exception log file in append mode
							error='Unknown error type'
							validation['warnings'] = validation['warnings'] + ': ' + str(error) + ': A Database error occurred, please contact admin'
							continue
					

				# Genomic type variants will need to be mapped to transcripts

				if (type == ':g.'):
					g_query = hp.parse_hgvs_variant(variant)
	
					# Genomic coordinates can be validated immediately
					error = 'false'
					try:
						vr.validate(g_query)
					except hgvs.exceptions.HGVSError as e:
						error = str(e)
					if error != 'false':
						reason = 'Invalid variant description'
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue
					else:
						pass
		
					# Set test to see if Norm alters the coords
					g_test = hn.normalize(g_query)
	
					# Perform test
					if g_query.posedit.pos != g_test.posedit.pos:
						validation['warnings'] = validation['warnings'] + ': ' + 'Input variant description normalized to ' + str(g_test)
						hgvs_genomic = g_test
					else:
						hgvs_genomic = g_query
	
					# Collect rel_var
					# rel_var is a keyworded list of relevant transcripts with associated coding variants 
					rel_var = va_func.relevant_transcripts(hgvs_genomic, evm, hdp, alt_aln_method)
					
					# Double check rel_vars have not been missed when mapping from a RefSeqGene
					if len(rel_var) != 0 and re.match('NG_', str(hgvs_genomic.ac)):
						alt_rel_var = []
						for found in rel_var:			
							hgvs_coding_variant = hp.parse_hgvs_variant(rel_var[0])
							try:
								hgvs_genomic = va_func.myevm_t_to_g(hgvs_coding_variant, evm, hdp, primary_assembly='GRCh37')	
							except hgvs.exceptions.HGVSDataNotAvailableError as e:
								try_rel_var = []
							else:
								try_rel_var = va_func.relevant_transcripts(hgvs_genomic, evm, hdp, alt_aln_method)	
							if len(try_rel_var) > len(rel_var) and len(try_rel_var) > alt_rel_var:
								alt_rel_var = try_rel_var
							# Re set
							try_rel_var = []		
							try:
								hgvs_genomic = va_func.myevm_t_to_g(hgvs_coding_variant, evm, hdp, primary_assembly='GRCh38')	
							except hgvs.exceptions.HGVSDataNotAvailableError as e:
								try_rel_var = []
							else:
								try_rel_var = va_func.relevant_transcripts(hgvs_genomic, evm, hdp, alt_aln_method)	
							if len(try_rel_var) > len(rel_var) and len(try_rel_var) > alt_rel_var:
								alt_rel_var = try_rel_var
						# re-set rel_var
						rel_var = alt_rel_var					

				
					# Tripple check this assumption by querying the gene position database
					if len(rel_var) == 0:
						vcf_dict = va_H2V.hgvs2vcf(hgvs_genomic)
						not_di = str(hgvs_genomic.ac) + ':g.' + str(vcf_dict['pos']) + '_' + str(int(vcf_dict['pos']) + (len(vcf_dict['ref']) -1)) + 'del' + vcf_dict['ref'] + 'ins' + vcf_dict['alt'] 
						hgvs_not_di = hp.parse_hgvs_variant(not_di)
						rel_var = va_func.relevant_transcripts(hgvs_not_di, evm, hdp, alt_aln_method)					

					# list return statements
					if len(rel_var) == 0:
	
						# Check for NG_
						rsg = re.compile('^NG_')
						if rsg.search(variant):
							# parse
							hgvs_refseqgene = hp.parse_hgvs_variant(variant)
							# Convert to chromosomal position
							refseqgene_data = va_g2g.rsg_to_chr(hgvs_refseqgene, primary_assembly, hn, vr)
							# There should only ever be one description returned
							refseqgene_data = refseqgene_data[0]

							# Extract data
							if refseqgene_data['valid'] == 'true':
								input = refseqgene_data['hgvs_genomic'] 
								# re_submit
								# Tag the line so that it is not written out
								validation['warnings'] = validation['warnings'] + ': ' + variant + ' automapped to genome position ' + str(input)
								query = {'quibble' : input, 'id' : validation['id'], 'warnings' : validation['warnings'], 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
								coding = 'intergenic'
								batch_list.append(query)
							else:
								error = 'Mapping unavailable for RefSeqGene ' + variant + ' using alignment method = ' + alt_aln_method
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue

						# Chromosome build is not supported or intergenic??? 
						else:
							sfm = va_scb.supported_for_mapping(hgvs_genomic.ac, primary_assembly)
							if sfm == 'true':
								try:
									vr.validate(hgvs_genomic)
								except hgvs.exceptions.HGVSError as e:
									error = str(e)
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									continue
								else:
									# Map to RefSeqGene if available
									refseqgene_data = va_g2g.chr_to_rsg(hgvs_genomic, hn, vr)
									rsg_data = ''
									# Example {'gene': 'NTHL1', 'hgvs_refseqgene': 'NG_008412.1:g.3455_3464delCAAACACACA', 'valid': 'true'}
									for data in refseqgene_data:
										if data['valid'] == 'true':
											data['hgvs_refseqgene'] = hp.parse_hgvs_variant(data['hgvs_refseqgene'])
											data['hgvs_refseqgene'] = valstr(data['hgvs_refseqgene'])
											rsg_data = rsg_data + data['hgvs_refseqgene'] + ' (' + data['gene'] + '), '
						
									error = 'Suspected intergenic region, No transcripts fully overlap the input genomic coordinates'
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									validation['genomic_g'] = valstr(hgvs_genomic)
									validation['genomic_r'] = str(rsg_data.split('(')[0])
									continue
							else:
								error = 'Please ensure the requested chromosome version relates to a supported genome build. Supported genome builds are: GRCh37, GRCh38 and NCBI36'
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue		
			
					else:
						# Tag the line so that it is not written out
						validation['write'] = 'false'

						# Set variables for problem specific warnings
						gapped_alignment_warning = ''
						corrective_action_taken = ''
						gapped_transcripts = ''
						auto_info = ''
			
						# Create a pseudo VCF so that normalization can be applied and a delins can be generated
						hgvs_genomic_variant = hgvs_genomic					
						# Reverse normalize hgvs_genomic_variant: NOTE will replace ref
						reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic_variant)
						hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
						
						# VCF
						vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic)
						chr = vcf_dict['chr']	
						pos = vcf_dict['pos']				
						ref = vcf_dict['ref']
						alt = vcf_dict['alt']									
			
						# Generate an end position
						end = str(int(pos) + len(ref) -1)
						pos = str(pos)
			
						# Store a not real deletion insertion
						stored_hgvs_not_delins = hp.parse_hgvs_variant(str(hgvs_genomic_5pr.ac) + ':' +  hgvs_genomic_5pr.type + '.' +  pos + '_' + end + 'del' + ref + 'ins' + alt)
						
			
						# Set non-valid caution to false
						non_valid_caution = 'false'
			
						# make an empty rel_var
						nw_rel_var = []
			
						# loop through rel_var and amend where required
						for var in rel_var:
							# Store the current hgvs:c. description
							saved_hgvs_coding = hp.parse_hgvs_variant(var)							

							# Remove un-selected transcripts
							if select_transcripts != 'all':
								tx_ac = saved_hgvs_coding.ac
								# If it's in the selected tx dict, keep it
								if tx_ac.split('.')[0] in select_transcripts_dict.keys():
									pass
								# If not get rid of it!
								else:
									continue							
							
							# print '\nstored_gen_pos'
							# print stored_hgvs_not_delins

							# Get orientation of the gene wrt genome and a list of exons mapped to the genome
							ori = va_func.tx_exons(tx_ac=saved_hgvs_coding.ac, alt_ac=hgvs_genomic_5pr.ac, alt_aln_method=alt_aln_method, hdp=hdp)
							orientation = int(ori[0]['alt_strand'])
							intronic_variant = 'false'

							if orientation == -1:
								# position genomic at its most 5 prime position
								try:
									query_genomic = reverse_normalize.normalize(hgvs_genomic)
								except:
									query_genomic = hgvs_genomic	
								# Map to the transcript ant test for movement
								try:
									hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
								except hgvs.exceptions.HGVSError as e:
									hgvs_seek_var = saved_hgvs_coding
								else:		
									seek_var = valstr(hgvs_seek_var)
									seek_ac = str(hgvs_seek_var.ac)
								if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (saved_hgvs_coding.posedit.pos.start.base +  saved_hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (saved_hgvs_coding.posedit.pos.end.base +  saved_hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
									pass
								else:
									hgvs_seek_var = saved_hgvs_coding			

							elif orientation != -1:
								# position genomic at its most 3 prime position
								try:
									query_genomic = hn.normalize(hgvs_genomic)
								except:
									query_genomic = hgvs_genomic	
							# Map to the transcript and test for movement
							try:
								hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
							except hgvs.exceptions.HGVSError as e:
								hgvs_seek_var = saved_hgvs_coding
							else:		
								seek_var = valstr(hgvs_seek_var)
								seek_ac = str(hgvs_seek_var.ac)
								if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (saved_hgvs_coding.posedit.pos.start.base +  saved_hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (saved_hgvs_coding.posedit.pos.end.base +  saved_hgvs_coding.posedit.pos.end.offset):
									pass
								else:
									hgvs_seek_var = saved_hgvs_coding	

							try:	
								intron_test = hn.normalize(hgvs_seek_var)	
							except hgvs.exceptions.HGVSUnsupportedOperationError as e:
								error = str(e)
								if re.match('Normalization of intronic variants is not supported', error):
									# Double check to see whether the variant is actually intronic?
									for exon in ori:
										genomic_start = int(exon['alt_start_i'])
										genomic_end = int(exon['alt_end_i'])
										if (hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
											intronic_variant = 'false'
											break							
										else:							
											intronic_variant = 'true'
				
							if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-', str(hgvs_seek_var.posedit.pos)):
								# Double check to see whether the variant is actually intronic?
								for exon in ori:
									genomic_start = int(exon['alt_start_i'])
									genomic_end = int(exon['alt_end_i'])
									if (hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
										intronic_variant = 'false'
										break							
									else:							
										intronic_variant = 'true'

							if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-', str(hgvs_seek_var.posedit.pos)):
								# Double check to see whether the variant is actually intronic?
								for exon in ori:
									genomic_start = int(exon['alt_start_i'])
									genomic_end = int(exon['alt_end_i'])
									if (hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
										intronic_variant = 'false'
										break							
									else:							
										intronic_variant = 'true'										

							# print hgvs_seek_var				
							# print 'intronic_variant = ' + intronic_variant


							# If exonic, process
							if intronic_variant != 'true':	
								# map form reverse normalized g. to c. 
								hgvs_from_5n_g = no_norm_evm.g_to_t(hgvs_genomic_5pr, saved_hgvs_coding.ac)
					
								# Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
								disparity_deletion_in = ['false', 'false']
								if stored_hgvs_not_delins != '':
									# Refresh hgvs_not_delins from stored_hgvs_not_delins
									hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
									# This test will only occur in dup of single base, insertion or substitution
									if not re.search('_', str(hgvs_not_delins.posedit.pos)):
										# For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
										plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
										plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
										plussed_hgvs_not_delins.posedit.edit.ref = ''
										transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins, str(saved_hgvs_coding.ac))
										if ((transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
											if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
												start = hgvs_not_delins.posedit.pos.start.base - 1
												end = hgvs_not_delins.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
												hgvs_not_delins.posedit.edit.ref = ref_bases
												hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]
											elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
											elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and not re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
												start = hgvs_not_delins.posedit.pos.start.base - 1
												end = hgvs_not_delins.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
												hgvs_not_delins.posedit.edit.ref = ref_bases
												hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]						
										else:
											if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
												start = hgvs_not_delins.posedit.pos.start.base - 1
												end = hgvs_not_delins.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
												hgvs_not_delins.posedit.edit.ref = ref_bases
												hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]
											elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
											elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and not re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
												start = hgvs_not_delins.posedit.pos.start.base - 1
												end = hgvs_not_delins.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
												hgvs_not_delins.posedit.edit.ref = ref_bases
												hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] # + ref_bases[1:]						
									else:
										pass	
				
									tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
									# Create normalized version of tx_hgvs_not_delins
									rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
									# Check for +ve base and adjust
									if (re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start))) and (re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) or re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end))):
										# Remove offsetting to span the gap
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										try:
											rn_tx_hgvs_not_delins.posedit.edit.alt = '' 
										except:
											pass
									elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
										# move tx end base to next available non-offset base
										rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
										rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										if re.match('NM_', str(rn_tx_hgvs_not_delins)):
											test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
										else:
											test_tx_var = rn_tx_hgvs_not_delins	
										# re-make genomic and tx
										hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
										rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
									elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
										# move tx start base to previous available non-offset base
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										if re.match('NM_', str(rn_tx_hgvs_not_delins)):
											test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
										else:
											test_tx_var = rn_tx_hgvs_not_delins	
										# re-make genomic and tx
										hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
										rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0							
									else:
										pass

									# Check for -ve base and adjust
									if re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
										# Remove offsetting to span the gap
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										try:
											rn_tx_hgvs_not_delins.posedit.edit.alt = '' 
										except:
											pass
									elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
										# move tx end base back to next available non-offset base
										rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
										# Delete the ref
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										# Add the additional base to the ALT
										start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
										end = rn_tx_hgvs_not_delins.posedit.pos.end.base
										ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac),start,end)
										rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
										if re.match('NM_', str(rn_tx_hgvs_not_delins)):
											test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
										else:
											test_tx_var = rn_tx_hgvs_not_delins	
										# re-make genomic and tx
										hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
										rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
									elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
										# move tx start base to previous available non-offset base
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										if re.match('NM_', str(rn_tx_hgvs_not_delins)):
											test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
										else:
											test_tx_var = rn_tx_hgvs_not_delins	
										# re-make genomic and tx
										hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
										rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0							
									else:
										pass

									# Logic
									if len(hgvs_not_delins.posedit.edit.ref) < len(rn_tx_hgvs_not_delins.posedit.edit.ref):
										gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(hgvs_not_delins.posedit.edit.ref)
										disparity_deletion_in = ['chromosome', gap_length]
									elif len(hgvs_not_delins.posedit.edit.ref) > len(rn_tx_hgvs_not_delins.posedit.edit.ref):
										gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(rn_tx_hgvs_not_delins.posedit.edit.ref)
										disparity_deletion_in = ['transcript', gap_length]
									else:
										pass
		
								# print disparity_deletion_in
								
								
								# GAP IN THE TRANSCRIPT	DISPARITY DETECTED	
								if disparity_deletion_in[0] == 'transcript':			
									gap_position = ''
									gapped_alignment_warning = str(hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + primary_assembly

									# ANY VARIANT WHOLLY WITHIN THE GAP
									if (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start))) and (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end))):
										non_valid_caution = 'true'
										auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos) + ' aligns to a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)						
										gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
										gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
										if str(tx_hgvs_not_delins.posedit.pos.start) == str(tx_hgvs_not_delins.posedit.pos.end):
											if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
												tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
												tx_hgvs_not_delins.posedit.pos.end.offset = 0
											elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
												tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
												tx_hgvs_not_delins.posedit.pos.start.offset = 0									
										if re.search('\=', str(tx_hgvs_not_delins.posedit)):
											middle = tx_hgvs_not_delins.posedit.edit.alt
										elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)):
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
										elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
											middle = tx_hgvs_not_delins.posedit.edit.alt#[1:]
										elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
											middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
											tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
										tx_hgvs_not_delins.posedit.pos.start.offset = 0
										tx_hgvs_not_delins.posedit.pos.end.offset = 0							
										beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base)
										end = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
										tx_hgvs_not_delins.posedit.edit.ref = beginning + end
										tx_hgvs_not_delins.posedit.edit.alt = beginning + middle + end
										# map to the genome and back						
										tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										else:
											c = tx_hgvs_not_delins	
										g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
										c = no_norm_evm.g_to_t(g, str(c.ac))
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
										else:
											tx_hgvs_not_delins = c	
										hgvs_refreshed_variant = tx_hgvs_not_delins														
										# Alignment position
										for_location_c = copy.deepcopy(hgvs_refreshed_variant)
										if re.match('NM_', str(for_location_c)):
											for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										gps = for_location_c.posedit.pos.start.base
										gpe = for_location_c.posedit.pos.start.base + 1
										gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
										auto_info = auto_info + '%s' %(gap_position)
									else:
										if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)):
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.start.base) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
											non_valid_caution = 'true'
											tx_hgvs_not_delins.posedit.pos.start.offset = 0
											beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base) 
											middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
											end = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base, tx_hgvs_not_delins.posedit.pos.end.base)
											tx_hgvs_not_delins.posedit.edit.ref = beginning + end
											tx_hgvs_not_delins.posedit.edit.alt = beginning + middle
											# map to the genome and back						
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
											if re.match('NM_', str(tx_hgvs_not_delins.ac)):
												c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											else:
												c = tx_hgvs_not_delins	
											g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
											c = no_norm_evm.g_to_t(g, str(c.ac))
											if re.match('NM_', str(tx_hgvs_not_delins.ac)):
												tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
											else:
												tx_hgvs_not_delins = c									
											hgvs_refreshed_variant = tx_hgvs_not_delins
											# Alignment position
											for_location_c = copy.deepcopy(hgvs_refreshed_variant)
											if re.match('NM_', str(for_location_c)):
												for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											gps = for_location_c.posedit.pos.start.base -1
											gpe = for_location_c.posedit.pos.start.base
											gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
											# Warn update
											auto_info = auto_info + '%s' %(gap_position)
										elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)						
											gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
											non_valid_caution = 'true'
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
											tx_hgvs_not_delins.posedit.pos.end.offset = 0
											beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.end.base - 1)
											middle = tx_hgvs_not_delins.posedit.edit.alt
											end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
											tx_hgvs_not_delins.posedit.edit.ref = beginning + end
											tx_hgvs_not_delins.posedit.edit.alt = middle + end
											# map to the genome and back						
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
											if re.match('NM_', str(tx_hgvs_not_delins.ac)):
												c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											else:
												c = tx_hgvs_not_delins	
											g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
											c = no_norm_evm.g_to_t(g, str(c.ac))
											if re.match('NM_', str(tx_hgvs_not_delins.ac)):
												tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
											else:
												tx_hgvs_not_delins = c	
											hgvs_refreshed_variant = tx_hgvs_not_delins				
											# Alignment position
											for_location_c = copy.deepcopy(hgvs_refreshed_variant)
											if re.match('NM_', str(for_location_c)):
												for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											gps = for_location_c.posedit.pos.end.base - 1
											gpe = for_location_c.posedit.pos.end.base
											gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
											# Warn update
											auto_info = auto_info + '%s' %(gap_position)
										elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)):
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.start.base) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
											non_valid_caution = 'true'
											tx_hgvs_not_delins.posedit.pos.start.offset = 0
											tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
											beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base)
											middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
											end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base, tx_hgvs_not_delins.posedit.pos.end.base)
											if not re.search('\=', str(tx_hgvs_not_delins.posedit)):
												tx_hgvs_not_delins.posedit.edit.alt = beginning + middle
											else:
												tx_hgvs_not_delins.posedit.edit.alt = beginning + middle + end
											tx_hgvs_not_delins.posedit.edit.ref = beginning + end
											# map to the genome and back						
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
											if re.match('NM_', str(tx_hgvs_not_delins.ac)):
												c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											else:
												c = tx_hgvs_not_delins	
											g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
											c = no_norm_evm.g_to_t(g, str(c.ac))
											if re.match('NM_', str(tx_hgvs_not_delins.ac)):
												tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
											else:
												tx_hgvs_not_delins = c	
											hgvs_refreshed_variant = tx_hgvs_not_delins								
											# Alignment position
											for_location_c = copy.deepcopy(hgvs_refreshed_variant)
											if re.match('NM_', str(for_location_c)):
												for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											gps = for_location_c.posedit.pos.start.base -1
											gpe = for_location_c.posedit.pos.start.base
											gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
											# Warn update
											auto_info = auto_info + '%s' %(gap_position)
										elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)						
											gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
											non_valid_caution = 'true'
											tx_hgvs_not_delins.posedit.pos.end.offset = 0
											beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.end.base - 1)
											if len(tx_hgvs_not_delins.posedit.edit.ref) < len(tx_hgvs_not_delins.posedit.edit.alt):
												middle = tx_hgvs_not_delins.posedit.edit.alt.replace(tx_hgvs_not_delins.posedit.edit.ref, '' , 1)
												middle = beginning[0] + middle
											elif len(tx_hgvs_not_delins.posedit.edit.ref) > len(tx_hgvs_not_delins.posedit.edit.alt):
												middle = tx_hgvs_not_delins.posedit.edit.alt
											elif len(tx_hgvs_not_delins.posedit.edit.ref) == len(tx_hgvs_not_delins.posedit.edit.alt):
												middle = tx_hgvs_not_delins.posedit.edit.alt																
											end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
											tx_hgvs_not_delins.posedit.edit.ref = beginning + end
											tx_hgvs_not_delins.posedit.edit.alt = middle + end
											# map to the genome and back						
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
											if re.match('NM_', str(tx_hgvs_not_delins.ac)):
												c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											else:
												c = tx_hgvs_not_delins	
											g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
											c = no_norm_evm.g_to_t(g, str(c.ac))
											if re.match('NM_', str(tx_hgvs_not_delins.ac)):
												tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
											else:
												tx_hgvs_not_delins = c	
											hgvs_refreshed_variant = tx_hgvs_not_delins								
											# Alignment position
											for_location_c = copy.deepcopy(hgvs_refreshed_variant)
											if re.match('NM_', str(for_location_c)):
												for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											gps = for_location_c.posedit.pos.end.base - 1
											gpe = for_location_c.posedit.pos.end.base
											gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
											# Warn update
											auto_info = auto_info + '%s' %(gap_position)
										else:
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos) + ' aligns across a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac) + '\n'
											hgvs_refreshed_variant = tx_hgvs_not_delins
											gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
					
								# GAP IN THE CHROMOSOME		
								elif disparity_deletion_in[0] == 'chromosome':
									# Set warning variables
									gap_position = ''
									gapped_alignment_warning = str(hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + primary_assembly
									hgvs_refreshed_variant = tx_hgvs_not_delins
									gapped_transcripts = gapped_transcripts + str(hgvs_refreshed_variant.ac) + ' '	
									# print hgvs_refreshed_variant
								else:
									# Keep the same by re-setting rel_var
									hgvs_refreshed_variant = saved_hgvs_coding
				
								# Edit the output
								if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(hgvs_refreshed_variant.type)):
									hgvs_refreshed_variant = evm.n_to_c(hgvs_refreshed_variant) 
								else:
									pass
								try:
									hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
									if hgvs_refreshed_variant.posedit.edit.type == 'delins' and hgvs_refreshed_variant.posedit.edit.ref[-1] == hgvs_refreshed_variant.posedit.edit.alt[-1]:
										hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[0:-1]
										hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[0:-1]
										hgvs_refreshed_variant.posedit.pos.end.base = hgvs_refreshed_variant.posedit.pos.end.base - 1
										hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
									elif hgvs_refreshed_variant.posedit.edit.type == 'delins' and hgvs_refreshed_variant.posedit.edit.ref[0] == hgvs_refreshed_variant.posedit.edit.alt[0]:
										hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[1:]
										hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[1:]
										hgvs_refreshed_variant.posedit.pos.start.base = hgvs_refreshed_variant.posedit.pos.start.base + 1
										hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)	
								except:
									pass
								
								# print 'in ' + str(var)
								# print 'out ' + str(hgvs_refreshed_variant)
								# Send to empty nw_rel_var
								nw_rel_var.append(hgvs_refreshed_variant)
				
							# Otherwise these variants need to be set
							else:
								corrective_action_taken = ''
								gapped_alignment_warning = ''
								# Send to empty nw_rel_var
								nw_rel_var.append(saved_hgvs_coding)
						
						# Warn the user that the g. description is not valid
						if gapped_alignment_warning != '':								
							if disparity_deletion_in[0] == 'transcript': 
								corrective_action_taken = 'Automap has deleted  ' + str(disparity_deletion_in[1]) + ' bp from chromosomal reference sequence ' + str(hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference sequence(s)' + gapped_transcripts
							if disparity_deletion_in[0] == 'chromosome':
								corrective_action_taken = 'Automap has added  ' + str(disparity_deletion_in[1]) + ' bp to chromosomal reference sequence ' + str(hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference sequence(s) ' + gapped_transcripts

						# Add additional data to the front of automap
						if auto_info != '':
							automap = auto_info + '\n' + automap

						#
						rel_var = copy.deepcopy(nw_rel_var)
						# print rel_var
						
						# Set the values and append to batch_list
						for c_description in rel_var:
# 							if select_transcripts != 'all':
# 								hgvs_c_description = hp.parse_hgvs_variant(c_description)
# 								tx_ac = hgvs_c_description.ac
# 								if tx_ac.split('.')[0] in select_transcripts_dict.keys():
# 									pass
# 								# If not get rid of it!
# 								else:
# 									continue
							query = {'quibble' : str(c_description), 'id' : validation['id'], 'warnings' : validation['warnings'], 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
							batch_list.append(query)		
		
						# Call next description
						continue

				# TYPE = :c.

				if type == ':c.' or type == ':n.':
					# Flag for validation 
					valid = 'false'
					# Collect information for genomic level validation
					obj = hp.parse_hgvs_variant(variant)	
					tx_ac = obj.ac
			
					# Do we keep it?
					if select_transcripts != 'all':
						if tx_ac in select_transcripts_dict_plus_version.keys():	
							pass
						# If not get rid of it!
						else:
							# By marking it as Do Not Write and continuing through the validation loop
							validation['write'] = 'false'
							continue
					else:
						pass	
			
					# Set a cross_variant object
					cross_variant = 'false'
					# Se rec_var to '' so it can be updated later
					rec_var = ''
					try:
						to_g = va_func.myevm_t_to_g(obj, evm, hdp, primary_assembly)
						genomic_ac = to_g.ac
					except hgvs.exceptions.HGVSDataNotAvailableError as e:
						try:
							gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(tx_ac)
						except:
							gene_symbol = 'none'	
						if gene_symbol == 'none':
							error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting '  + tx_ac + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
						else:
							error = 'Required information for ' + tx_ac + ' is missing from the Universal Transcript Archive, please select an alternative version of ' + tx_ac + ' by submitting ' + tx_ac + ' or ' + gene_symbol + ' to  https://variantvalidator.org/ref_finder/, or select an alternative genome build'
						validation['warnings'] = validation['warnings'] + ': ' + str(error)
						continue

					# Get orientation of the gene wrt genome and a list of exons mapped to the genome
					ori = va_func.tx_exons(tx_ac=tx_ac, alt_ac=genomic_ac, alt_aln_method=alt_aln_method, hdp=hdp)
					orientation = int(ori[0]['alt_strand'])
					intronic_variant = 'false'					
					
					# Collect variant sequence information via normalisation (normalizer) or if intronic via mapping
					# INTRONIC OFFSETS - Required for Exon table
					# Variable to collect offset to exon boundary
					ex_offset = 0
					plus = re.compile("\d\+\d") 	# finds digit + digit
					minus = re.compile("\d\-\d")	# finds digit - digit

					geno = re.compile(':g.')
					if plus.search(input) or minus.search(input):
						to_g = va_func.genomic(variant, evm, hp, hdp, primary_assembly)
						es = re.compile('error')
						if es.search(str(to_g)):
							if alt_aln_method != 'genebuild':
								error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
								reason = "An error has occurred"
								excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue

							else:
								error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
								reason = "An error has occurred"
								excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue			
	
						else:
							# Normalise the g variant
							to_g = hn.normalize(to_g)
							variant = str(va_func.myevm_g_to_t(hdp, evm, to_g, tx_ac))
							tx_ac = ''
					elif geno.search(input):
						if plus.search(variant) or minus.search(variant):
							to_g = va_func.genomic(variant, evm, hp, hdp, primary_assembly)
							es = re.compile('error')
							if es.search(str(to_g)):
								if alt_aln_method != 'genebuild':
									error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
									reason = "An error has occurred"
									excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									continue					
			
								else:
									error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
									reason = "An error has occurred"
									excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									continue			
					
						else:
							# Normalise the g variant
							to_g = hp.parse_hgvs_variant(to_g)
							variant = str(va_func.myevm_g_to_t(hdp, evm, to_g, tx_ac))
							tx_ac = ''	

					else:	
						# Normalize the variant
						error = 'false'
						try:
							h_variant = hn.normalize(obj)
						except hgvs.exceptions.HGVSError as e:
							error = str(e) 
							h_variant = obj
							variant = variant
							caution = 'This coding sequence variant description spans at least one intron'
							automap = 'Use of the coresponding genomic sequence variant descriptions may be invalid. Please refer to the Complex nomenclature tab'
						else:
							variant = str(h_variant)
							
						tx_ac = ''
						# Create a crosser (exon boundary crossed) variant
						crossed_variant = str(evm._maybe_normalize(obj))
						if variant == crossed_variant:
							cross_variant = 'false'
						else:
							hgvs_crossed_variant = evm._maybe_normalize(obj)
							cross_variant = ["Coding sequence allowing for exon boundary crossing (default = no crossing)",	crossed_variant, hgvs_crossed_variant.ac]
							cr_available = 'true'
	
						# control of cross_variant
						if boundary == 'false':
							cross_variant = 'false'
		
							error = va_func.validate(variant, hp=hp, vr=vr)
							if error == 'false':
								valid = 'true'
							else:
								excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
							
					# Tackle the plus intronic offset
					cck = 'false'
					if (plus.search(input)):
						# Regular expression catches the start of the interval only based on .00+00 pattern
						inv_start = re.compile("\.\d+\+\d")
						if (inv_start.search(input)):
							# Find pattern e.g. +0000 and assign to a variable
							off_value = re.search(r"(\+\d+)", input)
							off_value = off_value.group(1)
							# Integerise the value and assign to ex_offset
							ex_offset = int(off_value)
							cck = 'true'
					if (minus.search(input)):
						# Regular expression catches the start of the interval only based on .00-00 pattern
						inv_start = re.compile("\.\d+\-\d")
						if (inv_start.search(input)):
							# Find pattern e.g. -0000 and assign to a variable
							off_value = re.search(r"(\-\d+)", input)
							off_value = off_value.group(1)
							# Integerise the value and assign to ex_offset
							ex_offset = int(off_value)
							cck = 'true'

					# COORDINATE CHECKER
					# hgvs will handle incorrect coordinates so need to automap errors
					# Make sure any input intronic coordinates are correct
					# Get the desired transcript
					pat_r = re.compile(':r.')
					pat_g = re.compile(':g.')
					if cck == 'true':
						dl = re.compile('del')
						# This should only ever hit coding and RNA variants
						if dl.search(variant):
							# RNA
							if pat_r.search(trapped_input):

								coding = va_func.coding(variant, hp)
								trans_acc = coding.ac
								# c to Genome coordinates - Map the variant to the genome
								pre_var = va_func.genomic(variant, evm, hp, hdp, primary_assembly)
								# genome back to C coordinates
								post_var = va_func.myevm_g_to_t(hdp, evm, pre_var, trans_acc)
						
								test = hp.parse_hgvs_variant(input)
								if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
									caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
									automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
									# automapping of variant completed
									# Change to rna variant
									posedit = query.posedit
									posedit = posedit.lower()
									query.posedit = posedit
									query.type = 'r'
									post_var = str(query)
									automap = trapped_input + ' automapped to ' + str(post_var)
									validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(automap)
									relevant = "Select the automapped transcript and click Submit to analyse"
									rel_var = []
									rel_var.append(post_var)
									# Add gene symbols to the link
									cp_rel = copy.copy(rel_var)
									del rel_var[:]
									for accessions in cp_rel:
										error = 'false'
										hgvs_vt = hp.parse_hgvs_variant(str(accessions))
										try:
											tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
										except hgvs.exceptions.HGVSError as e:
											error = str(e)
										if error != 'false':
											accessions = ['', str(hgvs_vt)]
											rel_var.append(accessions)

										else:
										# Get hgnc Gene name from command
											data = va_func.hgnc_rest(path = "/search/prev_symbol/" + tx_id_info[6])
											if data['error'] != 'false':
												reason = 'Cannot currently display the required information:'
												error =  decoded['error']
												validation['warnings'] = validation['warnings'] + ': ' + str(error)
												continue
											else:
												# Set the hgnc name correctly
												# If the name is correct no record will be found
												if int(data['record']['response']['numFound']) == 0:
													current = tx_id_info[6]
												else:
													current = data['record']['response']['docs'][0]['symbol']
											accessions = [str(current), str(hgvs_vt)]
											rel_var.append(accessions)
											# Add gene symbols to the link
											cp_rel = copy.copy(rel_var)
											del rel_var[:]
											for accessions in cp_rel:
												error = 'false'
												hgvs_vt = hp.parse_hgvs_variant(str(accessions[1]))
												try:
													tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
												except hgvs.exceptions.HGVSError as e:
													error = str(e)
												if error != 'false':
													accessions = ['', str(hgvs_vt)]
													rel_var.append(accessions)
												else:
												# Get hgnc Gene name from command
													data = va_func.hgnc_rest(path = "/search/prev_symbol/" + tx_id_info[6])
													if data['error'] != 'false':
														reason = 'Cannot currently display the required information:'
														error =  decoded['error']
														validation['warnings'] = validation['warnings'] + ': ' + str(error)
														continue					
													else:
														# Set the hgnc name correctly
														# If the name is correct no record will be found
														if int(data['record']['response']['numFound']) == 0:
															current = tx_id_info[6]
														else:
															current = data['record']['response']['docs'][0]['symbol']
													accessions = [str(current), str(hgvs_vt)]
													rel_var.append(accessions)
									#Kill current line and append for re-submission
									# Tag the line so that it is not written out
									validation['write'] = 'false'
									# Set the values and append to batch_list
									query = {'quibble' : valstr(hgvs_vt), 'id' : validation['id'], 'warnings' : automap, 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
									batch_list.append(query)
		
							# Coding
							else:
								coding = va_func.coding(variant, hp)
								trans_acc = coding.ac
								# c to Genome coordinates - Map the variant to the genome
								pre_var = hp.parse_hgvs_variant(variant)
								try:
									pre_var = va_func.myevm_t_to_g(pre_var, evm, hdp, primary_assembly)
								except:
									e = sys.exc_info()[1]
									error = str(e)
									reason = 'Input coordinates may be invalid'
									if error == 'expected from_start_i <= from_end_i':
										error = 'Automap is unable to correct the input exon/intron boundary coordinates, please check your variant description'
										validation['warnings'] = validation['warnings'] + ': ' + str(error)
										continue
									else:
										pass	
								else:
									pass
								# genome back to C coordinates
								try:
									post_var = va_func.myevm_g_to_t(hdp, evm, pre_var, trans_acc)
								except hgvs.exceptions.HGVSError as error:
									validation['warnings'] = validation['warnings'] + ': ' +str(error)
									continue
								query = post_var
								test = hp.parse_hgvs_variant(input)
								if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
									caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
									automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
									# automapping of variant completed
									automap = trapped_input + ' automapped to ' + str(post_var)
									validation['warnings'] = str(validation['warnings']) + str(caution) + ': ' + str(automap)
									relevant = "Select the automapped transcript and click Submit to analyse"
									rel_var = []
									rel_var.append(post_var)
									# Add gene symbols to the link
									cp_rel = copy.copy(rel_var)
									del rel_var[:]
									for accessions in cp_rel:
										error = 'false'
										hgvs_vt = hp.parse_hgvs_variant(str(accessions))
										try:
											tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
										except hgvs.exceptions.HGVSError as e:
											error = str(e)
										if error != 'false':
											accessions = ['', str(hgvs_vt)]
											rel_var.append(accessions)
										else:
										# Get hgnc Gene name from command
											data = va_func.hgnc_rest(path = "/search/prev_symbol/" + tx_id_info[6])
											if data['error'] != 'false':
												reason = 'Cannot currently display the required information:'
												error =  decoded['error']
												validation['warnings'] = validation['warnings'] + ': ' + str(error)
												continue
				
											else:
												# Set the hgnc name correctly
												# If the name is correct no record will be found
												if int(data['record']['response']['numFound']) == 0:
													current = tx_id_info[6]
												else:
													current = data['record']['response']['docs'][0]['symbol']
											accessions = [str(current), str(hgvs_vt)]
											rel_var.append(accessions)
											# Add gene symbols to the link
											cp_rel = copy.copy(rel_var)
											del rel_var[:]
											for accessions in cp_rel:
												error = 'false'
												hgvs_vt = hp.parse_hgvs_variant(str(accessions[1]))
												try:
													tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
												except hgvs.exceptions.HGVSError as e:
													error = str(e)
												if error != 'false':
													accessions = ['', str(hgvs_vt)]
													rel_var.append(accessions)
												else:
												# Get hgnc Gene name from command
													data = va_func.hgnc_rest(path = "/search/prev_symbol/" + tx_id_info[6])
													if data['error'] != 'false':
														reason = 'Cannot currently display the required information:'
														error =  decoded['error']
														validation['warnings'] = validation['warnings'] + ': ' + str(error)
														continue
			
													else:
														# Set the hgnc name correctly
														# If the name is correct no record will be found
														if int(data['record']['response']['numFound']) == 0:
															current = tx_id_info[6]
														else:
															current = data['record']['response']['docs'][0]['symbol']
													accessions = [str(current), str(hgvs_vt)]
													rel_var.append(accessions)
									#Kill current line and append for re-submission
									# Tag the line so that it is not written out
									validation['write'] = 'false'
									# Set the values and append to batch_list
									query = {'quibble' : valstr(hgvs_vt), 'id' : validation['id'], 'warnings' : automap, 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
									batch_list.append(query)
	
						else:
							if pat_r.search(trapped_input):
								coding = va_func.coding(variant, hp)
								trans_acc = coding.ac
								# c to Genome coordinates - Map the variant to the genome
								pre_var = va_func.genomic(variant, evm, hp, hdp, primary_assembly)
								# genome back to C coordinates
								post_var = va_func.myevm_g_to_t(hdp, evm, pre_var, trans_acc)
						
								test = hp.parse_hgvs_variant(input)
								if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
									caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
									automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
									# automapping of variant completed
									# Change to rna variant
									posedit = query.posedit
									posedit = posedit.lower()
									query.posedit = posedit
									query.type = 'r'
									post_var = str(query)
									automap = input + ' automapped to ' + post_var
									validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(automap)
									relevant = "Select the automapped transcript and click Submit to analyse"
									rel_var = []
									rel_var.append(post_var)
									# Add gene symbols to the link
									cp_rel = copy.copy(rel_var)
									del rel_var[:]
									for accessions in cp_rel:
										error = 'false'
										hgvs_vt = hp.parse_hgvs_variant(str(accessions))
										try:
											tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
										except hgvs.exceptions.HGVSError as e:
											error = str(e)
										if error != 'false':
											accessions = ['', str(hgvs_vt)]
											rel_var.append(accessions)
										else:
										# Get hgnc Gene name from command
											data = va_func.hgnc_rest(path = "/search/prev_symbol/" + tx_id_info[6])
											if data['error'] != 'false':
												reason = 'Cannot currently display the required information:'
												error =  decoded['error']
												validation['warnings'] = validation['warnings'] + ': ' + str(error)
												continue
			
											else:
												# Set the hgnc name correctly
												# If the name is correct no record will be found
												if int(data['record']['response']['numFound']) == 0:
													current = tx_id_info[6]
												else:
													current = data['record']['response']['docs'][0]['symbol']
											accessions = [str(current), str(hgvs_vt)]
											rel_var.append(accessions)
									#Kill current line and append for re-submission
									# Tag the line so that it is not written out
									validation['write'] = 'false'
									# Set the values and append to batch_list
									query = {'quibble' : valstr(hgvs_vt), 'id' : validation['id'], 'warnings' : automap, 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
									batch_list.append(query)
		
							else:
								coding = va_func.coding(variant, hp)
								trans_acc = coding.ac
								# c to Genome coordinates - Map the variant to the genome
								pre_var = va_func.genomic(variant, evm, hp, hdp, primary_assembly)
								# genome back to C coordinates
								post_var = va_func.myevm_g_to_t(hdp, evm, pre_var, trans_acc)
						
								test = hp.parse_hgvs_variant(input)
								if post_var.posedit.pos.start.base != test.posedit.pos.start.base or post_var.posedit.pos.end.base != test.posedit.pos.end.base:
									caution = 'The entered coordinates do not agree with the intron/exon boundaries for the selected transcript:'
									automap = 'Automap has corrected the coordinates to match the intron/exon boundaries for the selected transcript'
									# automapping of variant completed
									automap = str(trapped_input) + ' automapped to ' + str(post_var)
									validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(automap)
									relevant = "Select the automapped transcript and click Submit to analyse"
									rel_var = []
									rel_var.append(post_var)
									# Add gene symbols to the link
									cp_rel = copy.copy(rel_var)
									del rel_var[:]
									for accessions in cp_rel:
										error = 'false'
										hgvs_vt = hp.parse_hgvs_variant(str(accessions))
										try:
											tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
										except hgvs.exceptions.HGVSError as e:
											error = str(e)
										if error != 'false':
											accessions = ['', str(hgvs_vt)]
											rel_var.append(accessions)
										else:
										# Get hgnc Gene name from command
											data = va_func.hgnc_rest(path = "/search/prev_symbol/" + tx_id_info[6])
											if data['error'] != 'false':
												reason = 'Cannot currently display the required information:'
												error =  decoded['error']
												validation['warnings'] = validation['warnings'] + ': ' + str(error)
												continue
				
											else:
												# Set the hgnc name correctly
												# If the name is correct no record will be found
												if int(data['record']['response']['numFound']) == 0:
													current = tx_id_info[6]
												else:
													current = data['record']['response']['docs'][0]['symbol']
											accessions = [str(current), str(hgvs_vt)]
											rel_var.append(accessions)
									#Kill current line and append for re-submission
									# Tag the line so that it is not written out
									validation['write'] = 'false'
									# Set the values and append to batch_list
									query = {'quibble' : valstr(hgvs_vt), 'id' : validation['id'], 'warnings' : automap, 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
									batch_list.append(query)


					# If cck not true
					elif pat_r.search(trapped_input):
						# set input hgvs object
						hgvs_rna_input = hp.parse_hgvs_variant(trapped_input) # Traps the hgvs variant of r. for further use
						inp = str(va_func.hgvs_r_to_c(hgvs_rna_input))
						# Regex
						plus = re.compile("\d\+\d") 	# finds digit + digit
						minus = re.compile("\d\-\d")	# finds digit - digit
						if plus.search(input) or minus.search(input):
							to_g = va_func.genomic(inp, evm, hp, hdp, primary_assembly)
							es = re.compile('error')
							if es.search(str(to_g)):
								if alt_aln_method != 'genebuild':
									error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
									reason = "An error has occurred"
									excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									continue

								else:
									error = "If the following error message does not address the issue and the problem persists please contact admin: " + to_g
									reason = "An error has occurred"
									excep = "%s -- %s -- %s\n" %(time.ctime(), reason, variant)
									validation['warnings'] = validation['warnings'] + ': ' + str(error)
									continue			
	
							else:
								#Set variants pre and post genomic norm
								hgvs_inp = va_func.myevm_g_to_t(hdp, evm,to_g, tx_ac=obj.ac)
								to_g = hn.normalize(to_g)
								hgvs_otp = va_func.myevm_g_to_t(hdp, evm,to_g, tx_ac=obj.ac)					
								tx_ac = ''
						else:	
							# Set variants pre and post RNA norm
							hgvs_inp = hp.parse_hgvs_variant(inp)
							hgvs_otp = hn.normalize(hgvs_inp)
							tx_ac = ''
	
						# Set remaining variables	
						redit = str(hgvs_otp.posedit.edit)
						redit = redit.lower()
						hgvs_otp.posedit.edit = redit
						otp = str(hgvs_otp)
						query = str(hgvs_otp.posedit.pos)
						test = str(hgvs_inp.posedit.pos)
						query = query.replace('T', 'U')
						query = query.replace('ENSU', 'ENST')
						test = test.replace('T', 'U')
						test = test.replace('ENSU', 'ENST')
						output = otp.replace(':c.', ':r.')
						# Apply coordinates test
						if query != test:
							caution = 'The variant description ' + input + ' requires alteration to comply with HGVS variant nomenclature:'
							automap = 'Automap has corrected the variant description'
							# automapping of variant completed
							automap = trapped_input + ' automapped to ' + output
							validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(automap)
							relevant = "Select the automapped transcript and click Submit to analyse"
							rel_var = []
							rel_var.append(output)
							# Add gene symbols to the link
							cp_rel = copy.copy(rel_var)
							del rel_var[:]
							for accessions in cp_rel:
								error = 'false'
								hgvs_vt = hp.parse_hgvs_variant(str(accessions))
								try:
									tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
								except hgvs.exceptions.HGVSError as e:
									error = str(e)
								if error != 'false':
									accessions = ['', str(hgvs_vt)]
									rel_var.append(accessions)
								else:
								# Get hgnc Gene name from command
									data = va_func.hgnc_rest(path = "/search/prev_symbol/" + tx_id_info[6])
									if data['error'] != 'false':
										reason = 'Cannot currently display the required information:'
										error =  decoded['error']
										validation['warnings'] = validation['warnings'] + ': ' + str(error)
										continue					
				
									else:
										# Set the hgnc name correctly
										# If the name is correct no record will be found
										if int(data['record']['response']['numFound']) == 0:
											current = tx_id_info[6]
										else:
											current = data['record']['response']['docs'][0]['symbol']
									accessions = [str(current), str(hgvs_vt)]
									rel_var.append(accessions)
							#Kill current line and append for re-submission
							# Tag the line so that it is not written out
							validation['write'] = 'false'
							# Set the values and append to batch_list
							query = {'quibble' : valstr(hgvs_vt), 'id' : validation['id'], 'warnings' : automap, 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
							batch_list.append(query)

					elif pat_g.search(input):
						pass
	
					else:
						query = hp.parse_hgvs_variant(variant)
						test = hp.parse_hgvs_variant(input)
						if query.posedit.pos != test.posedit.pos:
							caution = 'The variant description ' + input + ' requires alteration to comply with HGVS variant nomenclature:'
							automap = 'Automap has corrected the variant description'
							# automapping of variant completed
							automap = str(test) + ' automapped to ' + str(query)
							validation['warnings'] = validation['warnings'] + ': ' + str(caution) + ': ' + str(automap)
							relevant = "Select the automapped transcript and click Submit to analyse"
							rel_var = []
							rel_var.append(query)
							# Add gene symbols to the link
							cp_rel = copy.copy(rel_var)
							del rel_var[:]
							for accessions in cp_rel:
								error = 'false'
								hgvs_vt = hp.parse_hgvs_variant(str(accessions))
								try:
									tx_id_info = hdp.get_tx_identity_info(str(hgvs_vt.ac))
								except hgvs.exceptions.HGVSError as e:
									error = str(e)
								if error != 'false':
									accessions = ['', str(hgvs_vt)]
									rel_var.append(accessions)
								else:
								# Get hgnc Gene name from command
									data = va_func.hgnc_rest(path = "/search/prev_symbol/" + tx_id_info[6])
									if data['error'] != 'false':
										reason = 'Cannot currently display the required information:'
										error =  decoded['error']
										validation['warnings'] = validation['warnings'] + ': ' + str(error)
										continue
			
									else:
										# Set the hgnc name correctly
										# If the name is correct no record will be found
										if int(data['record']['response']['numFound']) == 0:
											current = tx_id_info[6]
										else:
											current = data['record']['response']['docs'][0]['symbol']
									accessions = [str(current), str(hgvs_vt)]
									rel_var.append(accessions)
							#Kill current line and append for re-submission
							# Tag the line so that it is not written out
							validation['write'] = 'false'
							# Set the values and append to batch_list
							query = {'quibble' : valstr(hgvs_vt), 'id' : validation['id'], 'warnings' : automap, 'description' : '', 'coding' : '', 'coding_g' : '', 'genomic_r' : '', 'genomic_g' : '', 'protein' : '', 'write' : 'true', 'primary_assembly' : primary_assembly, 'order' : ordering}
							batch_list.append(query)


					# VALIDATION of intronic variants
					pre_valid = hp.parse_hgvs_variant(input)
					post_valid = hp.parse_hgvs_variant(variant)
					if valid == 'false':
						error = 'false'
						genomic_validation = str(va_func.genomic(input, evm, hp, hdp, primary_assembly))
						del_end = re.compile('\ddel$')
						delins = re.compile('delins')
						inv = re.compile('inv')
						if valstr(pre_valid) != valstr(post_valid):
							if type != ':g.':
								if caution == '':
									caution = valstr(pre_valid) + ' automapped to ' + valstr(post_valid)
								else:
									pass				
								validation['warnings'] = validation['warnings'] + ': ' + str(caution)
							else:
								pass
						else:
							pass
	
						# Apply validation to intronic variant descriptions (should be valid but make sure)			
						error = va_func.validate(genomic_validation, hp=hp, vr=vr)
						if error == 'false':
							valid = 'true'
						else:

							excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
							validation['warnings'] = validation['warnings'] + ': ' + str(error)
							continue

					if valid == 'true':
						var_tab = 'true'
						cores = "HGVS-compliant variant descriptions"  + warning
	
						# v0.1a1 edit
						if valstr(pre_valid) != valstr(post_valid):
							if type == ':g.':
								if caution == '':
									caution = valstr(pre_valid) + ' automapped to ' + valstr(post_valid)
								else:
									pass				
								validation['warnings'] = validation['warnings'] + ': ' + str(caution)
							else:
								pass
						else:
							pass

						# COLLECT VARIANT DESCRIPTIONS
						##############################
	
						# Coding sequence - BASED ON NORMALIZED VARIANT IF EXONIC
						hgvs_coding = va_func.coding(variant, hp)
						boundary = re.compile('exon-intron boundary')
						spanning = re.compile('exon/intron')
						try:
							hgvs_coding = hn.normalize(hgvs_coding)
						except hgvs.exceptions.HGVSError as e:
							error = str(e)
	
						if boundary.search(str(error)) or spanning.search(str(error)):
							try:
								hgvs_coding = evm._maybe_normalize(hgvs_coding)
							except hgvs.exceptions.HGVSError as error:
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
						else:
							pass
	
						# END
				
						coding = valstr(hgvs_coding)
						# print 'hgvs_coding'
						# print hgvs_coding

	
						# RNA sequence
						hgvs_rna = copy.deepcopy(hgvs_coding)
						hgvs_rna = va_func.hgvs_c_to_r(hgvs_rna)
						rna = str(hgvs_rna)
	
						# Genomic sequence
						try:
							hgvs_genomic = va_func.myevm_t_to_g(hgvs_coding, evm, hdp, primary_assembly)
						except hgvs.exceptions.HGVSError as e:
							error = str(e)
						if error == 'false':
							hgvs_genomic = va_func.genomic(variant, evm, hp, hdp, primary_assembly)
						elif boundary.search(str(error)) or spanning.search(str(error)):
							hgvs_genomic = va_func.myevm_t_to_g(hgvs_coding, evm, hdp, primary_assembly)
						else:
							ens = re.compile('ENS')
							ref = re.compile('NM')
							if alt_aln_method == 'genebuild' and ref.search(variant):
								error = str(error)
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
							elif alt_aln_method != 'genebuild' and ens.search(variant):
								error = str(error)	
								validation['warnings'] = validation['warnings'] + ': ' + str(error)
								continue
							else:
								pass
		
						# genomic_possibilities
						# 1. take the simple 3 pr normalized hgvs_genomic
						# 2. Lock in hgvs_genomic at its most 5 prime position wrt genome
						hgvs_genomic_possibilities = []
						rn_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic)
						hgvs_genomic_possibilities.append(rn_hgvs_genomic)
						if orientation != -1:
							try:
								chromosome_normalized_hgvs_coding = reverse_normalize.normalize(hgvs_coding)
							except hgvs.exceptions.HGVSUnsupportedOperationError as e:
								error = str(e)
								if re.match('Normalization of intronic variants is not supported', error):
									chromosome_normalized_hgvs_coding = hgvs_coding	
						else:
							try:	
								chromosome_normalized_hgvs_coding = hn.normalize(hgvs_coding)	
							except hgvs.exceptions.HGVSUnsupportedOperationError as e:
								error = str(e)
								if re.match('Normalization of intronic variants is not supported', error):
									chromosome_normalized_hgvs_coding = hgvs_coding
			
						# print 'chromosome_normalized_hgvs_coding'
						# print chromosome_normalized_hgvs_coding

						most_3pr_hgvs_genomic = va_func.myevm_t_to_g(chromosome_normalized_hgvs_coding, no_norm_evm, hdp, primary_assembly)
						hgvs_genomic_possibilities.append(most_3pr_hgvs_genomic)

						# Set variables for problem specific warnings
						gapped_alignment_warning = ''
						corrective_action_taken = ''
						gapped_transcripts = ''
						auto_info = ''
			
						# Mark as not disparity detected
						disparity_deletion_in = ['false', 'false']
						# print hgvs_genomic_possibilities
						
						# Loop through to see if a gap can be located
						for possibility in hgvs_genomic_possibilities:
				
							# Use VCF generation code to push hgvs_genomic as for 5 prime as possible to uncover gaps
							hgvs_genomic_variant = possibility
							stored_hgvs_genomic_variant = copy.deepcopy(hgvs_genomic_variant)					

							# Reverse normalize hgvs_genomic_variant: NOTE will replace ref
							try:
								reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic_variant)
							except hgvs.exceptions.HGVSError as e:
								# Strange error caused by gap in genomic
								error = str(e)
								if re.search('base start position must be <= end position', error):
									if hgvs_genomic.posedit.edit.type == 'delins':
										start = hgvs_genomic.posedit.pos.start.base
										end = hgvs_genomic.posedit.pos.end.base
										lhb = sf.fetch_seq(str(hgvs_genomic.ac),end-1,end)
										rhb = sf.fetch_seq(str(hgvs_genomic.ac),start-1,start)						
										hgvs_genomic.posedit.edit.ref = lhb + rhb
										hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
										hgvs_genomic.posedit.pos.start.base = end
										hgvs_genomic.posedit.pos.end.base = start
										reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic)
									if hgvs_genomic.posedit.edit.type == 'del':
										start = hgvs_genomic.posedit.pos.start.base
										end = hgvs_genomic.posedit.pos.end.base
										lhb = sf.fetch_seq(str(hgvs_genomic.ac),end-1,end)
										rhb = sf.fetch_seq(str(hgvs_genomic.ac),start-1,start)						
										hgvs_genomic.posedit.edit.ref = lhb + rhb
										hgvs_genomic.posedit.edit.alt = lhb + rhb
										hgvs_genomic.posedit.pos.start.base = end
										hgvs_genomic.posedit.pos.end.base = start
										reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic)
								if re.search('insertion length must be 1', error):
									if hgvs_genomic.posedit.edit.type == 'ins':
										start = hgvs_genomic.posedit.pos.start.base
										end = hgvs_genomic.posedit.pos.end.base
										ref_bases = sf.fetch_seq(str(hgvs_genomic.ac),start-1,end)
										lhb = sf.fetch_seq(str(hgvs_genomic.ac),start-1,start)
										rhb = sf.fetch_seq(str(hgvs_genomic.ac),start,end)						
										hgvs_genomic.posedit.edit.ref = lhb + rhb
										hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
										reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic)						

							hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
							# Store a copy for later use
							stored_hgvs_genomic_5pr = copy.deepcopy(hgvs_genomic_5pr)
				
							# Create VCF
							vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic)
							chr = vcf_dict['chr']	
							pos = vcf_dict['pos']				
							ref = vcf_dict['ref']
							alt = vcf_dict['alt']
				
							# Look for exonic gaps within transcript or chromosome
							no_normalized_c = 'false' # Mark true to not produce an additional normalization of c.

							# Generate an end position
							end = str(int(pos) + len(ref) -1)
							pos = str(pos)
	
							# Store a not real deletion insertion to test for gapping
							stored_hgvs_not_delins = hp.parse_hgvs_variant(str(hgvs_genomic_5pr.ac) + ':' +  hgvs_genomic_5pr.type + '.' +  pos + '_' + end + 'del' + ref + 'ins' + alt)
							v = [chr, pos,ref,alt]

							# print 'stored_hgvs_not_delins = ' + str(stored_hgvs_not_delins)
				
				
							# print 'hash'
							v = [chr, pos,ref,alt]
				
				
							# print '-'.join(v)

							# Save a copy of current hgvs_coding
							saved_hgvs_coding = no_norm_evm.g_to_t(stored_hgvs_not_delins, hgvs_coding.ac)
		
							# TO BATCH AND API
							# Detect intronic variation using normalization
							intronic_variant = 'false'
							# Look for normalized variant options that do not match hgvs_coding
							if orientation == -1:
								# position genomic at its most 5 prime position
								try:
									query_genomic = reverse_normalize.normalize(hgvs_genomic)
								except:
									query_genomic = hgvs_genomic	
								# Map to the transcript and test for movement
								try:
									hgvs_seek_var = evm.g_to_t(query_genomic, hgvs_coding.ac)
								except hgvs.exceptions.HGVSError as e:
									hgvs_seek_var = saved_hgvs_coding
								else:		
									seek_var = valstr(hgvs_seek_var)
									seek_ac = str(hgvs_seek_var.ac)
								if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (hgvs_coding.posedit.pos.start.base +  hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (hgvs_coding.posedit.pos.end.base +  hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
									pass
								else:
									hgvs_seek_var = saved_hgvs_coding					

							elif orientation != -1:
								# position genomic at its most 3 prime position
								try:
									query_genomic = hn.normalize(hgvs_genomic)
								except:
									query_genomic = hgvs_genomic	
								# Map to the transcript ant test for movement
								try:
									hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
								except hgvs.exceptions.HGVSError as e:
									hgvs_seek_var = saved_hgvs_coding
								else:		
									seek_var = valstr(hgvs_seek_var)
									seek_ac = str(hgvs_seek_var.ac)
								if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (hgvs_coding.posedit.pos.start.base +  hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (hgvs_coding.posedit.pos.end.base +  hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
									pass
								else:
									hgvs_seek_var = saved_hgvs_coding		

							try:	
								intron_test = hn.normalize(hgvs_seek_var)	
							except hgvs.exceptions.HGVSUnsupportedOperationError as e:
								error = str(e)
								if re.match('Normalization of intronic variants is not supported', error):
									# Double check to see whether the variant is actually intronic?
									for exon in ori:
										genomic_start = int(exon['alt_start_i'])
										genomic_end = int(exon['alt_end_i'])
										if (hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
											intronic_variant = 'false'
											break
										else:							
											intronic_variant = 'true'

							if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-', str(hgvs_seek_var.posedit.pos)):
								# Double check to see whether the variant is actually intronic?
								for exon in ori:
									genomic_start = int(exon['alt_start_i'])
									genomic_end = int(exon['alt_end_i'])
									if (hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
										intronic_variant = 'false'
										break							
									else:							
										intronic_variant = 'true'

							if intronic_variant != 'true':
								# Flag RefSeqGene for ammendment
								amend_RefSeqGene = 'false'
								# Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
								if stored_hgvs_not_delins != '':
									# Refresh hgvs_not_delins from stored_hgvs_not_delins
									hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
									# This test will only occur in dup of single base, insertion or substitution
									if not re.search('_', str(hgvs_not_delins.posedit.pos)):
										# For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
										plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
										plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
										plussed_hgvs_not_delins.posedit.edit.ref = ''
										transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins, str(saved_hgvs_coding.ac))
										if ((transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
											if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
												start = hgvs_not_delins.posedit.pos.start.base - 1
												end = hgvs_not_delins.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
												hgvs_not_delins.posedit.edit.ref = ref_bases
												hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]
											elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
											elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and not re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
												start = hgvs_not_delins.posedit.pos.start.base - 1
												end = hgvs_not_delins.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
												hgvs_not_delins.posedit.edit.ref = ref_bases
												hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]						
										else:
											if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
												start = hgvs_not_delins.posedit.pos.start.base - 1
												end = hgvs_not_delins.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
												hgvs_not_delins.posedit.edit.ref = ref_bases
												hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] # + ref_bases[1:]
											elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
											elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and not re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
												hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
												start = hgvs_not_delins.posedit.pos.start.base - 1
												end = hgvs_not_delins.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
												hgvs_not_delins.posedit.edit.ref = ref_bases
												hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] # + ref_bases[1:]							
									else:
										pass						
									tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
									# Create normalized version of tx_hgvs_not_delins
									rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
									# Check for +1 base and adjust
									if re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
										# Remove offsetting to span the gap
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										try:
											rn_tx_hgvs_not_delins.posedit.edit.alt = '' 
										except:
											pass
					
									elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
										# move tx end base to next available non-offset base
										rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
										rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										if re.match('NM_', str(rn_tx_hgvs_not_delins)):
											test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
										else:
											test_tx_var = rn_tx_hgvs_not_delins	
										# re-make genomic and tx
										hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
										rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
										# tx_hgvs_not_delins = rn_tx_hgvs_not_delins
									elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
										# move tx start base to previous available non-offset base
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										if re.match('NM_', str(rn_tx_hgvs_not_delins)):
											test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
										else:
											test_tx_var = rn_tx_hgvs_not_delins	
										# re-make genomic and tx
										hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
										rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0							
									else:
										pass

									# Check for -ve base and adjust
									if re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
										# Remove offsetting to span the gap
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										try:
											rn_tx_hgvs_not_delins.posedit.edit.alt = '' 
										except:
											pass
									elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
										rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
										# Delete the ref
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										# Add the additional base to the ALT
										start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
										end = rn_tx_hgvs_not_delins.posedit.pos.end.base
										ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac),start,end)
										rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
										if re.match('NM_', str(rn_tx_hgvs_not_delins)):
											test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
										else:
											test_tx_var = rn_tx_hgvs_not_delins	
										# re-make genomic and tx
										hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
										rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
									elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
										# move tx start base to previous available non-offset base
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
										rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
										rn_tx_hgvs_not_delins.posedit.edit.ref = ''
										if re.match('NM_', str(rn_tx_hgvs_not_delins)):
											test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
										else:
											test_tx_var = rn_tx_hgvs_not_delins	
										# re-make genomic and tx
										hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
										rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
										rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0							
									else:
										pass

									# Logic
									if len(hgvs_not_delins.posedit.edit.ref) < len(rn_tx_hgvs_not_delins.posedit.edit.ref):
										gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(hgvs_not_delins.posedit.edit.ref)
										disparity_deletion_in = ['chromosome', gap_length]
									elif len(hgvs_not_delins.posedit.edit.ref) > len(rn_tx_hgvs_not_delins.posedit.edit.ref):
										gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(rn_tx_hgvs_not_delins.posedit.edit.ref)
										disparity_deletion_in = ['transcript', gap_length]
									else:
										pass

								# print 'At hgvs_genomic'
								# print disparity_deletion_in	
																		
								amend_RefSeqGene = 'false'
								# Recreate hgvs_genomic
								if disparity_deletion_in[0] == 'transcript':
									hgvs_genomic = hgvs_not_delins 	
						
								# GAP IN THE TRANSCRIPT	DISPARITY DETECTED	
								if disparity_deletion_in[0] == 'transcript':			
						
									amend_RefSeqGene = 'true'
									# ANY VARIANT WHOLLY WITHIN THE GAP
									if (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start))) and (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end))):
										auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
										non_valid_caution = 'true'
										if str(tx_hgvs_not_delins.posedit.pos.start) == str(tx_hgvs_not_delins.posedit.pos.end):
											if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
												tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
												tx_hgvs_not_delins.posedit.pos.end.offset = 0
											elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
												tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
												tx_hgvs_not_delins.posedit.pos.start.offset = 0									
										if re.search('\=', str(tx_hgvs_not_delins.posedit)):
											middle = tx_hgvs_not_delins.posedit.edit.alt
										elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)):
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
										elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
											middle = tx_hgvs_not_delins.posedit.edit.alt#[1:]
										elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
											middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
											tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
										tx_hgvs_not_delins.posedit.pos.start.offset = 0
										tx_hgvs_not_delins.posedit.pos.end.offset = 0							
										beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base)
										end = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
										tx_hgvs_not_delins.posedit.edit.ref = beginning + end
										tx_hgvs_not_delins.posedit.edit.alt = beginning + middle + end
										tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
										hgvs_refreshed_variant = tx_hgvs_not_delins
										# Alignment position
										for_location_c = copy.deepcopy(hgvs_refreshed_variant)
										if re.match('NM_', str(for_location_c)):
											for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										gps = for_location_c.posedit.pos.start.base -1
										gpe = for_location_c.posedit.pos.start.base
										gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'															
										# Warn update
										auto_info = auto_info + '%s' %(gap_position)
									else:
										if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)):
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.start.base) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
											non_valid_caution = 'true'
											tx_hgvs_not_delins.posedit.pos.start.offset = 0
											beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base) 
											middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
											end = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base, tx_hgvs_not_delins.posedit.pos.end.base)
											tx_hgvs_not_delins.posedit.edit.ref = beginning + end
											tx_hgvs_not_delins.posedit.edit.alt = beginning + middle								
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
											hgvs_refreshed_variant = tx_hgvs_not_delins
											# Alignment position
											for_location_c = copy.deepcopy(hgvs_refreshed_variant)
											if re.match('NM_', str(for_location_c)):
												for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											gps = for_location_c.posedit.pos.start.base -1
											gpe = for_location_c.posedit.pos.start.base
											gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
											# Warn update
											auto_info = auto_info + '%s' %(gap_position)
										elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
											non_valid_caution = 'true'
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
											tx_hgvs_not_delins.posedit.pos.end.offset = 0
											beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.end.base - 1)
											middle = tx_hgvs_not_delins.posedit.edit.alt
											end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
											tx_hgvs_not_delins.posedit.edit.ref = beginning + end
											tx_hgvs_not_delins.posedit.edit.alt = middle + end
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
											hgvs_refreshed_variant = tx_hgvs_not_delins				
											# Alignment position
											for_location_c = copy.deepcopy(hgvs_refreshed_variant)
											if re.match('NM_', str(for_location_c)):
												for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											gps = for_location_c.posedit.pos.start.base -1
											gpe = for_location_c.posedit.pos.start.base
											gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'			
											# Warn update
											auto_info = auto_info + '%s' %(gap_position)
										elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)):								
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.start.base) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
											non_valid_caution = 'true'
											tx_hgvs_not_delins.posedit.pos.start.offset = 0
											tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
											beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base)
											middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
											end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base, tx_hgvs_not_delins.posedit.pos.end.base)
											if not re.search('\=', str(tx_hgvs_not_delins.posedit)):
												tx_hgvs_not_delins.posedit.edit.alt = beginning + middle
											else:
												tx_hgvs_not_delins.posedit.edit.alt = beginning + middle + end
											tx_hgvs_not_delins.posedit.edit.ref = beginning + end
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
											hgvs_refreshed_variant = tx_hgvs_not_delins								
											# Alignment position
											for_location_c = copy.deepcopy(hgvs_refreshed_variant)
											if re.match('NM_', str(for_location_c)):
												for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											gps = for_location_c.posedit.pos.start.base -1
											gpe = for_location_c.posedit.pos.start.base
											gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
											# Warn update
											auto_info = auto_info + '%s' %(gap_position)
										elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
											non_valid_caution = 'true'
											tx_hgvs_not_delins.posedit.pos.end.offset = 0
											beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.end.base - 1)
											if len(tx_hgvs_not_delins.posedit.edit.ref) < len(tx_hgvs_not_delins.posedit.edit.alt):
												middle = tx_hgvs_not_delins.posedit.edit.alt.replace(tx_hgvs_not_delins.posedit.edit.ref, '' , 1)
												middle = beginning[0] + middle
											elif len(tx_hgvs_not_delins.posedit.edit.ref) > len(tx_hgvs_not_delins.posedit.edit.alt):
												middle = tx_hgvs_not_delins.posedit.edit.alt
											elif len(tx_hgvs_not_delins.posedit.edit.ref) == len(tx_hgvs_not_delins.posedit.edit.alt):
												middle = tx_hgvs_not_delins.posedit.edit.alt																
											end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
											tx_hgvs_not_delins.posedit.edit.ref = beginning + end
											tx_hgvs_not_delins.posedit.edit.alt = middle + end
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
											hgvs_refreshed_variant = tx_hgvs_not_delins								
											# Alignment position
											for_location_c = copy.deepcopy(hgvs_refreshed_variant)
											if re.match('NM_', str(for_location_c)):
												for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
											gps = for_location_c.posedit.pos.start.base -1
											gpe = for_location_c.posedit.pos.start.base
											gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'							
											# Warn update
											auto_info = auto_info + '%s' %(gap_position)
										else:
											auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos) + ' aligns across a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac) + '\n'
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
											hgvs_refreshed_variant = tx_hgvs_not_delins
						
								# GAP IN THE CHROMOSOME		
								elif disparity_deletion_in[0] == 'chromosome':
									amend_RefSeqGene = 'true'
									hgvs_refreshed_variant = chromosome_normalized_hgvs_coding
									# print hgvs_refreshed_variant
								else:
									# Keep the same by re-setting rel_var
									hgvs_refreshed_variant = hgvs_coding
									amend_RefSeqGene = 'false'

								# Edit the output
								if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(hgvs_refreshed_variant.type)):
									hgvs_refreshed_variant = no_norm_evm.n_to_c(hgvs_refreshed_variant)
								else:
									pass
								# Quick check to make sure the coding variant has not changed
								to_test = hn.normalize(hgvs_refreshed_variant)
								if str(to_test.posedit.edit) != str(hgvs_coding.posedit.edit):
									# Try the next available genomic option
									continue
								# Update hgvs_genomic
								hgvs_genomic = va_func.myevm_t_to_g(hgvs_refreshed_variant, no_norm_evm, hdp, primary_assembly)
								# print 'un-normalized genomic = ' + str(hgvs_genomic)				
								
								# Genomic gap warning
								if disparity_deletion_in[0] == 'chromosome':
									auto_info = auto_info + "Prior to 3' normalization, Transcript variant " + str(hgvs_refreshed_variant) + ' aligns across a ' + str(disparity_deletion_in[1]) + '-bp gap in genomic sequence ' + str(stored_hgvs_not_delins.ac) + ' between positions ' + str(hgvs_genomic.posedit.pos.start.base) + '_' + str(hgvs_genomic.posedit.pos.end.base) + '\n'
							
							
							# If it is intronic, these vairables will not have been set
							else:
								amend_RefSeqGene = 'false'
								no_normalized_c = 'false'	
				
							# Break if gap has been detected
							if disparity_deletion_in[0] != 'false':
								break

						# Warn user about gapping								
						if auto_info != '':
							auto_info = str(auto_info.replace('\n', ': '))
							validation['warnings'] = validation['warnings'] + ': ' + str(auto_info)
													
						# Normailse hgvs_genomic
						try:
							hgvs_genomic = hn.normalize(hgvs_genomic)
						except hgvs.exceptions.HGVSError as e:
							# Strange error caused by gap in genomic
							error = str(e)

							if re.search('base start position must be <= end position', error) and disparity_deletion_in[0] == 'chromosome':
								if hgvs_genomic.posedit.edit.type == 'delins':
									start = hgvs_genomic.posedit.pos.start.base
									end = hgvs_genomic.posedit.pos.end.base
									lhb = sf.fetch_seq(str(hgvs_genomic.ac),end-1,end)
									rhb = sf.fetch_seq(str(hgvs_genomic.ac),start-1,start)
									hgvs_genomic.posedit.edit.ref = lhb + rhb
									hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
									hgvs_genomic.posedit.pos.start.base = end
									hgvs_genomic.posedit.pos.end.base = start
									hgvs_genomic = hn.normalize(hgvs_genomic)
								if hgvs_genomic.posedit.edit.type == 'del':
									start = hgvs_genomic.posedit.pos.start.base
									end = hgvs_genomic.posedit.pos.end.base
									lhb = sf.fetch_seq(str(hgvs_genomic.ac),end-1,end)
									rhb = sf.fetch_seq(str(hgvs_genomic.ac),start-1,start)
									hgvs_genomic.posedit.edit.ref = lhb + rhb
									hgvs_genomic.posedit.edit.alt = lhb + rhb
									hgvs_genomic.posedit.pos.start.base = end
									hgvs_genomic.posedit.pos.end.base = start
									hgvs_genomic = hn.normalize(hgvs_genomic)
						genomic = valstr(hgvs_genomic)
			
						# print 'normalized_hgvs_genomic = ' + str(hgvs_genomic) + '\n'

						# Create pseudo VCF based on amended hgvs_genomic
						hgvs_genomic_variant = hgvs_genomic					
						# Reverse normalize hgvs_genomic_variant: NOTE will replace ref
						reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic_variant)

						hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
	
						# Create vcf
						vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic)
						chr = vcf_dict['chr']	
						pos = vcf_dict['pos']				
						ref = vcf_dict['ref']
						alt = vcf_dict['alt']
			
						# Create a VCF call
						vcf_component_list = [str(chr), str(pos), str(ref), (alt)]
						vcf_genomic = '-'.join(vcf_component_list)

						# DO NOT DELETE
						# Generate an end position
						end = str(int(pos) + len(ref) -1)
						pos = str(pos)
		
						# DO NOT DELETE
						stored_hgvs_not_delins = hp.parse_hgvs_variant(str(hgvs_genomic_5pr.ac) + ':' +  hgvs_genomic_5pr.type + '.' +  pos + '_' + end + 'del' + ref + 'ins' + alt)

						
						# Apply gap code to re-format hgvs_coding
						# Store the current hgvs:c. description
						saved_hgvs_coding = copy.deepcopy(hgvs_coding)

						# Get orientation of the gene wrt genome and a list of exons mapped to the genome
						ori = va_func.tx_exons(tx_ac=saved_hgvs_coding.ac, alt_ac=hgvs_genomic_5pr.ac, alt_aln_method=alt_aln_method, hdp=hdp)
						orientation = int(ori[0]['alt_strand'])
	
						# Look for normalized variant options that do not match hgvs_coding
						hgvs_genomic = copy.deepcopy(hgvs_genomic_variant)
						if orientation == -1:
							# position genomic at its most 5 prime position
							try:
								query_genomic = reverse_normalize.normalize(hgvs_genomic)
							except:
								query_genomic = hgvs_genomic	
							# Map to the transcript ant test for movement
							try:
								hgvs_seek_var = evm.g_to_t(query_genomic, hgvs_coding.ac)
							except hgvs.exceptions.HGVSError as e:
								hgvs_seek_var = saved_hgvs_coding
							else:		
								seek_var = valstr(hgvs_seek_var)
								seek_ac = str(hgvs_seek_var.ac)
							if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (hgvs_coding.posedit.pos.start.base +  hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (hgvs_coding.posedit.pos.end.base +  hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
								pass
							else:
								hgvs_seek_var = saved_hgvs_coding					
 
						elif orientation != -1:
							# position genomic at its most 3 prime position
							try:
								query_genomic = hn.normalize(hgvs_genomic)
							except:
								query_genomic = hgvs_genomic	
							# Map to the transcript ant test for movement
							try:
								hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
							except hgvs.exceptions.HGVSError as e:
								hgvs_seek_var = saved_hgvs_coding
							else:		
								seek_var = valstr(hgvs_seek_var)
								seek_ac = str(hgvs_seek_var.ac)
							if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (saved_hgvs_coding.posedit.pos.start.base +  saved_hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (saved_hgvs_coding.posedit.pos.end.base +  saved_hgvs_coding.posedit.pos.end.offset):
								pass
							else:
								hgvs_seek_var = saved_hgvs_coding	
	
						# is it in an exon?
						is_it_in_an_exon = 'no'
						for exon in ori:
							genomic_start = int(exon['alt_start_i'])
							genomic_end = int(exon['alt_end_i'])
							# Take from stored copy
							hgvs_genomic_5pr = copy.deepcopy(stored_hgvs_genomic_5pr)
							if (hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
								is_it_in_an_exon = 'yes'	
						if is_it_in_an_exon == 'yes':
							# map form reverse normalized g. to c. 
							hgvs_from_5n_g = no_norm_evm.g_to_t(hgvs_genomic_5pr, saved_hgvs_coding.ac)
				
							# Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
							disparity_deletion_in = ['false', 'false']
							if stored_hgvs_not_delins != '':
								# Refresh hgvs_not_delins from stored_hgvs_not_delins
								hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
								# This test will only occur in dup of single base, insertion or substitution
								if not re.search('_', str(hgvs_not_delins.posedit.pos)):
									# For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
									plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
									plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
									plussed_hgvs_not_delins.posedit.edit.ref = ''
									transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins, str(saved_hgvs_coding.ac))
									if ((transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
										if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
											hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
											start = hgvs_not_delins.posedit.pos.start.base - 1
											end = hgvs_not_delins.posedit.pos.end.base
											ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
											hgvs_not_delins.posedit.edit.ref = ref_bases
											hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]
										elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
											hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
										elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and not re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
											hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
											start = hgvs_not_delins.posedit.pos.start.base - 1
											end = hgvs_not_delins.posedit.pos.end.base
											ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
											hgvs_not_delins.posedit.edit.ref = ref_bases
											hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]						
									else:
										if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
											hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
											start = hgvs_not_delins.posedit.pos.start.base - 1
											end = hgvs_not_delins.posedit.pos.end.base
											ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
											hgvs_not_delins.posedit.edit.ref = ref_bases
											hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]
										elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
											hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
										elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and not re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
											hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
											start = hgvs_not_delins.posedit.pos.start.base - 1
											end = hgvs_not_delins.posedit.pos.end.base
											ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
											hgvs_not_delins.posedit.edit.ref = ref_bases
											hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] # + ref_bases[1:]						
								else:
									pass						
								tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
								# Create normalized version of tx_hgvs_not_delins
								rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
								# Check for +ve base and adjust
								if re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
									# Remove offsetting to span the gap
									rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
									rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
									rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
									rn_tx_hgvs_not_delins.posedit.edit.ref = ''
									try:
										rn_tx_hgvs_not_delins.posedit.edit.alt = '' 
									except:
										pass
					
								elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
									# move tx end base to next available non-offset base
									rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
									rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
									rn_tx_hgvs_not_delins.posedit.edit.ref = ''
									if re.match('NM_', str(rn_tx_hgvs_not_delins)):
										test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
									else:
										test_tx_var = rn_tx_hgvs_not_delins	
									# re-make genomic and tx
									hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
									rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
								elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
									rn_tx_hgvs_not_delins.posedit.edit.ref = ''
									# move tx start base to previous available non-offset base
									rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
									if re.match('NM_', str(rn_tx_hgvs_not_delins)):
										test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
									else:
										test_tx_var = rn_tx_hgvs_not_delins	
									# re-make genomic and tx
									hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
									rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
									rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0							
								else:
									pass

								# Check for -ve base and adjust
								if re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
									# Remove offsetting to span the gap
									rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
									rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
									rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
									rn_tx_hgvs_not_delins.posedit.edit.ref = ''
									try:
										rn_tx_hgvs_not_delins.posedit.edit.alt = '' 
									except:
										pass
								elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
									# move tx end base back to next available non-offset base
									# rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base - 1
									rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
									# Delete the ref
									rn_tx_hgvs_not_delins.posedit.edit.ref = ''
									# Add the additional base to the ALT
									start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
									end = rn_tx_hgvs_not_delins.posedit.pos.end.base
									ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac),start,end)
									rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
									if re.match('NM_', str(rn_tx_hgvs_not_delins)):
										test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
									else:
										test_tx_var = rn_tx_hgvs_not_delins	
									# re-make genomic and tx
									hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
									rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
								elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
									# move tx start base to previous available non-offset base
									rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
									rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
									rn_tx_hgvs_not_delins.posedit.edit.ref = ''
									if re.match('NM_', str(rn_tx_hgvs_not_delins)):
										test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
									else:
										test_tx_var = rn_tx_hgvs_not_delins	
									# re-make genomic and tx
									hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
									rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
									rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0							
								else:
									pass

								# Logic
								if len(hgvs_not_delins.posedit.edit.ref) < len(rn_tx_hgvs_not_delins.posedit.edit.ref):
									gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(hgvs_not_delins.posedit.edit.ref)
									disparity_deletion_in = ['chromosome', gap_length]
								elif len(hgvs_not_delins.posedit.edit.ref) > len(rn_tx_hgvs_not_delins.posedit.edit.ref):
									gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(rn_tx_hgvs_not_delins.posedit.edit.ref)
									disparity_deletion_in = ['transcript', gap_length]
								else:
									pass
																					
							# GAP IN THE TRANSCRIPT	DISPARITY DETECTED	
							if disparity_deletion_in[0] == 'transcript':			
								gap_position = ''
								gapped_alignment_warning = str(hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + primary_assembly

								# ANY VARIANT WHOLLY WITHIN THE GAP
								if (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start))) and (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end))):
									non_valid_caution = 'true'
									auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos) + ' aligns to a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)						
									gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
									if str(tx_hgvs_not_delins.posedit.pos.start) == str(tx_hgvs_not_delins.posedit.pos.end):
										if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
											tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
											tx_hgvs_not_delins.posedit.pos.end.offset = 0
										elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
											tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
											tx_hgvs_not_delins.posedit.pos.start.offset = 0									
									if re.search('\=', str(tx_hgvs_not_delins.posedit)):
										middle = tx_hgvs_not_delins.posedit.edit.alt
									elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)):
										tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
									elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
										middle = tx_hgvs_not_delins.posedit.edit.alt#[1:]
									elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
										middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
										tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
									tx_hgvs_not_delins.posedit.pos.start.offset = 0
									tx_hgvs_not_delins.posedit.pos.end.offset = 0							
									beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base)
									end = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
									tx_hgvs_not_delins.posedit.edit.ref = beginning + end
									tx_hgvs_not_delins.posedit.edit.alt = beginning + middle + end
									# map to the genome and back						
									tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
									if re.match('NM_', str(tx_hgvs_not_delins.ac)):
										c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
									else:
										c = tx_hgvs_not_delins	
									g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
									c = no_norm_evm.g_to_t(g, str(c.ac))
									if re.match('NM_', str(tx_hgvs_not_delins.ac)):
										tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
									else:
										tx_hgvs_not_delins = c	
									hgvs_refreshed_variant = tx_hgvs_not_delins								
									# Alignment position
									for_location_c = copy.deepcopy(hgvs_refreshed_variant)
									if re.match('NM_', str(for_location_c)):
										for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
									gps = for_location_c.posedit.pos.start.base
									gpe = for_location_c.posedit.pos.start.base + 1
									gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
									auto_info = auto_info + '%s' %(gap_position)

								else:
									if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)):
										auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.start.base) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
										non_valid_caution = 'true'
										tx_hgvs_not_delins.posedit.pos.start.offset = 0
										beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base) 
										middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
										end = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base, tx_hgvs_not_delins.posedit.pos.end.base)
										tx_hgvs_not_delins.posedit.edit.ref = beginning + end
										tx_hgvs_not_delins.posedit.edit.alt = beginning + middle
										# map to the genome and back						
										tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										else:
											c = tx_hgvs_not_delins	
										g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
										c = no_norm_evm.g_to_t(g, str(c.ac))
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
										else:
											tx_hgvs_not_delins = c									
										hgvs_refreshed_variant = tx_hgvs_not_delins
										# Alignment position
										for_location_c = copy.deepcopy(hgvs_refreshed_variant)
										if re.match('NM_', str(for_location_c)):
											for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										gps = for_location_c.posedit.pos.start.base -1
										gpe = for_location_c.posedit.pos.start.base
										gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'

									elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
										auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)						
										gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
										non_valid_caution = 'true'
										tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
										tx_hgvs_not_delins.posedit.pos.end.offset = 0
										beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.end.base - 1)
										middle = tx_hgvs_not_delins.posedit.edit.alt
										end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
										tx_hgvs_not_delins.posedit.edit.ref = beginning + end
										tx_hgvs_not_delins.posedit.edit.alt = middle + end
										# map to the genome and back						
										tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										else:
											c = tx_hgvs_not_delins	
										g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
										c = no_norm_evm.g_to_t(g, str(c.ac))
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
										else:
											tx_hgvs_not_delins = c	
										hgvs_refreshed_variant = tx_hgvs_not_delins								
										# Alignment position
										for_location_c = copy.deepcopy(hgvs_refreshed_variant)
										if re.match('NM_', str(for_location_c)):
											for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										gps = for_location_c.posedit.pos.end.base - 1
										gpe = for_location_c.posedit.pos.end.base
										gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'

									elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)):
										auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.start.base) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
										non_valid_caution = 'true'
										tx_hgvs_not_delins.posedit.pos.start.offset = 0
										tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
										beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base)
										middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
										end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base, tx_hgvs_not_delins.posedit.pos.end.base)
										if not re.search('\=', str(tx_hgvs_not_delins.posedit)):
											tx_hgvs_not_delins.posedit.edit.alt = beginning + middle
										else:
											tx_hgvs_not_delins.posedit.edit.alt = beginning + middle + end
										tx_hgvs_not_delins.posedit.edit.ref = beginning + end
										# map to the genome and back						
										tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										else:
											c = tx_hgvs_not_delins	
										g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
										c = no_norm_evm.g_to_t(g, str(c.ac))
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
										else:
											tx_hgvs_not_delins = c	
										hgvs_refreshed_variant = tx_hgvs_not_delins									
										# Alignment position
										for_location_c = copy.deepcopy(hgvs_refreshed_variant)
										if re.match('NM_', str(for_location_c)):
											for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										gps = for_location_c.posedit.pos.start.base -1
										gpe = for_location_c.posedit.pos.start.base
										gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'

									elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
										auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)						
										gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
										non_valid_caution = 'true'
										tx_hgvs_not_delins.posedit.pos.end.offset = 0
										beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.end.base - 1)
										if len(tx_hgvs_not_delins.posedit.edit.ref) < len(tx_hgvs_not_delins.posedit.edit.alt):
											middle = tx_hgvs_not_delins.posedit.edit.alt.replace(tx_hgvs_not_delins.posedit.edit.ref, '' , 1)
											middle = beginning[0] + middle
										elif len(tx_hgvs_not_delins.posedit.edit.ref) > len(tx_hgvs_not_delins.posedit.edit.alt):
											middle = tx_hgvs_not_delins.posedit.edit.alt
										elif len(tx_hgvs_not_delins.posedit.edit.ref) == len(tx_hgvs_not_delins.posedit.edit.alt):
											middle = tx_hgvs_not_delins.posedit.edit.alt																
										end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
										tx_hgvs_not_delins.posedit.edit.ref = beginning + end
										tx_hgvs_not_delins.posedit.edit.alt = middle + end
										# map to the genome and back						
										tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.start.base + len(tx_hgvs_not_delins.posedit.edit.ref) - 1
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										else:
											c = tx_hgvs_not_delins	
										g = va_func.myevm_t_to_g(c, no_norm_evm, hdp, primary_assembly)
										c = no_norm_evm.g_to_t(g, str(c.ac))
										if re.match('NM_', str(tx_hgvs_not_delins.ac)):
											tx_hgvs_not_delins = no_norm_evm.c_to_n(c)
										else:
											tx_hgvs_not_delins = c	
										hgvs_refreshed_variant = tx_hgvs_not_delins								
										# Alignment position
										for_location_c = copy.deepcopy(hgvs_refreshed_variant)
										if re.match('NM_', str(for_location_c)):
											for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
										gps = for_location_c.posedit.pos.end.base - 1
										gpe = for_location_c.posedit.pos.end.base
										gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
									else:
										auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos) + ' aligns across a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac) + '\n'
										hgvs_refreshed_variant = tx_hgvs_not_delins
										gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)	
						
							# GAP IN THE CHROMOSOME		

							elif disparity_deletion_in[0] == 'chromosome':
								# Set warning variables
								gap_position = ''
								gapped_alignment_warning = str(hgvs_genomic_5pr) + ' does not represent a true variant because it is an artefact of aligning the transcripts listed below with genome build ' + primary_assembly
								hgvs_refreshed_variant = tx_hgvs_not_delins
								gapped_transcripts = gapped_transcripts + str(hgvs_refreshed_variant.ac) + ' '	
							else:
								# Keep the same by re-setting rel_var
								hgvs_refreshed_variant = saved_hgvs_coding
				
							# Edit the output
							if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(hgvs_refreshed_variant.type)):
								hgvs_refreshed_variant = evm.n_to_c(hgvs_refreshed_variant) 
							else:
								pass
							try:
								hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
								if hgvs_refreshed_variant.posedit.edit.type == 'delins' and hgvs_refreshed_variant.posedit.edit.ref[-1] == hgvs_refreshed_variant.posedit.edit.alt[-1]:
									hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[0:-1]
									hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[0:-1]
									hgvs_refreshed_variant.posedit.pos.end.base = hgvs_refreshed_variant.posedit.pos.end.base - 1
									hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
								elif hgvs_refreshed_variant.posedit.edit.type == 'delins' and hgvs_refreshed_variant.posedit.edit.ref[0] == hgvs_refreshed_variant.posedit.edit.alt[0]:
									hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[1:]
									hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[1:]
									hgvs_refreshed_variant.posedit.pos.start.base = hgvs_refreshed_variant.posedit.pos.start.base + 1
									hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)	
							except:
								pass
							hgvs_coding = copy.deepcopy(hgvs_refreshed_variant) 
							coding = valstr(hgvs_coding)
							variant = coding
						
						# OBTAIN THE RefSeqGene coordinates
						# Attempt 1 = UTA
						sequences_for_tx = hdp.get_tx_mapping_options(hgvs_coding.ac)
						recovered_rsg = []
						for sequence in sequences_for_tx:
							if re.search('^NG_', sequence[1]):
								recovered_rsg.append(sequence[1])			
						recovered_rsg.sort()
						recovered_rsg.reverse()
						try:
							refseqgene_ac = recovered_rsg[0]
						except:
							refseqgene_ac = ''	
						# Setup gap amended RefSeqGene coding variant
						if amend_RefSeqGene == 'true':
							refseqgene_info = va_dbCrl.data.get_refseqgene_info(refseqgene_ac, primary_assembly)
							# If the RefSeqGene is an exact match to the selected genome build - map from gap corrected HGVS coding
							if refseqgene_info[0] != 'none':
								reverse_normalize_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic)
								refseqgene_coding_variant = no_norm_evm.g_to_t(reverse_normalize_hgvs_genomic, str(hgvs_coding.ac))
							else:
								refseqgene_coding_variant = hgvs_coding		
						else:
							refseqgene_coding_variant = hgvs_coding	
						ref_g_dict = va_func.refseq(str(refseqgene_coding_variant), vm, refseqgene_ac, hp, evm, hdp, primary_assembly)
						if ref_g_dict['error'] == 'false':
							hgvs_refseq = ref_g_dict['ref_g']	
							# Normalise the hgvs RefSeq variant
							try:
								hgvs_refseq = hn.normalize(hgvs_refseq)
							except hgvs.exceptions.HGVSError as e:
								error = str(e)
								hgvs_refseq = 'RefSeqGene record not available'
								refseq = 'RefSeqGene record not available'
								hgvs_refseq_ac = 'RefSeqGene record not available'
								if re.search('insertion length must be 1', error):
									hgvs_refseq = 'RefSeqGene record not available'
									refseq = 'RefSeqGene record not available'
									hgvs_refseq_ac = 'RefSeqGene record not available'
							else:
								refseq = valstr(hgvs_refseq)
								hgvs_refseq_ac = hgvs_refseq.ac
						else:
							try:
								# RETRIEVE THE REQUIRED INFORMATION FROM GENBANK database
								refseqgene_data = va_g2g.chr_to_rsg(hgvs_genomic, hn, vr)
			
								# Log the RefSeq Gene as unavailable
								hgvs_refseq = 'RefSeqGene record not available'
								refseq = 'RefSeqGene record not available'
								hgvs_refseq_ac = 'RefSeqGene record not available'
			
								if len(refseqgene_data) > 1:
									for data in refseqgene_data:
										if hgnc == data['gene'] and data['valid'] == 'true':
											refseq = data['hgvs_refseqgene']
											hgvs_refseq = hp.parse_hgvs_variant(refseq)
											refseq = valstr(hgvs_refseq)
											hgvs_refseq_ac = hgvs_refseq.ac
											break

								elif len(refseqgene_data) == 1:
									for data in refseqgene_data:
										if data['valid'] == 'true':
											refseq = data['hgvs_refseqgene']
											hgvs_refseq = hp.parse_hgvs_variant(refseq)
											refseq = valstr(hgvs_refseq)
											hgvs_refseq_ac = hgvs_refseq.ac
											break
		
								else:
									if alt_aln_method == 'GRCh37':
										evm38 = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh38', alt_aln_method=alt_aln_method)
										hgvs_38 = evm38.c_to_g(hgvs_coding)
										refseqgene_data = va_g2g.chr_to_rsg(hgvs_38, hn, vr)															
										if len(refseqgene_data) > 1:
											for data in refseqgene_data:
												if hgnc == data['gene'] and data['valid'] == 'true':
													refseq = data['hgvs_refseqgene']
													hgvs_refseq = hp.parse_hgvs_variant(refseq)
													refseq = valstr(hgvs_refseq)
													hgvs_refseq_ac = hgvs_refseq.ac
													break

										if len(refseqgene_data) == 1:
											for data in refseqgene_data:
												if data['valid'] == 'true':
													refseq = data['hgvs_refseqgene']
													hgvs_refseq = hp.parse_hgvs_variant(refseq)
													refseq = valstr(hgvs_refseq)
													hgvs_refseq_ac = hgvs_refseq.ac
													break

									if alt_aln_method == 'GRCh38':
										evm38 = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method=alt_aln_method)
										hgvs_37 = evm37.c_to_g(hgvs_coding)
										refseqgene_data = va_g2g.chr_to_rsg(hgvs_37, hn, vr)															
										if len(refseqgene_data) > 1:
											for data in refseqgene_data:
												if hgnc == data['gene'] and data['valid'] == 'true':
													refseq = data['hgvs_refseqgene']
													hgvs_refseq = hp.parse_hgvs_variant(refseq)
													refseq = valstr(hgvs_refseq)
													hgvs_refseq_ac = hgvs_refseq.ac
													break

										if len(refseqgene_data) == 1:									
											for data in refseqgene_data:
												if hgnc == data['gene'] and data['valid'] == 'true':
													refseq = data['hgvs_refseqgene']
													hgvs_refseq = hp.parse_hgvs_variant(refseq)
													refseq = valstr(hgvs_refseq)
													hgvs_refseq_ac = hgvs_refseq.ac
													break
							except KeyError as e:	
								error = str(e)
								hgvs_refseq = 'RefSeqGene record not available'
								refseq = 'RefSeqGene record not available'
								hgvs_refseq_ac = 'RefSeqGene record not available'
							
						# Predicted effect on protein
						# Translation of inversions currently supported by hgvs - So let's tackle it manually
						inversion = re.compile('inv')
						if inversion.search(variant):
							# SeqFetcher
							sf = hgvs.dataproviders.seqfetcher.SeqFetcher()
							# TO BATCH
							if re.search(':n.', variant):
								hgvs_protein = va_func.protein(variant, evm, hp)
								protein = str(hgvs_protein)
							else:
								# Convert positions to n. position
								hgvs_naughty = evm.c_to_n(hgvs_coding)
								# Collect the deleted sequence using fetch_seq
								del_seq = sf.fetch_seq(str(hgvs_naughty.ac), start_i=hgvs_naughty.posedit.pos.start.base -1, end_i=hgvs_naughty.posedit.pos.end.base)
								# Make the inverted sequence
								my_seq = Seq(del_seq)
								inv_seq = my_seq.reverse_complement()
								# Collect the associated protein
								ass_prot = hdp.get_pro_ac_for_tx_ac(hgvs_coding.ac)
								# This method sometimes fails
								if str(ass_prot) == 'None':
									cod = str(hgvs_coding)
									cod = cod.replace('inv', 'del')
									cod = hp.parse_hgvs_variant(cod)
									p = evm.c_to_p(cod)
									ass_prot = p.ac	
					
								# Intronic inversions go down as uncertain
								int_pl = re.compile('\+')
								int_mi = re.compile('\-')
								if int_pl.search(variant) or int_mi.search(variant) or re.search('\*', variant):
									# Make the variant
									hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=ass_prot, type='p', posedit='(?)')
									protein = str(hgvs_protein) 
								else:
									# Need to obtain the cds_start
									inf = va_func.tx_identity_info(variant, hdp)
									cds_start = inf[3]

									# Extract the reference coding sequence from the UTA database
									try:
										ref_seq = sf.fetch_seq(str(hgvs_naughty.ac))
									except hgvs.exceptions.HGVSError as e:
										error = str(e)
										excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)								
										validation['warnings'] = validation['warnings'] + ': ' + str(error)
										continue
									else: 
										pass
									# Create the variant coding sequence
									var_seq = variantanalyser.links.n_inversion(ref_seq, del_seq, inv_seq, hgvs_naughty.posedit.pos.start.base, hgvs_naughty.posedit.pos.end.base)
									prot_ref_seq = variantanalyser.links.translate(ref_seq, cds_start)
									prot_var_seq = variantanalyser.links.translate(var_seq, cds_start)
									if prot_ref_seq == 'error':

										error='Unable to generate protein variant description. Admin have been made aware of the issue'
										validation['warnings'] = validation['warnings'] + ': ' + str(error)
										continue
									elif prot_var_seq == 'error':
										# Does the edit affect the start codon?
										if hgvs_coding.posedit.pos.start.base >= 1 and hgvs_coding.posedit.pos.start.base <= 3:
											hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=ass_prot, type='p', posedit='(Met1?)')
											protein = str(hgvs_protein)
										else:
											error='Unable to generate protein variant description. Admin have been made aware of the issue'
											validation['warnings'] = validation['warnings'] + ': ' + str(error)
											continue

									else:				
										# Gather the required information regarding variant interval and sequences
										pro_inv_info = variantanalyser.links.pro_inv_info(prot_ref_seq, prot_var_seq)
				
										if pro_inv_info['error'] == 'true':
											error = 'Translation erroro occurred, please contact admin'
											validation['warnings'] = validation['warnings'] + ': ' + str(error)
											continue
										elif pro_inv_info['variant'] != 'true':
											# Make the variant
											hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=ass_prot, type='p', posedit='(=)')
											protein = str(hgvs_protein)
										else:
											if pro_inv_info['terminate'] == 'true':
												end = 'Ter' + str(pro_inv_info['ter_pos'])
												pro_inv_info['prot_ins_seq'].replace('*', end) 
											# Complete variant description
											iv = hgvs.location.Interval(start=pro_inv_info['edit_start'], end=pro_inv_info['edit_end'])
											# Note for hgvs to continue working, we need to take the format delXXXinsyyy
											# Need to recode the single letter del and ins sequences
											del_thr = variantanalyser.links.one_to_three(pro_inv_info['prot_del_seq'])
											ins_thr = variantanalyser.links.one_to_three(pro_inv_info['prot_ins_seq'])
											# Make the edit
											del_len = len(del_thr)
											from_aa = del_thr[0:3]
											to_aa = del_thr[del_len-3:]
											if pro_inv_info['edit_start'] != pro_inv_info['edit_end']:
												posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(pro_inv_info['edit_end']) + 'delins' + ins_thr + ')'
											else:
												posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + ins_thr + ')'
											no = 'false'
											try:
												hgvs_p = hgvs.sequencevariant.SequenceVariant(ac=ass_prot, type='p', posedit=posedit)
											except hgvs.exceptions.HGVSError as e:
												no = e
											if no == 'false':
												prot = str(hgvs_p)
												hgvs_protein = hp.parse_hgvs_variant(prot)
												protein = str(hgvs_protein)
											else:
												excep = "%s -- %s -- %s\n" %(time.ctime(), error, variant)
												validation['warnings'] = validation['warnings'] + ': ' + str(no)
												continue
		
						else: 
							try:
								hgvs_protein = va_func.protein(variant, evm, hp)
							except IndexError as e:
								error = str(e)
								if re.search('string index out of range', error) and re.search('dup', variant):
									hgvs_ins = hp.parse_hgvs_variant(variant)
									inst = hgvs_ins.ac + ':c.' + str(hgvs_ins.posedit.pos.start.base - 1) + '_' + str(hgvs_ins.posedit.pos.start.base) + 'ins' + hgvs_ins.posedit.edit.ref
									hgvs_protein = va_func.protein(inst, evm, hp)
							protein = str(hgvs_protein)

						# Gene orientation wrt genome
						ori = va_func.tx_exons(tx_ac=hgvs_coding.ac, alt_ac=hgvs_genomic.ac, alt_aln_method=alt_aln_method, hdp=hdp)
						ori = int(ori[0]['alt_strand'])
			
						# Look for normalized variant options that do not match hgvs_coding
						if ori == -1:
							# position genomic at its most 5 prime position
							try:
								query_genomic = reverse_normalize.normalize(hgvs_genomic)
							except:
								query_genomic = hgvs_genomic	
							# Map to the transcript and test for movement
							try:
								hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
							except hgvs.exceptions.HGVSError as e:
								hgvs_seek_var = saved_hgvs_coding
							else:
								seek_var = valstr(hgvs_seek_var)
								seek_ac = str(hgvs_seek_var.ac)
							if saved_hgvs_coding.posedit.edit.type != hgvs_seek_var.posedit.edit.type:
								rec_var = 'false'
								hgvs_seek_var = saved_hgvs_coding
								seek_var = valstr(hgvs_seek_var)
								seek_ac = str(hgvs_seek_var.ac)
							elif intronic_variant != 'true':
								rec_var = 'false'
								hgvs_seek_var = saved_hgvs_coding
								seek_var = valstr(hgvs_seek_var)
								seek_ac = str(hgvs_seek_var.ac)
							elif (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (saved_hgvs_coding.posedit.pos.start.base +  saved_hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (saved_hgvs_coding.posedit.pos.end.base +  saved_hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
								automap = valstr(saved_hgvs_coding) + ' normalized to ' + valstr(hgvs_seek_var) 
								hgvs_coding = hgvs_seek_var
								coding = valstr(hgvs_coding)
								validation['warnings'] = validation['warnings'] + ': ' + automap
							else:
								coding = valstr(hgvs_coding)					

						elif ori != -1:
							# position genomic at its most 3 prime position
							try:
								query_genomic = hn.normalize(hgvs_genomic)
							except:
								query_genomic = hgvs_genomic	
							# Map to the transcript and test for movement
							try:
								hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
							except hgvs.exceptions.HGVSError as e:
								hgvs_seek_var = saved_hgvs_coding
							else:
								seek_var = valstr(hgvs_seek_var)
								seek_ac = str(hgvs_seek_var.ac)
							if saved_hgvs_coding.posedit.edit.type != hgvs_seek_var.posedit.edit.type:
								rec_var = 'false'
								hgvs_seek_var = saved_hgvs_coding
								seek_var = valstr(hgvs_seek_var)
								seek_ac = str(hgvs_seek_var.ac)
							elif intronic_variant != 'true':
								rec_var = 'false'
								hgvs_seek_var = saved_hgvs_coding
								seek_var = valstr(hgvs_seek_var)
								seek_ac = str(hgvs_seek_var.ac)							
							elif (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (saved_hgvs_coding.posedit.pos.start.base +  saved_hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (saved_hgvs_coding.posedit.pos.end.base +  saved_hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
								automap = valstr(saved_hgvs_coding) + ' normalized to ' + valstr(hgvs_seek_var) 
								hgvs_coding = hgvs_seek_var
								coding = valstr(hgvs_coding)
								validation['warnings'] = validation['warnings'] + ': ' + automap
							else:
								coding = valstr(hgvs_coding)

 				# Set the data 
				validation['description'] = hgnc_gene_info
				validation['coding'] = str(hgvs_coding)
				validation['genomic_r'] = str(hgvs_refseq)
				validation['genomic_g'] = str(hgvs_genomic)
				validation['protein'] = str(hgvs_protein)
				validation['primary_assembly'] = primary_assembly

			# Report errors to User and VV admin
			except:
				if VALIDATOR_DEBUG is not None:
					import traceback
					exc_type, exc_value, last_traceback = sys.exc_info()
					te = traceback.format_exc()
					tbk = [str(exc_type), str(exc_value), str(te)]
					er = '\n'.join(tbk)
					print str(er)
					raise variantValidatorException('Validation error')
				else:	
					import warnings
					import traceback
					warnings.warn('Validation error')
					error = 'Validation error'
					validation['warnings'] = str(error)
					#exc_type, exc_value, last_traceback = sys.exc_info()
					#te = traceback.format_exc()
					#tbk = [str(exc_type), str(exc_value), str(te)]
					#er = '\n'.join(tbk)
					continue

		# Outside the for loop		
		######################

		# order the rows
		# from operator import itemgetter
		by_order = sorted(batch_list, key=itemgetter('order'))

		for valid in by_order:  
			if 'write' in valid.keys():
				if valid['write'] == 'true':
					# Blank VCF
					chr = ''
					pos = ''
					ref = ''
					alt = ''
					
					# Fromulate a json type response
					dict_out = {}
					
					# warngins 
					warnings = valid['warnings']
					warnings = re.sub('del[GATC][GATC][GATC][GATC]+', 'del', warnings)					
					warnings = re.sub('^: ', '', warnings)
					warnings = re.sub('::', ':', warnings)
					
					# Submitted variant
					submitted = valid['id']
					
					# Genomic sequence variation
					genomic_variant = valid['genomic_g']
					
					# genomic accession
					if genomic_variant != '':
						hgvs_genomic_variant = hp.parse_hgvs_variant(genomic_variant)
						genomic_variant = valstr(hgvs_genomic_variant)
						genomic_accession = hgvs_genomic_variant.ac
					else:
						genomic_accession = ''	
					
					# RefSeqGene variation
					refseqgene_variant = valid['genomic_r']
					if re.search('RefSeqGene', refseqgene_variant) or refseqgene_variant == '':
						warnings = warnings + ': ' + refseqgene_variant 
						refseqgene_variant = ''
						lrg_variant = ''	
					else:
						hgvs_refseqgene_variant = hp.parse_hgvs_variant(refseqgene_variant)
						rsg_ac = va_dbCrl.data.get_lrgID_from_RefSeqGeneID(str(hgvs_refseqgene_variant.ac))
						if rsg_ac[0] == 'none':
							lrg_variant = ''
						else:
							hgvs_lrg = copy.deepcopy(hgvs_refseqgene_variant)
							hgvs_lrg.ac = rsg_ac[0]
							lrg_variant = str(hgvs_lrg)
							if rsg_ac[1] == 'public':
								pass
							else:
								 warnings = warnings + ': The current status of ' + str(hgvs_lrg.ac) + ' is pending therefore changes may be made to the LRG reference sequence'					

					# Transcript sequence variation
					tx_variant = valid['coding']
					if tx_variant != '':	
						try:
							tx_variant = tx_variant.split('(')[1]
						except:
							pass
						tx_variant = tx_variant.replace(')', '')
					
						# transcript accession
						hgvs_transcript_variant = hp.parse_hgvs_variant(tx_variant)
						transcript_accession = hgvs_transcript_variant.ac	

						# Handle LRG	
						lrg_status = 'public'
						lrg_transcript = va_dbCrl.data.get_lrgTranscriptID_from_RefSeqTranscriptID(transcript_accession)
						if lrg_transcript == 'none':
							lrg_transcript_variant = ''
						else:
							hgvs_lrg_t = copy.deepcopy(hgvs_transcript_variant)
							hgvs_lrg_t.ac = lrg_transcript
							lrg_transcript_variant = str(hgvs_lrg_t)
					else:
						transcript_accession = ''
						lrg_transcript_variant = ''
					
					# Look for intronic variants
					if transcript_accession != '' and genomic_accession != '':
						# Remove del bases
						str_transcript = valstr(hgvs_transcript_variant)
						hgvs_transcript_variant = hp.parse_hgvs_variant(str_transcript)
						error = 'false'
						try: 
							vr.validate(hgvs_transcript_variant)
						except hgvs.exceptions.HGVSError as e:
							error = str(e)
							if re.search('intronic variant', error):
								genome_context_transcript_variant = genomic_accession + '(' + transcript_accession + '):c.' + str(hgvs_transcript_variant.posedit)
								if refseqgene_variant != '':
									hgvs_refseqgene_variant = hp.parse_hgvs_variant(refseqgene_variant)
									refseqgene_accession = hgvs_refseqgene_variant.ac
									RefSeqGene_context_transcript_variant = refseqgene_accession + '(' + transcript_accession + '):c.' + str(hgvs_transcript_variant.posedit.pos) + str(hgvs_refseqgene_variant.posedit.edit)
								else:
									RefSeqGene_context_transcript_variant = ''
							else:
								genome_context_transcript_variant = '' #transcript_variant	 
								RefSeqGene_context_transcript_variant = ''
						else:
							genome_context_transcript_variant = '' #transcript_variant	 
							RefSeqGene_context_transcript_variant = ''
					else:
						genome_context_transcript_variant = ''				
						RefSeqGene_context_transcript_variant = ''
					
					# Protein description
					pedicted_protein_variant = valid['protein']
					
					# Gene
					if transcript_accession != '':
						try:
							gene_symbol = va_dbCrl.data.get_gene_symbol_from_transcriptID(transcript_accession)
						except:
							gene_symbol = 'Unable to verify gene symbol for ' + str(transcript_accession)
					else:
						gene_symbol = ''
					
					# Transcript description
					transcript_description = valid['description']
					
					# Create VCF
					if genomic_variant != '':					
						# TO BATCH AND API AND VALIDATOR
						vcf_dict = va_H2V.report_hgvs2vcf(hgvs_genomic_variant)
						vcf_chr = vcf_dict['chr']	
						vcf_pos = vcf_dict['pos']				
						vcf_ref = vcf_dict['ref']
						vcf_alt = vcf_dict['alt']										

					# Multiple genomic variants
					multi_gen_vars = []
					if tx_variant != '':
						hgvs_coding = hp.parse_hgvs_variant(str(tx_variant))
						multi_g = []
						multi_g_tab = []
						multi_list = []
						mapping_options = hdp.get_tx_mapping_options(hgvs_coding.ac)	
						for alt_chr in mapping_options:
							if (re.match('NC_', alt_chr[1]) or re.match('NT_', alt_chr[1]) or re.match('NW_', alt_chr[1])) and alt_chr[2] == alt_aln_method:
								multi_list.append(alt_chr[1])
						for alt_chr in multi_list:
							try:
								# Re set ori
								ori = va_func.tx_exons(tx_ac=hgvs_coding.ac, alt_ac=alt_chr, alt_aln_method=alt_aln_method, hdp=hdp)
								orientation = int(ori[0]['alt_strand'])								
								hgvs_alt_genomic = va_func.myvm_t_to_g(hgvs_coding, alt_chr, vm, hn, hdp, primary_assembly)
								
								# genomic_possibilities
								# 1. take the simple 3 pr normalized hgvs_genomic
								# 2. Lock in hgvs_genomic at its most 5 prime position wrt genome
								hgvs_genomic_possibilities = []
								rn_hgvs_genomic = reverse_normalize.normalize(hgvs_alt_genomic)
								hgvs_genomic_possibilities.append(rn_hgvs_genomic)
								if orientation != -1:
									try:
										chromosome_normalized_hgvs_coding = reverse_normalize.normalize(hgvs_coding)
									except hgvs.exceptions.HGVSUnsupportedOperationError as e:
										error = str(e)
										if re.match('Normalization of intronic variants is not supported', error) or re.match('Unsupported normalization of variants spanning the exon-intron boundary', error):
											chromosome_normalized_hgvs_coding = hgvs_coding	
								else:
									try:	
										chromosome_normalized_hgvs_coding = hn.normalize(hgvs_coding)	
									except hgvs.exceptions.HGVSUnsupportedOperationError as e:
										error = str(e)
										if re.match('Normalization of intronic variants is not supported', error) or re.match('Unsupported normalization of variants spanning the exon-intron boundary', error):
											chromosome_normalized_hgvs_coding = hgvs_coding
			
								most_3pr_hgvs_genomic = va_func.myvm_t_to_g(chromosome_normalized_hgvs_coding, alt_chr, vm, hn, hdp, primary_assembly)
								hgvs_genomic_possibilities.append(most_3pr_hgvs_genomic)

								# Set variables for problem specific warnings
								gapped_alignment_warning = ''
								corrective_action_taken = ''
								gapped_transcripts = ''
								auto_info = ''
			
								# Mark as not disparity detected
								disparity_deletion_in = ['false', 'false']
								# Loop through to see if a gap can be located
								for possibility in hgvs_genomic_possibilities:
				
									# Use VCF generation code to push hgvs_genomic as for 5 prime as possible to uncover gaps
									hgvs_genomic_variant = possibility
									stored_hgvs_genomic_variant = copy.deepcopy(hgvs_genomic_variant)					
				
									# Reverse normalize hgvs_genomic_variant: NOTE will replace ref
									try:
										reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic_variant)
									except hgvs.exceptions.HGVSError as e:
										# Strange error caused by gap in genomic
										error = str(e)
										if re.search('base start position must be <= end position', error):
											if hgvs_genomic.posedit.edit.type == 'delins':
												start = hgvs_genomic.posedit.pos.start.base
												end = hgvs_genomic.posedit.pos.end.base
												lhb = sf.fetch_seq(str(hgvs_genomic.ac),end-1,end)
												rhb = sf.fetch_seq(str(hgvs_genomic.ac),start-1,start)						
												hgvs_genomic.posedit.edit.ref = lhb + rhb
												hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
												hgvs_genomic.posedit.pos.start.base = end
												hgvs_genomic.posedit.pos.end.base = start
												reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic)
											if hgvs_genomic.posedit.edit.type == 'del':
												start = hgvs_genomic.posedit.pos.start.base
												end = hgvs_genomic.posedit.pos.end.base
												lhb = sf.fetch_seq(str(hgvs_genomic.ac),end-1,end)
												rhb = sf.fetch_seq(str(hgvs_genomic.ac),start-1,start)						
												hgvs_genomic.posedit.edit.ref = lhb + rhb
												hgvs_genomic.posedit.edit.alt = lhb + rhb
												hgvs_genomic.posedit.pos.start.base = end
												hgvs_genomic.posedit.pos.end.base = start
												reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic)
										if re.search('insertion length must be 1', error):
											if hgvs_genomic.posedit.edit.type == 'ins':
												start = hgvs_genomic.posedit.pos.start.base
												end = hgvs_genomic.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(hgvs_genomic.ac),start-1,end)
												lhb = sf.fetch_seq(str(hgvs_genomic.ac),start-1,start)
												rhb = sf.fetch_seq(str(hgvs_genomic.ac),start,end)						
												hgvs_genomic.posedit.edit.ref = lhb + rhb
												hgvs_genomic.posedit.edit.alt = lhb + hgvs_genomic.posedit.edit.alt + rhb
												reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic)						

									hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)
									# Store a copy for later use
									stored_hgvs_genomic_5pr = copy.deepcopy(hgvs_genomic_5pr)
				
									# Make VCF
									vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic)
									chr = vcf_dict['chr']	
									pos = vcf_dict['pos']				
									ref = vcf_dict['ref']
									alt = vcf_dict['alt']
				
									# Look for exonic gaps within transcript or chromosome
									no_normalized_c = 'false' # Mark true to not produce an additional normalization of c.

									# Generate an end position
									end = str(int(pos) + len(ref) -1)
									pos = str(pos)
	
									# Store a not real deletion insertion to test for gapping
									stored_hgvs_not_delins = hp.parse_hgvs_variant(str(hgvs_genomic_5pr.ac) + ':' +  hgvs_genomic_5pr.type + '.' +  pos + '_' + end + 'del' + ref + 'ins' + alt)
									v = [chr, pos,ref,alt]

									# Save a copy of current hgvs_coding
									saved_hgvs_coding = no_norm_evm.g_to_t(stored_hgvs_not_delins, hgvs_coding.ac)
		
									# TO BATCH AND API
									# Detect intronic variation using normalization
									intronic_variant = 'false'
									# Look for normalized variant options that do not match hgvs_coding
									if orientation == -1:
										# position genomic at its most 5 prime position
										try:
											query_genomic = reverse_normalize.normalize(hgvs_genomic)
										except:
											query_genomic = hgvs_genomic	
										# Map to the transcript ant test for movement
										try:
											hgvs_seek_var = evm.g_to_t(query_genomic, hgvs_coding.ac)
										except hgvs.exceptions.HGVSError as e:
											hgvs_seek_var = saved_hgvs_coding
										else:		
											seek_var = valstr(hgvs_seek_var)
											seek_ac = str(hgvs_seek_var.ac)
										if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (hgvs_coding.posedit.pos.start.base +  hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (hgvs_coding.posedit.pos.end.base +  hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
											pass
										else:
											hgvs_seek_var = saved_hgvs_coding					

									elif orientation != -1:
										# position genomic at its most 3 prime position
										try:
											query_genomic = hn.normalize(hgvs_genomic)
										except:
											query_genomic = hgvs_genomic	
										# Map to the transcript and test for movement
										try:
											hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
										except hgvs.exceptions.HGVSError as e:
											hgvs_seek_var = saved_hgvs_coding
										seek_var = valstr(hgvs_seek_var)
										seek_ac = str(hgvs_seek_var.ac)
										if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (hgvs_coding.posedit.pos.start.base +  hgvs_coding.posedit.pos.start.offset) and (hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (hgvs_coding.posedit.pos.end.base +  hgvs_coding.posedit.pos.end.offset) and rec_var != 'false':
											pass
										else:
											hgvs_seek_var = saved_hgvs_coding		

									try:	
										intron_test = hn.normalize(hgvs_seek_var)	
									except hgvs.exceptions.HGVSUnsupportedOperationError as e:
										error = str(e)
										if re.match('Normalization of intronic variants is not supported', error) or re.match('Unsupported normalization of variants spanning the exon-intron boundary', error):
											# Double check to see whether the variant is actually intronic?
											for exon in ori:
												genomic_start = int(exon['alt_start_i'])
												genomic_end = int(exon['alt_end_i'])
												if (hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
													intronic_variant = 'false'
													break
												else:							
													intronic_variant = 'true'

									if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-', str(hgvs_seek_var.posedit.pos)):
										# Double check to see whether the variant is actually intronic?
										for exon in ori:
											genomic_start = int(exon['alt_start_i'])
											genomic_end = int(exon['alt_end_i'])
											if (hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
												intronic_variant = 'false'
												break							
											else:							
												intronic_variant = 'true'

									if intronic_variant != 'true':
										# Flag RefSeqGene for ammendment
										amend_RefSeqGene = 'false'
										# Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
										if stored_hgvs_not_delins != '':
											# Refresh hgvs_not_delins from stored_hgvs_not_delins
											hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
											# This test will only occur in dup of single base, insertion or substitution
											if not re.search('_', str(hgvs_not_delins.posedit.pos)):
												# For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
												plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
												plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
												plussed_hgvs_not_delins.posedit.edit.ref = ''
												transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins, str(saved_hgvs_coding.ac))
												if ((transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
													if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
														hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
														start = hgvs_not_delins.posedit.pos.start.base - 1
														end = hgvs_not_delins.posedit.pos.end.base
														ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
														hgvs_not_delins.posedit.edit.ref = ref_bases
														hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]
													elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
														hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
													elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and not re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
														hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
														start = hgvs_not_delins.posedit.pos.start.base - 1
														end = hgvs_not_delins.posedit.pos.end.base
														ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
														hgvs_not_delins.posedit.edit.ref = ref_bases
														hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] + ref_bases[1:]						
												else:
													if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
														hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
														start = hgvs_not_delins.posedit.pos.start.base - 1
														end = hgvs_not_delins.posedit.pos.end.base
														ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
														hgvs_not_delins.posedit.edit.ref = ref_bases
														hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] # + ref_bases[1:]
													elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
														hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
													elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and not re.search('del', str(hgvs_genomic_5pr.posedit.edit)):
														hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
														start = hgvs_not_delins.posedit.pos.start.base - 1
														end = hgvs_not_delins.posedit.pos.end.base
														ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac),start,end)
														hgvs_not_delins.posedit.edit.ref = ref_bases
														hgvs_not_delins.posedit.edit.alt = ref_bases[:1] + hgvs_not_delins.posedit.edit.alt[1:] # + ref_bases[1:]							
											else:
												pass						
											tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
											# Create normalized version of tx_hgvs_not_delins
											rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
											# Check for +1 base and adjust
											if re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
												# Remove offsetting to span the gap
												rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
												rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
												rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
												rn_tx_hgvs_not_delins.posedit.edit.ref = ''
												try:
													rn_tx_hgvs_not_delins.posedit.edit.alt = '' 
												except:
													pass
					
											elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
												# move tx end base to next available non-offset base
												rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
												rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
												rn_tx_hgvs_not_delins.posedit.edit.ref = ''
												if re.match('NM_', str(rn_tx_hgvs_not_delins)):
													test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
												else:
													test_tx_var = rn_tx_hgvs_not_delins	
												# re-make genomic and tx
												hgvs_not_delins = va_func.myvm_t_to_g(test_tx_var, alt_chr, vm, hn, hdp, primary_assembly)
												rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
											elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
												# move tx start base to previous available non-offset base
												rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
												rn_tx_hgvs_not_delins.posedit.edit.ref = ''
												if re.match('NM_', str(rn_tx_hgvs_not_delins)):
													test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
												else:
													test_tx_var = rn_tx_hgvs_not_delins	
												# re-make genomic and tx
												hgvs_not_delins = va_func.myvm_t_to_g(test_tx_var, alt_chr, vm, hn, hdp, primary_assembly)
												rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
												rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0							
											else:
												pass

											# Check for -ve base and adjust
											if re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
												# Remove offsetting to span the gap
												rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
												rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
												rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
												rn_tx_hgvs_not_delins.posedit.edit.ref = ''
												try:
													rn_tx_hgvs_not_delins.posedit.edit.alt = '' 
												except:
													pass
											elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
												# move tx end base back to next available non-offset base
												rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
												# Delete the ref
												rn_tx_hgvs_not_delins.posedit.edit.ref = ''
												# Add the additional base to the ALT
												start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
												end = rn_tx_hgvs_not_delins.posedit.pos.end.base
												ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac),start,end)
												rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
												if re.match('NM_', str(rn_tx_hgvs_not_delins)):
													test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
												else:
													test_tx_var = rn_tx_hgvs_not_delins	
												# re-make genomic and tx
												hgvs_not_delins = va_func.myvm_t_to_g(test_tx_var, alt_chr, vm, hn, hdp, primary_assembly)
												rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
											elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
												# move tx start base to previous available non-offset base
												rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
												rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
												rn_tx_hgvs_not_delins.posedit.edit.ref = ''
												if re.match('NM_', str(rn_tx_hgvs_not_delins)):
													test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
												else:
													test_tx_var = rn_tx_hgvs_not_delins	
												# re-make genomic and tx
												# hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, no_norm_evm, hdp, primary_assembly)
												hgvs_not_delins = va_func.myvm_t_to_g(test_tx_var, alt_chr, vm, hn, hdp, primary_assembly)
												rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, str(saved_hgvs_coding.ac))
												rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0							
											else:
												pass

											# Logic
											if len(hgvs_not_delins.posedit.edit.ref) < len(rn_tx_hgvs_not_delins.posedit.edit.ref):
												gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(hgvs_not_delins.posedit.edit.ref)
												disparity_deletion_in = ['chromosome', gap_length]
											elif len(hgvs_not_delins.posedit.edit.ref) > len(rn_tx_hgvs_not_delins.posedit.edit.ref):
												gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(rn_tx_hgvs_not_delins.posedit.edit.ref)
												disparity_deletion_in = ['transcript', gap_length]
											else:
												pass

										# Amend RefSeqGene = 'false'
										amend_RefSeqGene = 'false'
										# Recreate hgvs_genomic
										if disparity_deletion_in[0] == 'transcript':
											hgvs_genomic = hgvs_not_delins 	
						
										# GAP IN THE TRANSCRIPT	DISPARITY DETECTED	
										if disparity_deletion_in[0] == 'transcript':			
						
											amend_RefSeqGene = 'true'
											# ANY VARIANT WHOLLY WITHIN THE GAP
											if (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start))) and (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end))):
												auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
												non_valid_caution = 'true'
												if str(tx_hgvs_not_delins.posedit.pos.start) == str(tx_hgvs_not_delins.posedit.pos.end):
													if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
														tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
														tx_hgvs_not_delins.posedit.pos.end.offset = 0
													elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
														tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
														tx_hgvs_not_delins.posedit.pos.start.offset = 0									
												if re.search('\=', str(tx_hgvs_not_delins.posedit)):
													middle = tx_hgvs_not_delins.posedit.edit.alt
												elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)):
													tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
												elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
													middle = tx_hgvs_not_delins.posedit.edit.alt#[1:]
												elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
													middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
													tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
												tx_hgvs_not_delins.posedit.pos.start.offset = 0
												tx_hgvs_not_delins.posedit.pos.end.offset = 0							
												beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base)
												end = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
												tx_hgvs_not_delins.posedit.edit.ref = beginning + end
												tx_hgvs_not_delins.posedit.edit.alt = beginning + middle + end	
												hgvs_refreshed_variant = tx_hgvs_not_delins
												# Alignment position
												for_location_c = copy.deepcopy(hgvs_refreshed_variant)
												if re.match('NM_', str(for_location_c)):
													for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
												gps = for_location_c.posedit.pos.start.base -1
												gpe = for_location_c.posedit.pos.start.base
												gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'															
												# Warn update
												auto_info = auto_info + '%s' %(gap_position)
											else:
												if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)):
													auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.start.base) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
													non_valid_caution = 'true'
													tx_hgvs_not_delins.posedit.pos.start.offset = 0
													beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base) 
													middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
													end = sf.fetch_seq(str(tx_hgvs_not_delins.ac), tx_hgvs_not_delins.posedit.pos.start.base, tx_hgvs_not_delins.posedit.pos.end.base)
													tx_hgvs_not_delins.posedit.edit.ref = beginning + end
													tx_hgvs_not_delins.posedit.edit.alt = beginning + middle
													hgvs_refreshed_variant = tx_hgvs_not_delins
													# Alignment position
													for_location_c = copy.deepcopy(hgvs_refreshed_variant)
													if re.match('NM_', str(for_location_c)):
														for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
													gps = for_location_c.posedit.pos.start.base -1
													gpe = for_location_c.posedit.pos.start.base
													gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
													# Warn update
													auto_info = auto_info + '%s' %(gap_position)
												elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)):
													auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
													non_valid_caution = 'true'
													tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
													tx_hgvs_not_delins.posedit.pos.end.offset = 0
													beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.end.base - 1)
													middle = tx_hgvs_not_delins.posedit.edit.alt
													end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
													tx_hgvs_not_delins.posedit.edit.ref = beginning + end
													tx_hgvs_not_delins.posedit.edit.alt = middle + end
													hgvs_refreshed_variant = tx_hgvs_not_delins				
													# Alignment position
													for_location_c = copy.deepcopy(hgvs_refreshed_variant)
													if re.match('NM_', str(for_location_c)):
														for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
													gps = for_location_c.posedit.pos.start.base -1
													gpe = for_location_c.posedit.pos.start.base
													gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'			
													# Warn update
													auto_info = auto_info + '%s' %(gap_position)
												elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)):								
													auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.start.base) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
													non_valid_caution = 'true'
													tx_hgvs_not_delins.posedit.pos.start.offset = 0
													tx_hgvs_not_delins.posedit.pos.start.base = tx_hgvs_not_delins.posedit.pos.start.base - 1
													beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.start.base)
													middle = tx_hgvs_not_delins.posedit.edit.alt[1:]
													end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base, tx_hgvs_not_delins.posedit.pos.end.base)
													if not re.search('\=', str(tx_hgvs_not_delins.posedit)):
														tx_hgvs_not_delins.posedit.edit.alt = beginning + middle
													else:
														tx_hgvs_not_delins.posedit.edit.alt = beginning + middle + end
													tx_hgvs_not_delins.posedit.edit.ref = beginning + end
													hgvs_refreshed_variant = tx_hgvs_not_delins								
													# Alignment position
													for_location_c = copy.deepcopy(hgvs_refreshed_variant)
													if re.match('NM_', str(for_location_c)):
														for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
													gps = for_location_c.posedit.pos.start.base -1
													gpe = for_location_c.posedit.pos.start.base
													gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'
													# Warn update
													auto_info = auto_info + '%s' %(gap_position)
												elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search('\-', str(tx_hgvs_not_delins.posedit.pos.start)):
													auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac)
													non_valid_caution = 'true'
													tx_hgvs_not_delins.posedit.pos.end.offset = 0
													beginning = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.start.base - 1, tx_hgvs_not_delins.posedit.pos.end.base - 1)
													if len(tx_hgvs_not_delins.posedit.edit.ref) < len(tx_hgvs_not_delins.posedit.edit.alt):
														middle = tx_hgvs_not_delins.posedit.edit.alt.replace(tx_hgvs_not_delins.posedit.edit.ref, '' , 1)
														middle = beginning[0] + middle
													elif len(tx_hgvs_not_delins.posedit.edit.ref) > len(tx_hgvs_not_delins.posedit.edit.alt):
														middle = tx_hgvs_not_delins.posedit.edit.alt
													elif len(tx_hgvs_not_delins.posedit.edit.ref) == len(tx_hgvs_not_delins.posedit.edit.alt):
														middle = tx_hgvs_not_delins.posedit.edit.alt																
													end = sf.fetch_seq(str(tx_hgvs_not_delins.ac),tx_hgvs_not_delins.posedit.pos.end.base - 1, tx_hgvs_not_delins.posedit.pos.end.base)
													tx_hgvs_not_delins.posedit.edit.ref = beginning + end
													tx_hgvs_not_delins.posedit.edit.alt = middle + end
													hgvs_refreshed_variant = tx_hgvs_not_delins								
													# Alignment position
													for_location_c = copy.deepcopy(hgvs_refreshed_variant)
													if re.match('NM_', str(for_location_c)):
														for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
													gps = for_location_c.posedit.pos.start.base -1
													gpe = for_location_c.posedit.pos.start.base
													gap_position = ' at position ' + str(gps) + '_' + str(gpe) + '\n'							
													# Warn update
													auto_info = auto_info + '%s' %(gap_position)
												else:
													auto_info = auto_info + 'Genome position ' + str(stored_hgvs_not_delins.ac) + ':g.' + str(stored_hgvs_not_delins.posedit.pos) + ' aligns across a ' + str(disparity_deletion_in[1]) + '-bp gap in transcript ' + str(tx_hgvs_not_delins.ac) + '\n'
													hgvs_refreshed_variant = tx_hgvs_not_delins
						
										# GAP IN THE CHROMOSOME		
										elif disparity_deletion_in[0] == 'chromosome':
											amend_RefSeqGene = 'true'
											hgvs_refreshed_variant = chromosome_normalized_hgvs_coding
										else:
											# Keep the same by re-setting rel_var
											hgvs_refreshed_variant = hgvs_coding
											amend_RefSeqGene = 'false'

										# Edit the output
										if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(hgvs_refreshed_variant.type)):
											hgvs_refreshed_variant = no_norm_evm.n_to_c(hgvs_refreshed_variant)
										else:
											pass
										# Quick check to make sure the coding variant has not changed
										to_test = hn.normalize(hgvs_refreshed_variant)
										if str(to_test.posedit.edit) != str(hgvs_coding.posedit.edit):
											# Try the next available genomic option
											continue
										# Update hgvs_genomic
										hgvs_alt_genomic = va_func.myvm_t_to_g(hgvs_refreshed_variant, alt_chr, vm, hn, hdp, primary_assembly)
				
									# If it is intronic, these vairables will not have been set
									else:
										amend_RefSeqGene = 'false'
										no_normalized_c = 'false'	
				
									# Break if gap has been detected
									if disparity_deletion_in[0] != 'false':
										break		

								# TO BATCH AND API
								# Normailse hgvs_genomic
								try:
									hgvs_alt_genomic = hn.normalize(hgvs_alt_genomic)
								except hgvs.exceptions.HGVSError as e:
									# Strange error caused by gap in genomic
									error = str(e)
									if re.search('base start position must be <= end position', error) and disparity_deletion_in[0] == 'chromosome':
										if hgvs_alt_genomic.posedit.edit.type == 'delins':
											start = hgvs_alt_genomic.posedit.pos.start.base
											end = hgvs_alt_genomic.posedit.pos.end.base
											lhb = sf.fetch_seq(str(hgvs_alt_genomic.ac),end-1,end)
											rhb = sf.fetch_seq(str(hgvs_alt_genomic.ac),start-1,start)
											hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
											hgvs_alt_genomic.posedit.edit.alt = lhb + hgvs_alt_genomic.posedit.edit.alt + rhb
											hgvs_alt_genomic.posedit.pos.start.base = end
											hgvs_alt_genomic.posedit.pos.end.base = start
											hgvs_alt_genomic = hn.normalize(hgvs_alt_genomic)
										if hgvs_alt_genomic.posedit.edit.type == 'del':
											start = hgvs_alt_genomic.posedit.pos.start.base
											end = hgvs_alt_genomic.posedit.pos.end.base
											lhb = sf.fetch_seq(str(hgvs_alt_genomic.ac),end-1,end)
											rhb = sf.fetch_seq(str(hgvs_alt_genomic.ac),start-1,start)
											hgvs_alt_genomic.posedit.edit.ref = lhb + rhb
											hgvs_alt_genomic.posedit.edit.alt = lhb + rhb
											hgvs_alt_genomic.posedit.pos.start.base = end
											hgvs_alt_genomic.posedit.pos.end.base = start
											hgvs_alt_genomic = hn.normalize(hgvs_alt_genomic)

								# Refresh the :g. variant
								multi_g.append(str(hgvs_alt_genomic))
							except:
								import os
								import traceback
								exc_type, exc_value, last_traceback = sys.exc_info()
								te = traceback.format_exc()
								error = str(te) 
								continue
						if multi_g != []:
							multi_g.sort()				
							multi_gen_vars = multi_g #'|'.join(multi_g)
						else:
							multi_gen_vars = []		
					
					else:
						multi_gen_vars = []													
					
					# Warn not directly mapped to specified genome build		
					if genomic_variant != '':
						caution = ''
						if not re.match('NC_', str(hgvs_genomic.ac)):
							warnings = warnings + ': ' + str(hgvs_coding) + ' can not be mapped directly to genome build ' + primary_assembly + '. See Alternative genomic loci for aligned genomic positions' 
							caution = 'can not be mapped directly to genome build'
						chr_num = va_scb.supported_for_mapping(str(hgvs_genomic.ac), primary_assembly)
						if chr_num ==  'false':
							already_seen = 'can not be mapped directly to genome build'
							if re.search(already_seen, caution):
								pass
							else:	
								warnings = warnings + ': ' + str(hgvs_coding) + ' can not be mapped directly to genome build ' + primary_assembly + '. See alt_genomic_loci for aligned genomic positions'

					warn_list = warnings.split(': ')
					warnings_out = []
					for warning in warn_list:
						warning.strip()
						if warning == '':
							continue
						warnings_out.append(warning)
					
					# Populate the dictionary
					dict_out['submitted_variant'] = submitted
					dict_out['HGVS_genomic_variant'] = 	genomic_variant 
					dict_out['gene_symbol'] = gene_symbol
					dict_out['transcript_description'] = transcript_description
					dict_out['HGVS_transcript_variant'] = tx_variant
					dict_out['genome_context_intronic_sequence'] = genome_context_transcript_variant
					dict_out['RefSeqGene_context_intronic_sequence'] = RefSeqGene_context_transcript_variant
					dict_out['HGVS_RefSeqGene_variant']	= refseqgene_variant
					dict_out['HGVS_predicted_protein_consequence'] = pedicted_protein_variant
					dict_out['validation_warnings'] = warnings_out
					dict_out['HGVS_LRG_transcript_variant'] = lrg_transcript_variant
					dict_out['HGVS_LRG_variant'] = lrg_variant
					dict_out['alt_genomic_loci'] = multi_gen_vars
					# vcf
					try:
						dict_out['vcf'] = {'chr' : vcf_chr, 'pos' : vcf_pos, 'ref': vcf_ref, 'alt' : vcf_alt}
					except:
						dict_out['vcf'] = {'chr' : '', 'pos' : '', 'ref': '', 'alt' : ''}	 	
					
					# Append to a list for return	
					batch_out.append(dict_out)
				else:
					continue
			else:
				continue
	
		# Exit the script
		return batch_out

	# Bug catcher
	except:
		# Debug mode
		if VALIDATOR_DEBUG is not None:
			import traceback
			exc_type, exc_value, last_traceback = sys.exc_info()
			te = traceback.format_exc()
			# tr = ''.join(traceback.format_stack())
			tbk = [str(exc_type), str(exc_value), str(te)]
			er = '\n'.join(tbk)
			print str(er)
		# Report and raise error
		error = [{'validation_warnings': 'Validation error'}]
		raise variantValidatorException('Validation error')
 		# Return
 		return
 		
# Generates a list of transcript (UTA supported) and transcript names from a gene symbol or RefSeq transcript ID
def gene2transcripts(query):
	caution = ''
	input = query
	input = input.upper()
	# Quick check for blank form
	if input == '':
		caution = {'error' : 'Please enter HGNC gene name or transcript identifier (NM_, NR_, or ENST)'}
		return caution
	else:
		caution = ''
		hgnc = input
		if re.match('NM_', hgnc) or re.match('NR_', hgnc): #  or re.match('ENST', hgnc):
			try:
				tx_info = hdp.get_tx_identity_info(hgnc)
				hgnc = tx_info[6] 
			except hgvs.exceptions.HGVSError as e:
				caution = {'error' : str(e)}
				return caution		
		# Look up current name
		current = va_func.hgnc_rest(path = "/search/prev_symbol/" + hgnc)
		# Look for historic names
		# If historic names = 0
		if str(current['record']['response']['numFound']) == '0':
			current_sym = hgnc
		else:
			current_sym = current['record']['response']['docs'][0]['symbol']
		# Look up previous symbols and gene name
		previous = va_func.hgnc_rest(path = "/fetch/symbol/" + current_sym)
		try:
			previous_sym = previous['record']['response']['docs'][0]['prev_symbol'][0]
		except:
			previous_sym = current_sym
	
		# Get gene name
		try:
			gene_name = previous['record']['response']['docs'][0]['name']
		except:
			caution = current_sym + ' is not a valid HGNC gene symbol'
			gene_name = 'Not found in the HGNC database of human gene names www.genenames.org'  
		
		# Look up previous name
		try:
			previous_name = previous['record']['response']['docs'][0]['prev_name'][0]
		except:
			previous_name = gene_name		

		# Get transcripts
		tx_for_gene = hdp.get_tx_for_gene(current_sym)
		if len(tx_for_gene) == 0:
			tx_for_gene = hdp.get_tx_for_gene(previous_sym)
		if len(tx_for_gene) == 0:
			tx_for_gene = {'error' : 'Unable to retrieve data from the UTA, please contact admin'}
		
		# Loop through each transcript and get the relevant transcript description
		genes_and_tx = []
		recovered_dict = {}
		for line in tx_for_gene:
			if  re.match('^NM_', line[3]) or re.match('^NR_', line[3]):
				# Transcript ID
				tx = line[3]
				tx_description = va_dbCrl.data.get_transcript_description(tx)
				# Check for duplicates
				if tx in recovered_dict.keys():
					continue
				else:	
					try:
						# Add to recovered_dict
						recovered_dict[tx] = ''
						genes_and_tx.append([tx, tx_description, line[1] + 1, line[2]])
					except:
						# Add to recovered_dict
						recovered_dict[tx] = ''
						genes_and_tx.append([tx, tx_description, 'not applicable', 'not applicable'])	 		
					# LRG information
					lrg_transcript = va_dbCrl.data.get_lrgTranscriptID_from_RefSeqTranscriptID(tx)
					if lrg_transcript == 'none':
						pass
					else:		
						genes_and_tx.append([lrg_transcript, tx_description, line[1] + 1, line[2]]) 		

 		cp_genes_and_tx = copy.deepcopy(genes_and_tx)
 		genes_and_tx = []
 		for tx in cp_genes_and_tx:
 			tx_d = {'reference' : tx[0],
 					'description' : tx[1],
 					'coding_start' : tx[2],
 					'coding_end' : tx[3]}
 			genes_and_tx.append(tx_d)		
 			
 		
 		# Return data table
 		g2d_data = {'current_symbol' : current_sym,
 					'previous_symbol' : previous_sym,
 					'current_name' : gene_name,
 					'previous_name' : previous_name,
 					'transcripts' :  genes_and_tx,
 					'error' : ''}  
 		
 		return g2d_data
 		
 		
 		
 		
 		
 		
 		
 		
 		