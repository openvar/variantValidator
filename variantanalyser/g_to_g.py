import os
import re
import functions
import hgvs
import hgvs.parser
import hgvs.exceptions
import dbControls
import dbControls.data as database_data

hp = hgvs.parser.Parser()
alt_aln_method = 'splign'

# From the hgvs parser import, create an instance of hgvs.parser.Parser
hp = hgvs.parser.Parser() 			

# Set file root
# Set up os paths data and log folders 
# FILE_ROOT = os.path.dirname(os.path.abspath(__file__))

# Covert chromosomal HGVS description to RefSeqGene
def chr_to_rsg(hgvs_genomic, hn, vr):
	# print 'chr_to_rsg triggered'
	hgvs_genomic = hn.normalize(hgvs_genomic)
	# split the description
	# Accessions
	chr_ac = hgvs_genomic.ac
	# Positions
	chr_start_pos = int(hgvs_genomic.posedit.pos.start.base)
	chr_end_pos = int(hgvs_genomic.posedit.pos.end.base)
	# edit
	chr_edit = hgvs_genomic.posedit.edit

	# Pre set variable, note there could be several
	rsg_data_set = []
	
	# Recover table from MySql
	all_info = database_data.get_g_to_g_info()
	for line in all_info:
		# # print line
		# Logic to identify the correct RefSeqGene
		# # print str(data)
		rsg_data = {}
		if chr_ac == line[1] and chr_start_pos >= int(line[2]) and chr_end_pos <= int(line[3]):
		#query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol FROM refSeqGene_loci"
			# (u'NG_034189.1', u'NC_000004.12', 190173122, 190177845, u'+', u'DUX4L1')
			# Set the values of the data dictionary
			rsg_data['rsg_ac'] = line[0]
			rsg_data['chr_ac'] = line[1]
			rsg_data['rsg_start'] = line[2]
			rsg_data['rsg_end'] = line[3]
			rsg_data['ori'] = line[4]
			rsg_data['gene'] = line[5]
			rsg_data_set.append(rsg_data)		
		else:
			continue
	
	# Compile descriptions and validate
	descriptions = []
	for rsg_data in rsg_data_set:
		rsg_ac = rsg_data['rsg_ac']
		rsg_start =	rsg_data['rsg_start']
		rsg_end = rsg_data['rsg_end']
		ori = rsg_data['ori']
		gene = rsg_data['gene']
		# String the description
		if ori == '+':
			rsg_description = rsg_ac + ':g.' + str(chr_start_pos - int(rsg_start)+1) + '_' + str(chr_end_pos - int(rsg_start)+1) + str(chr_edit)
			# # print '\n\n' + rsg_description + '\n\n'
			hgvs_refseqgene = hp.parse_hgvs_variant(rsg_description)
			try:
				hgvs_refseqgene = hn.normalize(hgvs_refseqgene)
			except:
				error = 'Not in SeqRepo'
				data = {'hgvs_refseqgene' : str(hgvs_refseqgene), 'gene' : gene, 'valid' : str(error)}
				descriptions.append(data)
				continue 	
			try:
				vr.validate(hgvs_refseqgene)
			except hgvs.exceptions.HGVSError as e:
				error = str(e)
				if re.search('does not agree with reference sequence', error):
					match = re.findall('\(([GATC]+)\)', error)
					new_ref = match[1]
					#if re.search('=', str(hgvs_refseqgene.posedit.edit)):
					#	hgvs_refseqgene.posedit.edit.ref = new_ref
					#	hgvs_refseqgene.posedit.edit.alt = new_ref
					#else: 
					hgvs_refseqgene.posedit.edit.ref = new_ref
					error = 'true'	
				else:
					pass	
				data = {'hgvs_refseqgene' : str(hgvs_refseqgene), 'gene' : gene, 'valid' : str(error)}
			else:
				data = {'hgvs_refseqgene' : str(hgvs_refseqgene), 'gene' : gene, 'valid' : 'true'}
			descriptions.append(data)
		if ori == '-':
			# Reverse complement of bases may be required. Let normalizer do the lifting for strings of bases
			# Look for scenarios with RC needed bases and extract the bases from the edit
			if re.search(r"((del[GATCUgatcu]+))", str(chr_edit)):
				bases = re.search(r"((del[GATCUgatcu]+))", str(chr_edit))
				bases = bases.group(1)
				chr_edit = 'del' + str(chr_edit).replace(bases, '')
			if re.search(r"((ins[GATCUgatcu]+))", str(chr_edit)):
				bases = re.search(r"((ins[GATCUgatcu]+))", str(chr_edit))
				bases = bases.group(1)
				ins_revcomp = functions.revcomp(bases)
				chr_edit = str(chr_edit).replace(bases, '') + 'ins' + ins_revcomp
			if re.search(r"((dup[GATCUgatcu]+))", str(chr_edit)):
				bases = re.search(r"((dup[GATCUgatcu]+))", str(chr_edit))
				bases = bases.group(1)
				chr_edit = 'dup' + str(chr_edit).replace(bases, '')
			if re.search(r"((inv[GATCUgatcu]+))", str(chr_edit)):
				bases = re.search(r"((inv[GATCUgatcu]+))", str(chr_edit))
				bases = bases.group(1)
				chr_edit = 'inv' + str(chr_edit).replace(bases, '')
			if re.search('>', str(chr_edit)):
				chr_edit = str(chr_edit)
				chr_edit = chr_edit.replace('A>', 't>')
				chr_edit = chr_edit.replace('T>', 'a>')
				chr_edit = chr_edit.replace('G>', 'c>')
				chr_edit = chr_edit.replace('C>', 'g>')
				chr_edit = chr_edit.replace('>A', '>t')
				chr_edit = chr_edit.replace('>T', '>a')
				chr_edit = chr_edit.replace('>G', '>c')
				chr_edit = chr_edit.replace('>C', '>g')
				chr_edit = chr_edit.replace('C=', 'g=')
				chr_edit = chr_edit.replace('G=', 'c=')
				chr_edit = chr_edit.replace('A=', 't=')
				chr_edit = chr_edit.replace('T=', 'a=')
				chr_edit = chr_edit.upper()
				
			rsg_description = rsg_ac + ':g.' + str( (int(rsg_end) - int(rsg_start)) -  (chr_end_pos - int(rsg_start))+1) + '_' + str( (int(rsg_end) - int(rsg_start)) -  (chr_start_pos - int(rsg_start))+1) + str(chr_edit)
			hgvs_refseqgene = hp.parse_hgvs_variant(rsg_description)
			try:
				hgvs_refseqgene = hn.normalize(hgvs_refseqgene)
			except:
				error = 'Not in SeqRepo'
				data = {'hgvs_refseqgene' : str(hgvs_refseqgene), 'gene' : gene, 'valid' : str(error)}
				descriptions.append(data)
				continue 
			try:
				vr.validate(hgvs_refseqgene)
			except hgvs.exceptions.HGVSError as e:
				error = str(e)
				if re.search('does not agree with reference sequence', error):
					match = re.findall('\(([GATC]+)\)', error)
					new_ref = match[1]
					#if re.search('=', str(hgvs_refseqgene.posedit.edit)):
					#	hgvs_refseqgene.posedit.edit.ref = new_ref
					#	hgvs_refseqgene.posedit.edit.alt = new_ref
					#else: 
					hgvs_refseqgene.posedit.edit.ref = new_ref
					error = 'true'	
				else:
					pass		
				data = {'hgvs_refseqgene' : str(hgvs_refseqgene), 'gene' : gene, 'valid' : str(error)}
			else:
				data = {'hgvs_refseqgene' : str(hgvs_refseqgene), 'gene' : gene, 'valid' : 'true'}
			descriptions.append(data)

	# Return the required data. This is a dictionary containing the rsg description, validation status and gene ID
	return descriptions

# Covert RefSeqGene HGVS description to Chromosomal
def rsg_to_chr(hgvs_refseqgene, primary_assembly, hn, vr):
	# normalize
	try:
		hgvs_refseqgene = hn.normalize(hgvs_refseqgene)
	except:
		pass
	# split the description
	# Accessions
	rsg_ac = hgvs_refseqgene.ac
	# Positions
	rsg_start_pos = int(hgvs_refseqgene.posedit.pos.start.base)
	rsg_end_pos = int(hgvs_refseqgene.posedit.pos.end.base)
	# edit
	rsg_edit = hgvs_refseqgene.posedit.edit
	
	# Pre set variable, note there could be several
	chr_data_set = []
	# Recover table from MySql
	all_info = database_data.get_g_to_g_info()
	for line in all_info:
		# Logic to identify the correct RefSeqGene
		chr_data = {}
		if rsg_ac == line[0] and primary_assembly == line[6]:
		#query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol FROM refSeqGene_loci"
			# (u'NG_034189.1', u'NC_000004.12', 190173122, 190177845, u'+', u'DUX4L1')
			# Set the values of the data dictionary
			chr_data['rsg_ac'] = line[0]
			chr_data['chr_ac'] = line[1]
			chr_data['rsg_start'] = line[2]
			chr_data['rsg_end'] = line[3]
			chr_data['ori'] = line[4]
			chr_data['gene'] = line[5]
			chr_data_set.append(chr_data)
		else:
			continue
	
	# Compile descriptions and validate
	descriptions = []
	for chr_data in chr_data_set:
		chr_ac = chr_data['chr_ac']
		rsg_ac = chr_data['rsg_ac']
		chr_start =	int(chr_data['rsg_start'])
		chr_end = int(chr_data['rsg_end'])
		ori = chr_data['ori']
		gene = chr_data['gene']
		# String the description
		if ori == '+':
			chr_description = chr_ac + ':g.' + str(chr_start + rsg_start_pos -1)  + '_' + str(chr_start + rsg_end_pos -1) + str(rsg_edit)
			hgvs_genomic = hp.parse_hgvs_variant(chr_description)
			hgvs_genomic = hn.normalize(hgvs_genomic)
			try:
				vr.validate(hgvs_genomic)
			except hgvs.exceptions.HGVSError as e:
				error = str(e)
				if re.search('does not agree with reference sequence', error):
					match = re.findall('\(([GATC]+)\)', error)
					new_ref = match[1]
					#if re.search('=', str(hgvs_genomic.posedit.edit)):
					#	hgvs_genomic.posedit.edit.ref = new_ref
					#	hgvs_genomic.posedit.edit.alt = new_ref
					#else: 
					hgvs_genomic.posedit.edit.ref = new_ref
					error = 'true'	
				else:
					pass	
				# # print str(e) + '\n3.'	
				data = {'hgvs_genomic' : str(hgvs_genomic), 'gene' : gene, 'valid' : str(error)}
			else:
				data = {'hgvs_genomic' : str(hgvs_genomic), 'gene' : gene, 'valid' : 'true'}
			descriptions.append(data)
		if ori == '-':
			# Reverse complement of bases may be required. Let normalizer do the lifting for strings of bases
			# Look for scenarios with RC needed bases and extract the bases from the edit
			if re.search(r"((del[GATCUgatcu]+))", str(rsg_edit)):
				bases = re.search(r"((del[GATCUgatcu]+))", str(rsg_edit))
				bases = bases.group(1)
				rsg_edit = 'del' + str(rsg_edit).replace(bases, '')
			if re.search(r"((ins[GATCUgatcu]+))", str(rsg_edit)):
				bases = re.search(r"((ins[GATCUgatcu]+))", str(rsg_edit))
				bases = bases.group(1)
				ins_revcomp = functions.revcomp(bases)
				rsg_edit = str(rsg_edit).replace(bases, '') + 'ins' + ins_revcomp
			if re.search(r"((dup[GATCUgatcu]+))", str(rsg_edit)):
				bases = re.search(r"((dup[GATCUgatcu]+))", str(rsg_edit))
				bases = bases.group(1)
				rsg_edit = 'dup' + str(rsg_edit).replace(bases, '')
			if re.search(r"((inv[GATCUgatcu]+))", str(rsg_edit)):
				bases = re.search(r"((inv[GATCUgatcu]+))", str(rsg_edit))
				bases = bases.group(1)
				rsg_edit = 'inv' + str(rsg_edit).replace(bases, '')
			if re.search('>', str(rsg_edit)):
				rsg_edit = str(rsg_edit)
				rsg_edit = rsg_edit.replace('A>', 't>')
				rsg_edit = rsg_edit.replace('T>', 'a>')
				rsg_edit = rsg_edit.replace('G>', 'c>')
				rsg_edit = rsg_edit.replace('C>', 'g>')
				rsg_edit = rsg_edit.replace('>A', '>t')
				rsg_edit = rsg_edit.replace('>T', '>a')
				rsg_edit = rsg_edit.replace('>G', '>c')
				rsg_edit = rsg_edit.replace('>C', '>g')
				rsg_edit = rsg_edit.replace('C=', 'g=')
				rsg_edit = rsg_edit.replace('G=', 'c=')
				rsg_edit = rsg_edit.replace('A=', 't=')
				rsg_edit = rsg_edit.replace('T=', 'a=')
				rsg_edit = rsg_edit.upper()
				
			chr_description =  chr_ac + ':g.' + str( int(chr_start) + (int(chr_end) - int(chr_start)) - rsg_end_pos + 1 ) + '_' + str( int(chr_start) + (int(chr_end) - int(chr_start)) - rsg_start_pos + 1 ) + str(rsg_edit)
		 
			hgvs_genomic = hp.parse_hgvs_variant(chr_description)
			hgvs_genomic = hn.normalize(hgvs_genomic)
			try:
				vr.validate(hgvs_genomic)
			except hgvs.exceptions.HGVSError as e:
				error = str(e)
				if re.search('does not agree with reference sequence', error):
					match = re.findall('\(([GATC]+)\)', error)
					new_ref = match[1]
					#if re.search('=', str(hgvs_genomic.posedit.edit)):
					#	hgvs_genomic.posedit.edit.ref = new_ref
					#	hgvs_genomic.posedit.edit.alt = new_ref
					#else: 
					hgvs_genomic.posedit.edit.ref = new_ref
					error = 'true'
				data = {'hgvs_genomic' : str(hgvs_genomic), 'gene' : gene, 'valid' : str(error)}
			else:
				data = {'hgvs_genomic' : str(hgvs_genomic), 'gene' : gene, 'valid' : 'true'}
			descriptions.append(data)

	# Return the required data. This is a dictionary containing the rsg description, validation status and gene ID
	return descriptions			
					