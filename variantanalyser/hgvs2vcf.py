# Import  modules
import re
import copy
import hgvs
import hgvs.dataproviders
import hgvs.normalizer
import supported_chromosome_builds

# Import Biopython modules
from Bio.Seq import Seq

# Set variables
hdp = hgvs.dataproviders.uta.connect(pooling=True)

# Reverse normalizer (5 prime)
reverse_normalize = hgvs.normalizer.Normalizer(hdp, 
		cross_boundaries=False, 
		shuffle_direction=5, 
		alt_aln_method='splign'
		)

# normalizer (5 prime)
normalize = hgvs.normalizer.Normalizer(hdp, 
		cross_boundaries=False, 
		shuffle_direction=3, 
		alt_aln_method='splign'
		)

# SeqFetcher
sf = hgvs.dataproviders.seqfetcher.SeqFetcher()

def hgvs2vcf(hgvs_genomic):		

	hgvs_genomic_variant = hgvs_genomic					
	# Reverse normalize hgvs_genomic_variant: NOTE will replace ref
	reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic_variant)
	hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

	# Chr
	chr = supported_chromosome_builds.to_chr_num(reverse_normalized_hgvs_genomic.ac)
	if chr is not None:
		pass
	else:
		chr = reverse_normalized_hgvs_genomic.ac		
	
	
	if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
		pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
		ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
		alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

	# Insertions	
	elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(reverse_normalized_hgvs_genomic.posedit))):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base) 
		alt_start = start - 1 #
		# Recover sequences
		ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),alt_start,end-1)
		ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		# Assemble  
		pos = start
		ref = ref_seq
		alt = ref_seq + ins_seq						

	# Substitutions
	elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
		ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
		alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
		pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

	# Deletions
	elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2
		start = start - 1
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,start)
		# Assemble  
		pos = str(start)
		ref = pre_base + hgvs_del_seq
		alt = pre_base

	
	# inv
	elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start -1
		start = start
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = ins_seq
		if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
			my_seq = Seq(vcf_del_seq)
			alt = str(my_seq.reverse_complement()) 
	
	# Delins
	elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(reverse_normalized_hgvs_genomic.posedit))):
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base -1)
		adj_start = start -1
		start = start
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = vcf_del_seq[:1] + ins_seq			
	
	
	# Duplications								
	elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base) #
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2 #
		start = start - 1 #
		# Recover sequences
		dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
		vcf_ref_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_ref_seq
		alt = vcf_ref_seq + dup_seq
	else:
		chr = ''
		ref = ''
		alt = ''
		pos = ''
		
	
	# ensure as 3' as possible
	if chr != '' and pos != '' and ref != '' and alt != '':
		if len(ref) > 1:
			rsb = list(str(ref))
			# if rsb[0] == rsb[1] and reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
			if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
				pos = int(pos) - 1
				prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),pos-1,pos)
				pos = str(pos)				
				ref = prev + ref
				alt = prev + alt
	
	# Dictionary the VCF
	vcf_dict = {'chr' : chr, 'pos' : pos, 'ref' : ref, 'alt' : alt}
	return vcf_dict		
	
def report_hgvs2vcf(hgvs_genomic):		

	hgvs_genomic_variant = hgvs_genomic					
	# Reverse normalize hgvs_genomic_variant: NOTE will replace ref
	reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic_variant)
	hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

	# Chr
	chr = supported_chromosome_builds.to_chr_num(reverse_normalized_hgvs_genomic.ac)
	if chr is not None:
		pass
	else:
		chr = reverse_normalized_hgvs_genomic.ac

	# UCSC Chr
	ucsc_chr = supported_chromosome_builds.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac)
	if ucsc_chr is not None:
		pass
	else:
		ucsc_chr = reverse_normalized_hgvs_genomic.ac

	# GRC Chr
	grc_chr = supported_chromosome_builds.to_chr_num_grc(reverse_normalized_hgvs_genomic.ac)
	if grc_chr is not None:
		pass
	else:
		grc_chr = reverse_normalized_hgvs_genomic.ac			
				


	# TO BATCH AND API AND VALIDATOR
	if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
		pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
		ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
		alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

	# Insertions	
	elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(reverse_normalized_hgvs_genomic.posedit))):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base) 
		alt_start = start - 1 #
		# Recover sequences
		ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),alt_start,end-1)
		ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		# Assemble  
		pos = start
		ref = ref_seq
		alt = ref_seq + ins_seq						

	# Substitutions
	elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
		ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
		alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
		pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

	# Deletions
	elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2
		start = start - 1
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,start)
		# Assemble  
		pos = str(start)
		ref = pre_base + hgvs_del_seq
		alt = pre_base

	
	# inv
	elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start -1
		start = start
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = ins_seq
		if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
			my_seq = Seq(vcf_del_seq)
			alt = str(my_seq.reverse_complement()) 
	
	# Delins
	elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(reverse_normalized_hgvs_genomic.posedit))):
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base -1)
		adj_start = start -1
		start = start
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = vcf_del_seq[:1] + ins_seq			
	
	
	# Duplications								
	elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base) #
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2 #
		start = start - 1 #
		# Recover sequences
		dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
		vcf_ref_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_ref_seq
		alt = vcf_ref_seq + dup_seq
	else:
		chr = ''
		ucsc_chr = ''
		grc_chr = ''
		ref = ''
		alt = ''
		pos = ''
	
	# Dictionary the VCF
	vcf_dict = {'pos' : pos, 'ref' : ref, 'alt' : alt, 'ucsc_chr' : ucsc_chr, 'grc_chr' : grc_chr}
	return vcf_dict		
	
def pos_lock_hgvs2vcf(hgvs_genomic):		
					
	# DO NOT Reverse normalize
	reverse_normalized_hgvs_genomic = hgvs_genomic
	hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

	# Chr
	if re.match('NC_', reverse_normalized_hgvs_genomic.ac):
		chr = supported_chromosome_builds.to_chr_num(reverse_normalized_hgvs_genomic.ac)
		if chr is not None:
			pass
		else:
			chr = reverse_normalized_hgvs_genomic.ac		
	else:
		chr = reverse_normalized_hgvs_genomic.ac	

	# TO BATCH AND API AND VALIDATOR
	if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
		pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
		ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
		alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

	# Insertions	
	elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(reverse_normalized_hgvs_genomic.posedit))):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base) 
		alt_start = start - 1 #
		# Recover sequences
		ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),alt_start,end-1)
		ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		# Assemble  
		pos = start
		ref = ref_seq
		alt = ref_seq + ins_seq						

	# Substitutions
	elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
		ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
		alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
		pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

	# Deletions
	elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2
		start = start - 1
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,start)
		# Assemble  
		pos = str(start)
		ref = pre_base + hgvs_del_seq
		alt = pre_base

	
	# inv
	elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start -1
		start = start
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = ins_seq
		if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
			my_seq = Seq(vcf_del_seq)
			alt = str(my_seq.reverse_complement()) 
	
	# Delins
	elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(reverse_normalized_hgvs_genomic.posedit))):
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base -1)
		adj_start = start -1
		start = start
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = vcf_del_seq[:1] + ins_seq			
	
	
	# Duplications								
	elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base) #
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2 #
		start = start - 1 #
		# Recover sequences
		dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
		vcf_ref_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_ref_seq
		alt = vcf_ref_seq + dup_seq
	else:
		chr = ''
		ref = ''
		alt = ''
		pos = ''
		
	vcf_dict = {'chr' : chr, 'pos' : pos, 'ref' : ref, 'alt' : alt}
	return vcf_dict	



def hard_right_hgvs2vcf(hgvs_genomic):		
	hgvs_genomic_variant = hgvs_genomic					
	# Reverse normalize hgvs_genomic_variant: NOTE will replace ref
	normalized_hgvs_genomic = normalize.normalize(hgvs_genomic_variant)
	hgvs_genomic_5pr = copy.deepcopy(normalized_hgvs_genomic)

	# Chr
	if re.match('NC_', normalized_hgvs_genomic.ac):
		chr = supported_chromosome_builds.to_chr_num(normalized_hgvs_genomic.ac)
		if chr is not None:
			pass
		else:
			chr = normalized_hgvs_genomic.ac		
	else:
		chr = normalized_hgvs_genomic.ac	

	# TO BATCH AND API AND VALIDATOR
	if re.search('[GATC]+\=', str(normalized_hgvs_genomic.posedit)):
		pos = str(normalized_hgvs_genomic.posedit.pos.start)
		ref = normalized_hgvs_genomic.posedit.edit.ref
		alt = normalized_hgvs_genomic.posedit.edit.ref

	# Insertions	
	elif (re.search('ins', str(normalized_hgvs_genomic.posedit)) and not re.search('del', str(normalized_hgvs_genomic.posedit))):						
		end = int(normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(normalized_hgvs_genomic.posedit.pos.start.base) 
		alt_start = start - 1 #
		# Recover sequences
		ref_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac),alt_start,end-1)
		ins_seq = normalized_hgvs_genomic.posedit.edit.alt
		# Assemble  
		pos = start
		ref = ref_seq
		alt = ref_seq + ins_seq						

	# Substitutions
	elif re.search('>', str(normalized_hgvs_genomic.posedit)):
		ref = normalized_hgvs_genomic.posedit.edit.ref
		alt = normalized_hgvs_genomic.posedit.edit.alt
		pos = str(normalized_hgvs_genomic.posedit.pos)

	# Deletions
	elif re.search('del', str(normalized_hgvs_genomic.posedit)) and not re.search('ins', str(normalized_hgvs_genomic.posedit)):						
		end = int(normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2
		start = start - 1
		try:
			ins_seq = normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac),start,end)
		pre_base = sf.fetch_seq(str(normalized_hgvs_genomic.ac),adj_start,start)
		# Assemble  
		pos = str(start)
		ref = pre_base + hgvs_del_seq
		alt = pre_base
	
	# inv
	elif re.search('inv', str(normalized_hgvs_genomic.posedit)):						
		end = int(normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start -1
		start = start
		try:
			ins_seq = normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = ins_seq
		if re.search('inv', str(normalized_hgvs_genomic.posedit)):
			my_seq = Seq(vcf_del_seq)
			alt = str(my_seq.reverse_complement()) 
	
	# Delins
	elif (re.search('del', str(normalized_hgvs_genomic.posedit)) and re.search('ins', str(normalized_hgvs_genomic.posedit))):
		end = int(normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(normalized_hgvs_genomic.posedit.pos.start.base -1)
		adj_start = start -1
		start = start
		try:
			ins_seq = normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = vcf_del_seq[:1] + ins_seq			
	
	
	# Duplications								
	elif (re.search('dup', str(normalized_hgvs_genomic.posedit))):
		end = int(normalized_hgvs_genomic.posedit.pos.end.base) #
		start = int(normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2 #
		start = start - 1 #
		# Recover sequences
		dup_seq = normalized_hgvs_genomic.posedit.edit.ref
		vcf_ref_seq	= sf.fetch_seq(str(normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_ref_seq
		alt = vcf_ref_seq + dup_seq
	else:
		chr = ''
		ref = ''
		alt = ''
		pos = ''
		
	
	# ADD SURROUNDING BASES
	if chr != '' and pos != '' and ref != '' and alt != '':
        # Add 2 post bases
		pos = int(pos)
		pre_end_pos = pos + len(ref)
		end_pos = pre_end_pos + 1
		post = sf.fetch_seq(str(normalized_hgvs_genomic.ac),pre_end_pos-1,end_pos)
		ref = ref + post
		alt = alt + post
# 		Add another pre-base
# 		pre_pos = int(pos) - 1
# 		pre_pos 
# 		prev = sf.fetch_seq(str(normalized_hgvs_genomic.ac),pre_pos-1,pre_pos)
# 		pos = str(pre_pos)				
# 		ref = prev + ref
# 		alt = prev + alt
	
	# Dictionary the VCF
	vcf_dict = {'chr' : chr, 'pos' : pos, 'ref' : ref, 'alt' : alt}
	return vcf_dict	

def hard_left_hgvs2vcf(hgvs_genomic):		
	hgvs_genomic_variant = hgvs_genomic					
	# Reverse normalize hgvs_genomic_variant: NOTE will replace ref
	reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic_variant)
	hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

	# Chr
	if re.match('NC_', reverse_normalized_hgvs_genomic.ac):
		chr = supported_chromosome_builds.to_chr_num(reverse_normalized_hgvs_genomic.ac)
		if chr is not None:
			pass
		else:
			chr = reverse_normalized_hgvs_genomic.ac		
	else:
		chr = reverse_normalized_hgvs_genomic.ac	

	if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
		pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
		ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
		alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

	# Insertions	
	elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(reverse_normalized_hgvs_genomic.posedit))):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base) 
		alt_start = start - 1 #
		# Recover sequences
		ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),alt_start,end-1)
		ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		# Assemble  
		pos = start
		ref = ref_seq
		alt = ref_seq + ins_seq						

	# Substitutions
	elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
		ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
		alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
		pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

	# Deletions
	elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2
		start = start - 1
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,start)
		# Assemble  
		pos = str(start)
		ref = pre_base + hgvs_del_seq
		alt = pre_base

	
	# inv
	elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):						
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start -1
		start = start
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = ins_seq
		if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
			my_seq = Seq(vcf_del_seq)
			alt = str(my_seq.reverse_complement()) 
	
	# Delins
	elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(reverse_normalized_hgvs_genomic.posedit))):
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base -1)
		adj_start = start -1
		start = start
		try:
			ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
		except:
			ins_seq = ''
		else:
			if str(ins_seq) == 'None':
				ins_seq = ''		
		# Recover sequences
		hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
		vcf_del_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_del_seq
		alt = vcf_del_seq[:1] + ins_seq			
	
	
	# Duplications								
	elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
		end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base) #
		start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
		adj_start = start - 2 #
		start = start - 1 #
		# Recover sequences
		dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
		vcf_ref_seq	= sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
		# Assemble  
		pos = str(start)
		ref = vcf_ref_seq
		alt = vcf_ref_seq + dup_seq
	else:
		chr = ''
		ref = ''
		alt = ''
		pos = ''
		
	# ADD SURROUNDING BASES
	if chr != '' and pos != '' and ref != '' and alt != '':
#       Add 2 post bases
# 		pos = int(pos)
# 		pre_end_pos = pos + len(ref)
# 		end_pos = pre_end_pos + 1
# 		post = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),pre_end_pos-1,end_pos)
# 		ref = ref + post
# 		alt = alt + post
		# Add another pre-base
		pre_pos = int(pos) - 1
		pre_pos 
		prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),pre_pos-1,pre_pos)
		pos = str(pre_pos)				
		ref = prev + ref
		alt = prev + alt
	
	# Dictionary the VCF
	vcf_dict = {'chr' : chr, 'pos' : pos, 'ref' : ref, 'alt' : alt}
	return vcf_dict			
