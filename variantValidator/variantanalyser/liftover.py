# -*- coding: utf-8 -*-
"""
liftover.py
 
Can be used to liftover genomic positions between builds 
"""

from pyliftover import LiftOver
import re
import copy
import hgvs
import hgvs.posedit
import hgvs.edit

# Assembly data
current_assemblies = {
	"hg19" : "GRCh37_p13",
	"hg38" : "GRCh38_p2"
	}

def chr_id(chr, assembly):
	# Available genome builds
	NCBI36 = {
			"1":"NC_000001.9",
			"2":"NC_000002.10",
			"3":"NC_000003.10",
			"4":"NC_000004.10",
			"5":"NC_000005.8",
			"6":"NC_000006.10",
			"7":"NC_000007.12",
			"8":"NC_000008.9",
			"9":"NC_000009.10",
			"10":"NC_000010.9",
			"11":"NC_000011.8",
			"12":"NC_000012.10",
			"13":"NC_000013.9",
			"14":"NC_000014.7",
			"15":"NC_000015.8",
			"16":"NC_000016.8",
			"17":"NC_000017.9",
			"18":"NC_000018.8",
			"19":"NC_000019.8",
			"20":"NC_000020.9",
			"21":"NC_000021.7",
			"22":"NC_000022.9",
			"23":"NC_000023.9",
			"24":"NC_000024.8",
			"x":"NC_000023.9",
			"y":"NC_000024.8",
			"X":"NC_000023.9",
			"Y":"NC_000024.8"	
			}

	GRCh37_p13 = {
				"1" : "NC_000001.10",
				"2" : "NC_000002.11",
				"3" : "NC_000003.11",
				"4" : "NC_000004.11",
				"5" : "NC_000005.9",
				"6" : "NC_000006.11",
				"7" : "NC_000007.13",
				"8" : "NC_000008.10",
				"9" : "NC_000009.11",
				"10" : "NC_000010.10",
				"11" : "NC_000011.9",
				"12" : "NC_000012.11",
				"13" : "NC_000013.10",
				"14" : "NC_000014.8",
				"15" : "NC_000015.9",
				"16" : "NC_000016.9",
				"17" : "NC_000017.10",
				"18" : "NC_000018.9",
				"19" : "NC_000019.9",
				"20" : "NC_000020.10",
				"21" : "NC_000021.8",
				"22" : "NC_000022.10",
				"23" : "NC_000023.10",
				"24" : "NC_000024.9",	
				"x" : "NC_000023.10",
				"y" : "NC_000024.9",
				"X" : "NC_000023.10",
				"Y" : "NC_000024.9"
				}

	GRCh38_p2 = {
				"1" : "NC_000001.11",
				"2" : "NC_000002.12",
				"3" : "NC_000003.12",
				"4" : "NC_000004.12",
				"5" : "NC_000005.10",
				"6" : "NC_000006.12",
				"7" : "NC_000007.14",
				"8" : "NC_000008.11",
				"9" : "NC_000009.12",
				"10" : "NC_000010.11",
				"11" : "NC_000011.10",
				"12" : "NC_000012.12",
				"13" : "NC_000013.11",
				"14" : "NC_000014.9",
				"15" : "NC_000015.10",
				"16" : "NC_000016.10",
				"17" : "NC_000017.11",
				"18" : "NC_000018.10",
				"19" : "NC_000019.10",
				"20" : "NC_000020.11",
				"21" : "NC_000021.9",
				"22" : "NC_000022.11",
				"23" : "NC_000023.11",
				"24" : "NC_000024.10",
				"x" : "NC_000023.11",
				"y" : "NC_000024.10",
				"X" : "NC_000023.11",
				"Y" : "NC_000024.10"			
    			}
	# Copy the call data to rs
	rs = copy.deepcopy(chr) 
	# extract the chromosome number string
	chromosome = chr
	# Convert call line to rs line
	if assembly == 'GRCh37_p13':
		ac = GRCh37_p13.get(chromosome)
	if assembly == 'GRCh38_p2':
		ac = GRCh38_p2.get(chromosome)
	if assembly == 'NCBI36':
		ac = NCBI36.get(chromosome)
	# Return ac 
	return ac

# HGVS variant object liftover
def hgvs_liftover(hgvs_variant, primary_assembly, target_assembly, hn):
	# Note hn will be set to Target assembly. This is as it should be
	
	#print primary_assembly
	#print target_assembly
	
	# Convert primary_assembly names
	if primary_assembly == 'NCBI36':
		primary_assembly = 'hg18'
	if primary_assembly == 'GRCh37':
		primary_assembly = 'hg19'
	if primary_assembly == 'GRCh38':
		primary_assembly = 'hg38'
	# Convert target_assembly names
	if target_assembly == 'NCBI36':
		target_assembly = 'hg18'
	if target_assembly == 'GRCh37':
		target_assembly = 'hg19'
	if target_assembly == 'GRCh38':
		target_assembly = 'hg38'
		
	# Import chain log
	lo = LiftOver(primary_assembly, target_assembly)
	
	# extract chromosome number and convert to UCSC chr tagged number
	accession = str(hgvs_variant.ac)
	unversion = accession.split('.')[0]
	chr = unversion[-2:]
	lead0 = re.compile('^0')
	if lead0.search(chr):
		chr = chr.replace('0','')
	else:
		pass
	chr = 'chr' + chr
	
	# check for ranges and assign position 1 and 2
	pos_1 = 'NA'
	pos_2 = 'NA'
	range = re.compile('_')
	if range.search(str(hgvs_variant.posedit.pos)):
		pos_1 = str(hgvs_variant.posedit.pos.start.base)
		pos_2 = str(hgvs_variant.posedit.pos.end.base)
	else:
		pos_1 = str(hgvs_variant.posedit.pos)
		
	# Get converted positions	
	lo_1 = ''
	lo_2 = ''
	lo_1 = str(lo.convert_coordinate(chr, int(pos_1)))
	if pos_2 != 'NA':
		lo_2 = str(lo.convert_coordinate(chr, int(pos_2)))
	else:
		lo_2 = lo_1
	lo_1 = lo_1.replace(',', '\t')
	lo_2 = lo_2.replace(',', '\t')
	lo_1 = lo_1.split()
	lo_2 = lo_2.split()
	#print str(lo_1)
	#print str(lo_2)
	lo_1 = lo_1[1]
	lo_2 = lo_2[1]
	
	# Get converted accession number
	chr_no = chr.replace('chr', '')
	current_target_assembly = current_assemblies.get(target_assembly)
	target_ac = chr_id(chr_no, current_target_assembly)
	# Convert hgvs object to new variant description	
	target_hgvs_variant = copy.deepcopy(hgvs_variant)
	target_hgvs_variant.ac = target_ac
	target_hgvs_variant.posedit.pos.start.base = int(lo_1)
	target_hgvs_variant.posedit.pos.end.base = int(lo_2)
	target_hgvs_variant = hn.normalize(target_hgvs_variant)
	return target_hgvs_variant		
		
# <LICENSE>

# </LICENSE>		
		
		
		