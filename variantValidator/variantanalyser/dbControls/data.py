# -*- coding: utf-8 -*-

# import database modules
import dbquery
import dbinsert
import dbupdate
import dbfetchone
import dbfetchall

# Import python modules
import os
import sys
import re

# Needs functions from variantanalyser - directory above, unless in a single directory
try:
	import variantanalyser.functions
except ImportError:	
	parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	os.sys.path.insert(0,parentdir)
	import functions	

# Retrieve transcript information
def in_entries(entry, table): 
	
	# Use dbquery.py to connect to mysql and return the necessary data
	# import dbquery
	
	# GeneNamesLoci
	if table == 'genePos37' or table == 'genePos38':
		row = dbquery.query_with_fetchone(entry.split('.')[0], table)
		
		if row[0] == 'error':
			data = {
				'error' : 'false',
				'description': 'false'
				}
		
			data['error'] = row[0]
			data['description'] = row[1]
	
		elif row[0] == 'none':
			data = {
				'none' : 'false',
				'description': 'false'
				}
		
			data['none'] = row[0]
			data['description'] = row[1]
	
		else:
			data = {}
			data['hgncID'] = row[0]
			data['symbol'] = row[1]
			data['name'] = row[2]
			data['prevSymbol'] = row[3]
			data['reference'] = row[4]
			data['assembly'] = row[5]
			data['chr'] = row[6]
			data['start'] = row[7]
			data['end'] = row[8]
			data['refSeqTranscriptID'] = row[9]
			data['refSeqGeneID'] = row[10]
		
	# Transcript ID
	if table == 'transcript_id':
		row = dbquery.query_with_fetchone(entry.split('.')[0], table)
	
		if row[0] == 'error':
			data = {
				'error' : 'false',
				'description': 'false'
				}
		
			data['error'] = row[0]
			data['description'] = row[1]
	
		elif row[0] == 'none':
			data = {
				'none' : 'false',
				'description': 'false'
				}
		
			data['none'] = row[0]
			data['description'] = row[1]
	
		else:
			data = {
				'accession' : 'false',
				'description': 'false',
				'updated' : 'false',
				'expiry' : 'false'
				}
		
			data['accession'] = row[0]
			data['description'] = row[1]
			data['updated'] = row[2]
			data['expiry'] = row[3]

	if table == 'transcript_info':
		row = dbquery.query_with_fetchone(entry.split('.')[0], table)
	
		if row[0] == 'error':
			data = {
				'error' : 'false',
				'description': 'false'
				}
		
			data['error'] = row[0]
			data['description'] = row[1]
	
		elif row[0] == 'none':
			data = {
				'none' : 'false',
				'description': 'false'
				}
		
			data['none'] = row[0]
			data['description'] = row[1]
	
		else:
			data = {
				'accession' : 'false',
				'description': 'false',
				'updated' : 'false',
				'expiry' : 'false'
				}
		
			data['accession'] = row[0]
			data['description'] = row[1]
			data['variant'] = row[2]
			data['version'] = row[3]
			data['hgnc_symbol'] = row[4]
			data['uta_symbol'] = row[5]
			data['updated'] = row[6]
			data['expiry'] = row[7]
	
	return data		

# Add new entry  	
def add_entry(entry, data, table):
	# import dbinsert
	success = dbinsert.insert(entry, data, table)
	return success
    
def insert_transcript_loci(add_data, primary_assembly):
	# import dbinsert
	success = dbinsert.insert_transcript_loci(add_data, primary_assembly)
	return success


# Update entries
def update_entry(entry, data, table):
	# import dbupdate
	success = dbupdate.update(entry, data, table)
	return success

def update_transcript_info_record(accession, hdp):
	accession = accession.split('.')[0] # list[3].split('.')[0]
		
# 	import from 1 level above
# 	import os
# 	import sys
# 	parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 	os.sys.path.insert(0,parentdir)
# 	import functions	
	
	# Search Entrez for corresponding record
	record = functions.entrez_efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
	version = record.id 
	description = record.description
	variant = '0'
	
	# import re
	if re.search('transcript variant', description):				
		tv = re.search('transcript variant \w+', description)
		tv = str(tv.group(0))
		tv = tv.replace('transcript variant', '')
		variant = tv.strip()
		variant = variant.upper() # Some tv descriptions are a or A
	else:
		variant = '0'
	
	# Get information from UTA
	try:
		uta_info = hdp.get_tx_identity_info(version)
	except:
		version_ac_ver = version.split('.')
		version = version_ac_ver[0] + '.' + str(int(version_ac_ver[1]) - 1)
		uta_info = hdp.get_tx_identity_info(version)
	
	uta_symbol = str(uta_info[6])	
	# Search hgnc rest to see if symbol is out of date
	rest_data = functions.hgnc_rest(path = "/search/prev_symbol/" + uta_symbol)
	# If the name is correct no record will be found
	if rest_data['error'] == 'false':
		if int(rest_data['record']['response']['numFound']) == 0:
			hgnc_symbol = uta_info[6]
		else:
			hgnc_symbol = rest_data['record']['response']['docs'][0]['symbol']
	else:
		hgnc_symbol = 'unassigned'
			
	# Query information
	query_info = [accession, description, variant, version, hgnc_symbol, uta_symbol]
	table='transcript_info'

	# Update the transcript_info table (needs plugging in) 
	returned_data = in_entries(accession, table)
	# If the entry is not in the database add it
	if 'none' in returned_data:
		add_entry(accession, query_info, table)
	# If the data in the entry has changed, update it
	else:
		update_entry(accession, query_info, table)
	return			

def update_refSeqGene_loci(rsg_data):
	# First query the database
	# import dbfetchone
	entry_exists = dbfetchone.get_refSeqGene_data_by_refSeqGeneID(rsg_data[0], rsg_data[2])
 	if entry_exists[0] == 'none':
 		# import dbinsert
 		dbinsert.insert_refSeqGene_data(rsg_data)
	else:
		# import dbupdate
		dbupdate.update_refSeqGene_data(rsg_data)
	return

def update__transcript_loci(update_data, primary_assembly):
	# import dbupdate
	success = dbupdate.update_transcript_loci(update_data, primary_assembly)
	return success	

def update_lrg_rs_lookup(lrg_rs_lookup):
	# First query the database
	# import dbfetchone
	rsgID = dbfetchone.get_RefSeqGeneID_from_lrgID(lrg_rs_lookup[0])
 	if rsgID == 'none':
 		# import dbinsert
 		dbinsert.insert_RefSeqGeneID_from_lrgID(lrg_rs_lookup)
		return
	else:
		return

def update_lrgt_rst(lrgtx_to_rstID):
	# First query the database
	# import dbfetchone
	rstID = dbfetchone.get_RefSeqTranscriptID_from_lrgTranscriptID(lrgtx_to_rstID[0])
 	if rstID == 'none':
 		# import dbinsert
 		dbinsert.insert_LRG_transcript_data(lrgtx_to_rstID)
		return
	else:
		return

# Direct methods (GET)
def get_transcript_info_for_gene(gene_symbol):
	# import dbfetchall
	rows = dbfetchall.get_transcript_info_for_gene(gene_symbol)
	return rows
	
def get_uta_symbol(gene_symbol):
	# returns the UTA gene symbol when HGNC gene symbol is input
	# import dbfetchone
	utaSymbol = str(dbfetchone.get_utaSymbol(gene_symbol)[0])
	return utaSymbol

def get_hgnc_symbol(gene_symbol):
	# returns the HGNC gene symbol when UTA gene symbol is input
	# import dbfetchone
	hgncSymbol = str(dbfetchone.get_hgncSymbol(gene_symbol)[0])
	return hgncSymbol

def get_current_hgnc_symbol(gene_symbol, primary_assembly):
	# returns current HGNC gene symbol when previous gene symbol is input
	# import dbfetchone
	hgncSymbol = dbfetchone.get_current_hgnc_symbol(gene_symbol, primary_assembly)
	return hgncSymbol
	
def get_transcript_description(transcript_id):
	# returns the transcript description for a given transcript
	# import dbfetchone
	tx_description = dbfetchone.get_transcript_description(transcript_id)			
	return tx_description

def get_gene_symbol_from_transcriptID(transcript_id):
	# returns gene symbol for a given transcript ID
	# import dbfetchone
	gene_symbol = dbfetchone.get_gene_symbol_from_transcriptID(transcript_id)
	return gene_symbol

def get_gene_symbol_from_refSeqGeneID(refSeqGeneID):
	# Returns the databases most up-to-date gene symbol for a given NG_ ID
	# import dbfetchone
	gene_symbol = dbfetchone.get_gene_symbol_from_refSeqGeneID(refSeqGeneID)
	return gene_symbol

def get_transcribed_span_for_transcript(transcript_id, primary_assembly):
	# Returns the chromosome and span for the selected transcript
	# import dbfetchone
	span = dbfetchone.get_transcribed_span_for_transcript(transcript_id, primary_assembly)
	return span

def get_transcribed_span_for_gene(gene_symbol, primary_assembly):
	# Returns the chromosome and span for the selected gene_symbol
	# import dbfetchall
	span = dbfetchall.get_transcribed_span_for_gene(gene_symbol, primary_assembly)
	span_out = {}
	chr = span[0][0]
	start = 0
	end = 0
	for row in span:
		if start == 0:
			start = row[1]
		if end == 0:
			end = row[2]
		if row[1] < start:
			start = row[1]
		if row[2] > end:
			end = row[2]			
	span_out['chr'] = chr
	span_out['start'] = start
	span_out['end'] = end						
	return span_out	

def get_g_to_g_info():
	# Recovers the g_to_g data table
	# import dbfetchall
	table = dbfetchall.get_g_to_g_info()
	return table

def get_all_transcriptID():
	# Returns a list of transcript IDs in our database
	# import dbfetchall	
	table = dbfetchall.get_all_transcriptID()
	return table

def get_RefSeqGeneID_from_lrgID(lrgID):
	# Get the relevant RefSeqGeneID for a given LRG ID
	# import dbfetchone
	rsgID = dbfetchone.get_RefSeqGeneID_from_lrgID(lrgID)
	return rsgID

def get_RefSeqTranscriptID_from_lrgTranscriptID(lrg_txID):
	# import dbfetchone
	rstID = dbfetchone.get_RefSeqTranscriptID_from_lrgTranscriptID(lrg_txID)
	return rstID

def get_lrgTranscriptID_from_RefSeqTranscriptID(rstID):
	# import dbfetchone
	lrg_tx = dbfetchone.get_lrgTranscriptID_from_RefSeqTranscriptID(rstID)
	return lrg_tx
	
def get_lrgID_from_RefSeqGeneID(rsgID):	
	# import dbfetchone
	lrgID = dbfetchone.get_lrgID_from_RefSeqGeneID(rsgID)
	return lrgID	

def get_refseqgene_info(refseqgene_id, primary_assembly):
	# import dbfetchone
	refseqgene_info = dbfetchone.get_refseqgene_info(refseqgene_id, primary_assembly)
	return refseqgene_info

