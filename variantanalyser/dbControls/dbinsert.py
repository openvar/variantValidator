# -*- coding: utf-8 -*-

"""
dbinsert.py

Functions which make MySQL insert statements
"""

import os
import dbConnection

# Set up os paths data and log folders
ROOT = os.path.dirname(os.path.abspath(__file__))

def insert(entry, data, table):	
	conn = dbConnection.get_connection().get_connection() # MySQLConnection(**db_config)
	cursor = conn.cursor()
	# MySQL queries
	
# 	if table == 'genePos37' or table == 'genePos38':
# 		query = "INSERT INTO " +  table + "(hgncID, symbol, name, prevSymbol, reference, assembly, chr, start, end, refSeqTranscriptID, refSeqGeneID, updated) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, NOW())"
# 		cursor.execute(query, (data['hgncID'], 
# 			data['symbol'],
# 			data['name'],
# 			data['prevSymbol'],
# 			data['reference'],
# 			data['assembly'],
# 			data['chr'],
# 			data['start'],
# 			data['end'],
# 			data['refSeqTranscriptID'],
# 			data['refSeqGeneID']
# 			))
	
# 	if table == 'transcript_id':
# 		accession = entry
# 		desc = data
# 		query = "INSERT INTO transcript_id(accession, description, updated) VALUES(%s,%s,NOW())" 
# 		cursor.execute(query, (accession, desc))

	if table == 'transcript_info':	
		accession = entry
		description = data[1]
		variant = data[2]
		version = data[3]
		hgnc_symbol = data[4]
		uta_symbol	= data[5]
		query = "INSERT INTO transcript_info(refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated) VALUES(%s,%s, %s, %s, %s, %s, NOW())" 
		cursor.execute(query, (accession, description, variant, version, hgnc_symbol, uta_symbol))
		
	# Query report
	if cursor.lastrowid:
		success = 'true'
	else:
		success = 'Unknown error'

	# Commit and close connection
	conn.commit()
	cursor.close()
	conn.close()
	return success
	
def insert_refSeqGene_data(rsg_data):
	query = "INSERT INTO refSeqGene_loci(refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, orientation, totalLength, chrPos, rsgPos, entrezID, hgncSymbol, updated) VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW())"
	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor()
	cursor.execute(query, (rsg_data[0], rsg_data[1], rsg_data[2], rsg_data[3], rsg_data[4], rsg_data[5], rsg_data[6], rsg_data[7], rsg_data[8], rsg_data[9], rsg_data[10]))
	# Query report
	if cursor.lastrowid:
		success = 'true'
	else:
		success = 'Unknown error'

	# Commit and close connection
	conn.commit()
	cursor.close()
	conn.close()
	return success

"""
marked for removal
"""
# def insert_transcript_loci(add_data, primary_assembly):
# 	data = add_data
# 	table = 'genePos' + str(primary_assembly.replace('GRCh', ''))
# 	query = "INSERT INTO " +  table + "(hgncID, symbol, name, prevSymbol, reference, assembly, chr, start, end, refSeqTranscriptID, refSeqGeneID, updated) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, NOW())"
# 	conn = dbConnection.get_connection().get_connection()
# 	cursor = conn.cursor()
# 	cursor.execute(query, (data['hgncID'], data['symbol'], data['name'], data['prevSymbol'], data['reference'], data['assembly'], data['chr'], data['start'], data['end'], data['refSeqTranscriptID'], data['refSeqGeneID']))			
# 	# Query report
# 	if cursor.lastrowid:
# 		success = 'true'
# 	else:
# 		success = 'Unknown error'
# 
# 	# Commit and close connection
# 	conn.commit()
# 	cursor.close()
# 	conn.close()
# 	return success

def insert_RefSeqGeneID_from_lrgID(lrg_rs_lookup):
	query = "INSERT INTO LRG_RSG_lookup(lrgID, hgncSymbol, RefSeqGeneID, status) VALUES(%s,%s,%s,%s)"
	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor()
	cursor.execute(query, (lrg_rs_lookup[0], lrg_rs_lookup[1], lrg_rs_lookup[2], lrg_rs_lookup[3]))			
	# Query report
	if cursor.lastrowid:
		success = 'true'
	else:
		success = 'Unknown error'

	# Commit and close connection
	conn.commit()
	cursor.close()
	conn.close()
	return success	
	
def insert_LRG_transcript_data(lrgtx_to_rstID):
	query = "INSERT INTO LRG_transcripts(LRGtranscriptID, RefSeqTranscriptID) VALUES(%s,%s)"
	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor()
	cursor.execute(query, (lrgtx_to_rstID[0], lrgtx_to_rstID[1]))			
	# Query report
	if cursor.lastrowid:
		success = 'true'
	else:
		success = 'Unknown error'

	# Commit and close connection
	conn.commit()
	cursor.close()
	conn.close()
	return success		

if __name__ == '__main__':
	insert()

# <LICENSE>

# </LICENSE>