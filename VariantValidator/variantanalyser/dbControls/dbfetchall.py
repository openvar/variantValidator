# -*- coding: utf-8 -*-

"""
dbfetchall.py

Functions which make MySQL fetchall queries
"""

import dbConnection

def execute(query):
	conn = dbConnection.get_connection().get_connection()# MySQLConnection(**db_config)
	cursor = conn.cursor()
	cursor.execute(query)
	# Commit query
	rows = []
	rows = cursor.fetchall()
	# if rows is not None:
	if rows != []:
		pass
	else:
		rows = ['none', 'No data']
	cursor.close()
	conn.close()
	return rows

# Methods

def get_transcript_info_for_gene(gene_symbol):
	query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated, IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_info WHERE hgncSymbol = '%s'" %(gene_symbol)
	rows = execute(query)
	return rows

"""
Marked for removal
"""
# def get_transcribed_span_for_gene(gene_symbol, primary_assembly):
# 	table = 'genePos' + primary_assembly.replace('GRCh', '')
# 	query = "SELECT chr, start, end FROM " + table + " WHERE symbol = '%s'" %(gene_symbol) 
# 	rows = execute(query)
# 	return rows

def get_g_to_g_info():
	query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol, genomeBuild FROM refSeqGene_loci"
	table = execute(query)
	return table

def get_all_transcriptID():
	query = "SELECT refSeqID FROM transcript_info"
	table = execute(query)
	return table	

if __name__ == '__main__':
	query_with_fetchall()
	
# <LICENSE>

# </LICENSE>	
