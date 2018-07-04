# -*- coding: utf-8 -*-

"""
dbupdate.py

Functions which make MySQL update statements
"""
import os
import dbConnection

# Set up os paths data and log folders
ROOT = os.path.dirname(os.path.abspath(__file__))

def update(entry, data, table):	
	# MySQL queries

# 	if table == 'genePos37' or table == 'genePos38':
# 		query = "UPDATE " + table +  " SET symbol=%s, name=%s, prevSymbol=%s, reference=%s, assembly=%s, chr=%s, start=%s, end=%s, refSeqTranscriptID=%s, refSeqGeneID=%s, updated=NOW()   WHERE hgncID = %s"	
# 		conn = dbConnection.get_connection().get_connection()
# 		cursor = conn.cursor()
# 		cursor.execute(query, (data['symbol'],
# 			data['name'],
# 			data['prevSymbol'],
# 			data['reference'],
# 			data['assembly'],
# 			data['chr'],
# 			data['start'],
# 			data['end'],
# 			data['refSeqTranscriptID'],
# 			data['refSeqGeneID'],
# 			data['hgncID']))
# 		success = 'true'
# 		conn.commit()			
# 		
# 	if table == 'transcript_id':
# 		accession = entry
# 		desc = data
# 		query = "UPDATE transcript_id SET description = %s, updated = NOW() WHERE accession = %s"
# 		conn = dbConnection.get_connection().get_connection()
# 		cursor = conn.cursor()
# 		cursor.execute(query, (desc, accession))
# 		success = 'true'
# 		conn.commit()

	if table == 'transcript_info':	
		accession = entry
		description = data[1]
		variant = data[2]
		version = data[3]
		hgnc_symbol = data[4]
		uta_symbol	= data[5]
		query = "UPDATE transcript_info SET description=%s, transcriptVariant=%s, currentVersion=%s, hgncSymbol=%s, utaSymbol=%s, updated=NOW() WHERE refSeqID = %s"
		conn = dbConnection.get_connection().get_connection()
		cursor = conn.cursor()
		cursor.execute(query, (description, variant, version, hgnc_symbol, uta_symbol, accession))
		success = 'true'
		conn.commit()
	cursor.close()
	conn.close()
	return success
	

def update_refSeqGene_data(rsg_data):	
	query = "UPDATE refSeqGene_loci SET hgncSymbol=%s, updated=NOW() WHERE refSeqGeneID=%s"
	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor()
	cursor.execute(query, (rsg_data[10], rsg_data[0]))
	success = 'true'
	conn.commit()
	cursor.close()
	conn.close()
	return success

"""
mark for removal
"""
# def update_transcript_loci(update_data, primary_assembly):
# 	data = update_data
# 	data['name'] = str(data['name']).replace("'", "\'")
# 	table = 'genePos' + str(primary_assembly.replace('GRCh', ''))
# 	query = "UPDATE " + table +  " SET hgncID=%s, symbol=%s, name=%s, prevSymbol=%s, reference=%s, assembly=%s, chr=%s, start=%s, end=%s, refSeqGeneID=%s, updated=NOW() WHERE refSeqTranscriptID=%s" 
# 	conn = dbConnection.get_connection().get_connection()
# 	cursor = conn.cursor()
# 	cursor.execute(query, (data['hgncID'],data['symbol'],data['name'],data['prevSymbol'],data['reference'],data['assembly'],data['chr'],data['start'],data['end'],data['refSeqGeneID'],data['refSeqTranscriptID']))
# 	success = 'true'
# 	conn.commit()
# 	cursor.close()
# 	conn.close()
# 	return success

if __name__ == '__main__':
	update()
# <LICENSE>

# </LICENSE>