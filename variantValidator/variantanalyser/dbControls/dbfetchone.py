# -*- coding: utf-8 -*-

"""
dbfetchone.py

Functions which make MySQL fetchone queries
"""

import dbConnection

def execute(query):	
	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor(buffered=True)
	cursor.execute(query)

	row = []
	row = cursor.fetchone()

	if row is not None:
		pass
	else:
		# print('No Data...')
		row = ['none', 'No data']
	cursor.close()
	conn.close()
	return row

# Methods
def get_utaSymbol(gene_symbol):
	query = "SELECT utaSymbol FROM transcript_info WHERE hgncSymbol = '%s'" %(gene_symbol)
	row = execute(query)
	return row
	
def get_hgncSymbol(gene_symbol):
	query = "SELECT hgncSymbol FROM transcript_info WHERE utaSymbol = '%s'" %(gene_symbol)
	row = execute(query)
	return row

"""
marked for removal
"""	
# def get_current_hgnc_symbol(gene_symbol, primary_assembly):
# 	strip_assembly = primary_assembly.replace('GRCh', '')
# 	query = "SELECT symbol FROM genePos%s WHERE symbol LIKE '%s' OR prevSymbol LIKE '%s'" %(strip_assembly, gene_symbol, gene_symbol)
# 	symbol = str(execute(query)[0])
# 	hgncSymbol = symbol.replace('|', '')
# 	return hgncSymbol

def get_transcript_description(transcript_id):
	transcript_id = transcript_id.split('.')[0]
	query = "SELECT description FROM transcript_info WHERE refSeqID = '%s'" %(transcript_id)
	tx_description = str(execute(query)[0])
	return tx_description	

def get_gene_symbol_from_transcriptID(transcript_id):
	transcript_id = transcript_id.split('.')[0]
	query = "SELECT hgncSymbol FROM transcript_info WHERE refSeqID = '%s'" %(transcript_id)
	gene_symbol = str(execute(query)[0])
	return gene_symbol

def get_refSeqGene_data_by_refSeqGeneID(refSeqGeneID, genomeBuild):
	query = "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, orientation, totalLength, chrPos, rsgPos, entrezID, hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = '%s' AND genomeBuild = '%s'" %(refSeqGeneID, genomeBuild)
	refSeqGene_data = execute(query)
	return refSeqGene_data

def get_gene_symbol_from_refSeqGeneID(refSeqGeneID):
	refseqgene_id = refSeqGeneID #.split('.')[0]
	query = "SELECT hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = '%s'" %(refseqgene_id)
	gene_symbol = str(execute(query)[0])
	return gene_symbol	

"""
marked for removal
"""
# def get_transcribed_span_for_transcript(transcript_id, primary_assembly):
# 	transcript_id = transcript_id.split('.')[0]
# 	table_identifier = primary_assembly.replace('GRCh', '')
# 	table = 'genePos' + table_identifier
# 	query = "SELECT chr, start, end FROM " + table + " WHERE refSeqTranscriptID = '%s'" %(transcript_id)
# 	span = execute(query)
# 	return span
	
def get_RefSeqGeneID_from_lrgID(lrgID):
	query = "SELECT RefSeqGeneID FROM LRG_RSG_lookup WHERE lrgID = '%s'" %(lrgID)
	rsgID = execute(query)
	rsgID = rsgID[0]
	return rsgID 
	
def get_RefSeqTranscriptID_from_lrgTranscriptID(lrgtxID):
	query = "SELECT RefSeqTranscriptID FROM LRG_transcripts WHERE LRGtranscriptID = '%s'" %(lrgtxID)
	rstID = execute(query)
	rstID = rstID[0]
	return rstID
	
def	get_lrgTranscriptID_from_RefSeqTranscriptID(rstID):		
	query = "SELECT LRGtranscriptID FROM LRG_transcripts WHERE RefSeqTranscriptID = '%s'" %(rstID)
	lrg_tx = execute(query)
	lrg_tx = lrg_tx[0]
	return lrg_tx

def get_lrgID_from_RefSeqGeneID(rsgID):
	query = "SELECT lrgID, status FROM LRG_RSG_lookup WHERE RefSeqGeneID = '%s'" %(rsgID)	
	lrgID = execute(query)
	lrgID = lrgID
	return lrgID
	
def get_refseqgene_info(refseqgene_id, primary_assembly):
	query = "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos FROM refSeqGene_loci WHERE refSeqGeneID = '%s' AND genomeBuild = '%s'" %(refseqgene_id, primary_assembly)	
	refseqgene_info = execute(query)
	return refseqgene_info	
	
def get_RefSeqProteinID_from_lrgProteinID(lrg_p):
	query = "SELECT RefSeqProteinID FROM LRG_proteins WHERE LRGproteinID = '%s'" %(lrg_p)
	rspID = execute(query)
	rspID = rspID[0]
	return rspID

def get_lrgProteinID_from_RefSeqProteinID(rs_p):
	query = "SELECT LRGproteinID FROM LRG_proteins WHERE  RefSeqProteinID = '%s'" %(rs_p)
	lrpID = execute(query)
	lrpID = lrpID[0]
	return lrpID	
	
if __name__ == '__main__':
	query_with_fetchone()
	
# <LICENSE>

# </LICENSE>