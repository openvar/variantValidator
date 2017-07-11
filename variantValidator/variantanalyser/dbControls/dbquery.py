import dbConnection

def query_with_fetchone(entry, table):
	# """ Connect to MySQL database """
	
	# MySQL queries
	if table == 'genePos37' or table == 'genePos38':
		import re
		if re.search('HGNC:', entry):			
			query = "SELECT hgncID, symbol, name, prevSymbol, reference, assembly, chr, start, end, refSeqTranscriptID, refSeqGeneID FROM " + table + " WHERE hgncID = '%s'" %(entry)
		if re.search('sym:', entry):
			symbol = entry.replace('sym:', '')
			query = "SELECT hgncID, symbol, name, prevSymbol, reference, assembly, chr, start, end, refSeqTranscriptID, refSeqGeneID FROM " + table + " WHERE symbol LIKE '%s' OR prevSymbol LIKE '|%s|'" %(symbol, symbol)
	if table == 'transcript_id':
		query = "SELECT accession, description, updated, IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_id WHERE accession = '%s'" %(entry)
	
	if table == 'transcript_info':
		query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated, IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_info WHERE refSeqID = '%s'" %(entry)

	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor(buffered=True)
	cursor.execute(query)

	# Blank list for row
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
	

if __name__ == '__main__':
	query_with_fetchone()
