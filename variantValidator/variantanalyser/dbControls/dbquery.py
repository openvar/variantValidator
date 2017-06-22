#from mysql.connector import MySQLConnection, Error
#from dbconfig import read_db_config
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
		# print query
	# configure the database
	# db_config = read_db_config()
	# print 'Database Query...'
#	try:
	# print('Connecting to MySQL database...')
	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor(buffered=True)
	cursor.execute(query)

	# print ('Connection established.')
	# Blank list for row
	row = []
	row = cursor.fetchone()

	if row is not None:
		pass
	else:
		# print('No Data...')
		row = ['none', 'No data']
		
#	except Error as e:
#		# print ('Connection failed.')
#		row = ['error', e]
		# print e

#	finally:
	#conn.close()
	cursor.close()
	conn.close()
	# print ('Connection closed...')
	# print str(row)
	return row
	

if __name__ == '__main__':
	query_with_fetchone()
