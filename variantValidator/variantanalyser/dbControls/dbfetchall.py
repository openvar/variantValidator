#from mysql.connector import MySQLConnection, Error
#from dbconfig import read_db_config
import dbConnection

def execute(query):
	# db_config = read_db_config()
	# try:
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
		# print('No Data...')
		rows = ['none', 'No data']
	#except Error as e:
		# print ('Connection failed.')
		#rows = ['error', e]
		# print e
	#finally:
	cursor.close()
	conn.close()
	# print ('Connection closed...')
	return rows

# Methods

def get_transcript_info_for_gene(gene_symbol):
	#""" Connect to MySQL database """
	query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated, IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_info WHERE hgncSymbol = '%s'" %(gene_symbol)
	rows = execute(query)
	return rows

def get_transcribed_span_for_gene(gene_symbol, primary_assembly):
	#""" Connect to MySQL database """
	table = 'genePos' + primary_assembly.replace('GRCh', '')
	query = "SELECT chr, start, end FROM " + table + " WHERE symbol = '%s'" %(gene_symbol) 
	rows = execute(query)
	return rows

def get_g_to_g_info():
	#""" Connect to MySQL database """
	query = "SELECT refSeqGeneID, refSeqChromosomeID, startPos, endPos, orientation, hgncSymbol, genomeBuild FROM refSeqGene_loci"
	table = execute(query)
	return table

def get_all_transcriptID():
	query = "SELECT refSeqID FROM transcript_info"
	table = execute(query)
	return table	

if __name__ == '__main__':
	query_with_fetchall()
	
	
