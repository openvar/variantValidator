# from mysql.connector import MySQLConnection, Error
# from dbconfig import read_db_config
import os
import dbConnection

# Set up os paths data and log folders
ROOT = os.path.dirname(os.path.abspath(__file__))

def update(entry, data, table):	
	# MySQL queries
	if table == 'genePos37' or table == 'genePos38':
		query = "UPDATE " + table +  " SET symbol=%s, name=%s, prevSymbol=%s, reference=%s, assembly=%s, chr=%s, start=%s, end=%s, refSeqTranscriptID=%s, refSeqGeneID=%s, updated=NOW()   WHERE hgncID = %s"
			
		conn = dbConnection.get_connection().get_connection()
		cursor = conn.cursor()
		cursor.execute(query, (data['symbol'],
			data['name'],
			data['prevSymbol'],
			data['reference'],
			data['assembly'],
			data['chr'],
			data['start'],
			data['end'],
			data['refSeqTranscriptID'],
			data['refSeqGeneID'],
			data['hgncID']))
		success = 'true'
		conn.commit()			
		
	if table == 'transcript_id':
		accession = entry
		desc = data
		query = "UPDATE transcript_id SET description = %s, updated = NOW() WHERE accession = %s"
	
		conn = dbConnection.get_connection().get_connection()
		cursor = conn.cursor()
		cursor.execute(query, (desc, accession))
		success = 'true'
		conn.commit()

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
	
#except Error as e:
#	error = str(e)
#	success = e		
#	me = open(os.path.join(ROOT, 'mysql_error.txt'), 'a')
#	me.write(error + '\n')
#	me.write(str(query) + '\n')		
#	me.close()		
	# print error

#finally:
	cursor.close()
	conn.close()
	# conn.close()
	return success
	

def update_refSeqGene_data(rsg_data):	
	# Configure database
	# db_config = read_db_config()
	
	# query = "UPDATE refSeqGene_loci SET refSeqChromosomeID='%s', genomeBuild='%s', startPos='%s', endPos='%s', orientation='%s', totalLength='%s', chrPos='%s', rsgPos='%s', entrezID='%s', hgncSymbol='%s', updated=NOW() WHERE refSeqGeneID='%s'" %(rsg_data[1], rsg_data[2], rsg_data[3], rsg_data[4], rsg_data[5], rsg_data[6], rsg_data[7], rsg_data[8], rsg_data[9], rsg_data[10], rsg_data[0])
	# query = "UPDATE refSeqGene_loci SET startPos=%s, endPos=%s, orientation=%s, totalLength=%s, chrPos=%s, rsgPos=%s, entrezID=%s, hgncSymbol=%s, updated=NOW() WHERE refSeqGeneID=%s AND refSeqChromosomeID='%s' AND genomeBuild='%s'"
	query = "UPDATE refSeqGene_loci SET hgncSymbol=%s, updated=NOW() WHERE refSeqGeneID=%s"

	# print query
#try:
	# print('Connecting to MySQL database...')
	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor()
	# cursor.execute(query, (rsg_data[3], rsg_data[4], rsg_data[5], rsg_data[6], rsg_data[7], rsg_data[8], rsg_data[9], rsg_data[10], rsg_data[0], rsg_data[1], rsg_data[2]))
	cursor.execute(query, (rsg_data[10], rsg_data[0]))
	success = 'true'
	conn.commit()
	
#except Error as e:
#	error = str(e)
#	success = e		
#	me = open(os.path.join(ROOT, 'mysql_error.txt'), 'a')
#	me.write(error + '\n')
#	me.write(str(query) + '\n')		
#	me.close()		
	# print error

#finally:
	cursor.close()
	conn.close()
	# conn.close()
	return success

def update_transcript_loci(update_data, primary_assembly):
	# Configure database
	# db_config = read_db_config()
	data = update_data
	data['name'] = str(data['name']).replace("'", "\'")
	table = 'genePos' + str(primary_assembly.replace('GRCh', ''))
	query = "UPDATE " + table +  " SET hgncID=%s, symbol=%s, name=%s, prevSymbol=%s, reference=%s, assembly=%s, chr=%s, start=%s, end=%s, refSeqGeneID=%s, updated=NOW() WHERE refSeqTranscriptID=%s" 
#try:
	# print('Connecting to MySQL database...')
	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor()
	cursor.execute(query, (data['hgncID'],data['symbol'],data['name'],data['prevSymbol'],data['reference'],data['assembly'],data['chr'],data['start'],data['end'],data['refSeqGeneID'],data['refSeqTranscriptID']))
	success = 'true'
	conn.commit()
	
#except Error as e:
#	error = str(e)
#	success = e		
#	me = open(os.path.join(ROOT, 'mysql_error.txt'), 'a')
#	me.write(error + '\n')
#	me.write(str(query) + '\n')		
#	me.close()		
	# print error

#finally:
	cursor.close()
	conn.close()
	# conn.close()
	return success

if __name__ == '__main__':
	update()
