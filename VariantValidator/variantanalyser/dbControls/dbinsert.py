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

def insert_LRG_protein_data(lrg_p, rs_p):
	query = "INSERT INTO LRG_proteins(LRGproteinID, RefSeqProteinID) VALUES(%s,%s)"
	conn = dbConnection.get_connection().get_connection()
	cursor = conn.cursor()
	cursor.execute(query, (lrg_p, rs_p))			
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
# Copyright (C) 2018  Peter Causey-Freeman, University of Leicester
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>