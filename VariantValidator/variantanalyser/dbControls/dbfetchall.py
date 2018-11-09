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
