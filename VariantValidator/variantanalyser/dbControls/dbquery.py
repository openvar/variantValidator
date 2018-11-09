# -*- coding: utf-8 -*-

"""
dbquery.py

Functions which fetch data from transcript_info
"""

import dbConnection

def query_with_fetchone(entry, table):
	# """ Connect to MySQL database """
	
	# MySQL queries
	
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
