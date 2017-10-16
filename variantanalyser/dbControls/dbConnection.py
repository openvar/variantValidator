# -*- coding: utf-8 -*-

"""
dbconnect.py

Connects to MySQL and returns a connection pool
"""

import mysql.connector
from mysql.connector.pooling import MySQLConnectionPool
import dbconfig
from dbconfig import read_db_config
import os

_connection_pool = None

def get_connection():
    global _connection_pool
    if not _connection_pool:
  		VALIDATOR_DB_URL = os.environ.get('VALIDATOR_DB_URL')
  		if VALIDATOR_DB_URL is not None:			
  			configurations = VALIDATOR_DB_URL.replace('mysqlx://', '')
  			user_pass,host_database = configurations.split('@')
  			user,password = user_pass.split(':')
  			host,database = host_database.split('/')
 			db_config = {
			   'user': user,
			   'password': password,
			   'host': host,
			   'database': database,
			   'raise_on_warnings': True,
			 }  			
  		else:	
  			db_config = read_db_config()
		_connection_pool = mysql.connector.pooling.MySQLConnectionPool(pool_size=10, **db_config) # MySQLConnection(**db_config)
    return _connection_pool

__all__ = [ 'getConnection' ]

# <LICENSE>

# </LICENSE>