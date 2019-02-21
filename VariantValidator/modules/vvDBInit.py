import mysql.connector
from mysql.connector.pooling import MySQLConnectionPool
import os

class Mixin():
    '''
    A mixin containing the database initialisation routines.
    '''
    def __init__(self,val,dbConfig):
        self.conn = None
        # self.cursor will be none UNLESS you're wrapping a function in @handleCursor, which automatically opens and
        # closes connections for you.
        self.cursor=None
        self.dbConfig=dbConfig
        # Construct database URL
        #'mysqlx://vvadmin:var1ant@127.0.0.1/validator'
        self.path="mysqlx://"+dbConfig["user"]+":"+dbConfig["password"]+"@"+dbConfig["host"]+"/"+dbConfig["database"]
        os.environ["VALIDATOR_DB_URL"]=self.path
        self.val=val
        self.pool=mysql.connector.pooling.MySQLConnectionPool(pool_size=10, **self.dbConfig)
        self.conn=self.pool.get_connection()

    def __del__(self):
        if self.conn:
            self.conn.close()
            self.conn=None
        if self.pool:
            self.pool.close()
            self.pool=None
        if self.val:
            self.val=None
