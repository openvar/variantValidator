import mysql.connector
from mysql.connector.pooling import MySQLConnectionPool


class Mixin:
    """
    A mixin containing the database initialisation routines.
    """
    def __init__(self, db_config):
        self.conn = None
        # self.cursor will be none UNLESS you're wrapping a function in @handleCursor, which automatically opens and
        # closes connections for you.
        self.cursor = None
        self.dbConfig = db_config

        self.pool = mysql.connector.pooling.MySQLConnectionPool(pool_size=10, connect_timeout=1209600, **self.dbConfig)
        self.conn = self.pool.get_connection()

    def __del__(self):
        if self.conn.is_connected():
            try:
                self.conn.close()
            except mysql.connector.errors.NotSupportedError:
                pass
            self.conn = None
        if self.pool:
            self.pool = None

# <LICENSE>
# Copyright (C) 2019 VariantValidator Contributors
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
