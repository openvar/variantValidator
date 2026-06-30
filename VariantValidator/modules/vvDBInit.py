import random

try:
    import mariadb
    MariaDBConnectionPool = mariadb.ConnectionPool
    MariaDBProgrammingError = mariadb.ProgrammingError
except ModuleNotFoundError:
    mariadb = None
    MariaDBConnectionPool = None

    class MariaDBProgrammingError(Exception):
        """Fallback exception when mariadb is unavailable."""
        pass

try:
    from mysql.connector.pooling import MySQLConnectionPool
except ModuleNotFoundError:
    MySQLConnectionPool = None


class Mixin:
    """
    A mixin containing the database initialisation routines.
    """

    def __init__(self, db_config):
        self.pool = None
        self.dbConfig = db_config
        self.init_db()

    def __del__(self):
        if getattr(self, "pool", None):
            self.pool = None

    def init_db(self):
        """
        Initialise MySQL or MariaDB connection pool.

        NOTE:
        - We keep default unicode behaviour for ALL VV modules.
        - The dbSNP loader overrides cursor behaviour only for itself.
        """

        # Prefer MySQL Connector/Python when available.
        if MySQLConnectionPool is not None:
            self.pool = MySQLConnectionPool(
                pool_size=5,
                **self.dbConfig,
            )
            return

        # Otherwise fall back to MariaDB.
        if MariaDBConnectionPool is not None:
            pool_kwargs = {
                "pool_size": 5,
                "pool_reset_connection": False,
                "host": self.dbConfig["host"],
                "user": self.dbConfig["user"],
                "port": int(self.dbConfig["port"]),
                "password": self.dbConfig["password"],
                "database": self.dbConfig["database"],
            }

            try:
                self.pool = MariaDBConnectionPool(
                    pool_name=f"pool{random.random()}",
                    **pool_kwargs,
                )

            except MariaDBProgrammingError:
                # Retry with a different pool name.
                self.pool = MariaDBConnectionPool(
                    pool_name=f"pool{random.random()}",
                    **pool_kwargs,
                )

            return

        raise ModuleNotFoundError(
            "Neither mysql.connector nor mariadb is installed."
        )

    def get_conn(self):
        """
        Get a connection from the pool.
        """
        try:
            conn = self.pool.get_connection()
        except Exception:
            self.init_db()
            conn = self.pool.get_connection()
        return conn

    def get_cursor(self, conn):
        """
        Default VV cursor:
        - Unicode/UTF-8
        - Buffered

        dbSNP loader will override this independently.
        """
        try:
            cursor = conn.cursor(buffered=True)
        except Exception:
            self.init_db()
            conn = self.get_conn()
            cursor = conn.cursor(buffered=True)
        return cursor

# <LICENSE>
# Copyright (C) 2016-2026 VariantValidator Contributors
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
