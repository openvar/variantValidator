import logging
import random
import time


logger = logging.getLogger(__name__)


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
        Get a live connection from the pool with retry + backoff.

        Handles:
        - stale pooled connections
        - MySQL timeouts
        - transient network issues
        """
        delays = (0, 0.5, 2, 5)
        last_exception = None

        for delay in delays:
            if delay:
                time.sleep(delay)

            try:
                conn = self.pool.get_connection()

                # Ensure the connection is alive.
                try:
                    conn.ping(
                        reconnect=True,
                        attempts=1,
                        delay=0,
                    )
                except Exception:
                    conn.close()
                    raise

                return conn

            except Exception as e:
                last_exception = e

                # Rebuild pool on failure.
                self.init_db()

                logger.exception(
                    "Database connection failed health check; "
                    "retrying with fresh connection"
                )

        raise last_exception

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


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
