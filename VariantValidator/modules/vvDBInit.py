import logging
import random
import re
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
    from mysql.connector.pooling import (
        MySQLConnectionPool,
    )
except ModuleNotFoundError:
    MySQLConnectionPool = None


class Mixin:
    """
    A mixin containing the database initialisation routines.

    Database write capability is detected when the connection pool is
    initialised. The result is stored on the object so write operations
    can safely operate in read-only deployments.
    """

    WRITE_PRIVILEGES = (
        "INSERT",
        "UPDATE",
        "DELETE",
    )

    def __init__(self, db_config):
        self.pool = None
        self.dbConfig = db_config

        self.write_permissions = {
            privilege: False
            for privilege in self.WRITE_PRIVILEGES
        }

        self.write_enabled = False

        self.init_db()

    def __del__(self):
        if getattr(
            self,
            "pool",
            None,
        ):
            self.pool = None

    def init_db(self):
        """
        Initialise MySQL or MariaDB connection pool.

        NOTE:
        - We keep default unicode behaviour for ALL VV modules.
        - The dbSNP loader overrides cursor behaviour only for itself.
        - Write permissions are detected after the pool is created.
        """

        # Prefer MySQL Connector/Python when available.
        if MySQLConnectionPool is not None:
            self.pool = MySQLConnectionPool(
                pool_size=5,
                **self.dbConfig,
            )

            self._detect_write_permissions()

            return

        # Otherwise fall back to MariaDB.
        if MariaDBConnectionPool is not None:
            pool_kwargs = {
                "pool_size": 5,
                "pool_reset_connection": False,
                "host": self.dbConfig["host"],
                "user": self.dbConfig["user"],
                "port": int(
                    self.dbConfig["port"],
                ),
                "password": self.dbConfig["password"],
                "database": self.dbConfig["database"],
            }

            try:
                self.pool = MariaDBConnectionPool(
                    pool_name=(
                        f"pool{random.random()}"
                    ),
                    **pool_kwargs,
                )
            except MariaDBProgrammingError:
                # Retry with a different pool name.
                self.pool = MariaDBConnectionPool(
                    pool_name=(
                        f"pool{random.random()}"
                    ),
                    **pool_kwargs,
                )

            self._detect_write_permissions()

            return

        raise ModuleNotFoundError(
            "Neither mysql.connector nor mariadb is installed."
        )

    def _detect_write_permissions(self):
        """
        Determine whether the current database account can perform
        INSERT, UPDATE and DELETE operations on the configured database.

        Detection is based entirely on MySQL/MariaDB grants. No write is
        attempted.

        The result is stored in:

            self.write_permissions

        and:

            self.write_enabled

        ``write_enabled`` is true only when all three normal VV DML
        privileges are available.
        """

        permissions = {
            privilege: False
            for privilege in self.WRITE_PRIVILEGES
        }

        conn = None
        cursor = None

        try:
            conn = self.pool.get_connection()

            try:
                conn.ping(
                    reconnect=True,
                    attempts=1,
                    delay=0,
                )
            except Exception:
                pass

            cursor = conn.cursor(
                buffered=True,
            )

            cursor.execute(
                "SHOW GRANTS"
            )

            grants = [
                row[0]
                for row in cursor.fetchall()
            ]

            permissions = (
                self._parse_write_permissions(
                    grants,
                    str(
                        self.dbConfig.get(
                            "database",
                            "",
                        )
                    ),
                )
            )

        except Exception:
            logger.exception(
                "Unable to determine database write "
                "permissions; assuming read-only mode"
            )

            permissions = {
                privilege: False
                for privilege in self.WRITE_PRIVILEGES
            }

        finally:
            if cursor is not None:
                try:
                    cursor.close()
                except Exception:
                    pass

            if conn is not None:
                try:
                    conn.close()
                except Exception:
                    pass

        self.write_permissions = permissions

        self.write_enabled = all(
            permissions.values()
        )

        if self.write_enabled:
            logger.info(
                "Database write access enabled for database %s",
                self.dbConfig.get(
                    "database",
                ),
            )
        else:
            disabled = [
                privilege
                for privilege, enabled
                in permissions.items()
                if not enabled
            ]

            logger.warning(
                "Database write access is restricted for "
                "database %s; disabled privileges: %s. "
                "VariantValidator will operate in read-only "
                "database mode where required.",
                self.dbConfig.get(
                    "database",
                ),
                ", ".join(disabled),
            )

    @classmethod
    def _parse_write_permissions(
        cls,
        grants,
        database,
    ):
        """
        Parse SHOW GRANTS output and determine the effective INSERT,
        UPDATE and DELETE capabilities for the configured database.

        Supports:

        - global grants such as ON *.*
        - database grants such as ON `validator`.*
        - partial revokes such as:

              REVOKE INSERT, UPDATE
              ON `validator`.*
              FROM `root`@`localhost`

        The latter is important for MySQL accounts that retain global
        privileges but are explicitly restricted for the VV database.
        """

        database_normalised = (
            database.strip(
                "`"
            ).lower()
        )

        permissions = {
            privilege: False
            for privilege in cls.WRITE_PRIVILEGES
        }

        partial_revokes = {
            privilege: False
            for privilege in cls.WRITE_PRIVILEGES
        }

        database_grants = {
            privilege: False
            for privilege in cls.WRITE_PRIVILEGES
        }

        global_grants = {
            privilege: False
            for privilege in cls.WRITE_PRIVILEGES
        }

        for grant in grants:
            if not isinstance(
                grant,
                str,
            ):
                continue

            grant_text = grant.strip()

            if not grant_text:
                continue

            upper = grant_text.upper()

            # ----------------------------------------------------------
            # Explicit partial REVOKE
            # ----------------------------------------------------------

            if upper.startswith(
                "REVOKE "
            ):
                if not cls._grant_applies_to_database(
                    grant_text,
                    database_normalised,
                ):
                    continue

                privileges = (
                    cls._extract_privileges(
                        grant_text,
                        "REVOKE",
                    )
                )

                for privilege in privileges:
                    if (
                        privilege
                        in partial_revokes
                    ):
                        partial_revokes[
                            privilege
                        ] = True

                continue

            # ----------------------------------------------------------
            # GRANT
            # ----------------------------------------------------------

            if not upper.startswith(
                "GRANT "
            ):
                continue

            scope = (
                cls._extract_grant_scope(
                    grant_text,
                )
            )

            if scope is None:
                continue

            privileges = (
                cls._extract_privileges(
                    grant_text,
                    "GRANT",
                )
            )

            if scope == "GLOBAL":
                for privilege in privileges:
                    if (
                        privilege
                        in global_grants
                    ):
                        global_grants[
                            privilege
                        ] = True

            elif (
                scope
                == database_normalised
            ):
                for privilege in privileges:
                    if (
                        privilege
                        in database_grants
                    ):
                        database_grants[
                            privilege
                        ] = True

        # --------------------------------------------------------------
        # Database-level grants are direct permission.
        # Global grants apply unless a partial revoke overrides them.
        # --------------------------------------------------------------

        for privilege in cls.WRITE_PRIVILEGES:
            if database_grants[
                privilege
            ]:
                permissions[
                    privilege
                ] = True
                continue

            if (
                global_grants[
                    privilege
                ]
                and not partial_revokes[
                    privilege
                ]
            ):
                permissions[
                    privilege
                ] = True

        return permissions

    @staticmethod
    def _extract_privileges(
        grant_text,
        statement,
    ):
        """
        Extract explicit privileges from a GRANT/REVOKE statement.

        ``ALL PRIVILEGES`` is treated as providing the normal DML
        privileges used by VariantValidator.
        """

        upper = grant_text.upper()

        if statement == "GRANT":
            body = re.sub(
                r"^GRANT\s+",
                "",
                grant_text,
                count=1,
                flags=re.IGNORECASE,
            )
        else:
            body = re.sub(
                r"^REVOKE\s+",
                "",
                grant_text,
                count=1,
                flags=re.IGNORECASE,
            )

        body = re.split(
            r"\s+ON\s+",
            body,
            maxsplit=1,
            flags=re.IGNORECASE,
        )[0]

        if "ALL PRIVILEGES" in upper:
            return [
                "INSERT",
                "UPDATE",
                "DELETE",
            ]

        privileges = []

        for privilege in Mixin.WRITE_PRIVILEGES:
            pattern = (
                rf"\b{privilege}\b"
            )

            if re.search(
                pattern,
                body,
                flags=re.IGNORECASE,
            ):
                privileges.append(
                    privilege
                )

        return privileges

    @staticmethod
    def _extract_grant_scope(
        grant_text,
    ):
        """
        Return:

            GLOBAL
            database name

        for the grant target.

        Unsupported scopes are returned as None.

        GRANT statements use:

            ON ... TO ...

        while REVOKE statements use:

            ON ... FROM ...
        """

        match = re.search(
            r"\bON\s+(.+?)\s+(?:TO|FROM)\s+",
            grant_text,
            flags=re.IGNORECASE,
        )

        if not match:
            return None

        target = (
            match.group(1)
            .strip()
        )

        if target == "*.*":
            return "GLOBAL"

        database_match = re.fullmatch(
            r"`?([^`.*]+)`?\.\*",
            target,
        )

        if database_match:
            return (
                database_match.group(1)
                .strip()
                .lower()
            )

        return None

    @staticmethod
    def _grant_applies_to_database(
        grant_text,
        database,
    ):
        """
        Return True when a GRANT or REVOKE explicitly applies to the
        configured database.
        """

        scope = (
            Mixin._extract_grant_scope(
                grant_text,
            )
        )

        if scope is None:
            return False

        if scope == "GLOBAL":
            return False

        return (
            scope.lower()
            == database.lower()
        )

    def can_write(
        self,
        privilege="INSERT",
    ):
        """
        Return whether the configured database account currently has
        the requested write privilege.
        """

        privilege = (
            privilege.upper()
        )

        if privilege not in (
            self.write_permissions
        ):
            return False

        return bool(
            self.write_permissions[
                privilege
            ]
        )

    def get_conn(self):
        """
        Get a live connection from the pool with retry + backoff.

        Handles:
        - stale pooled connections
        - MySQL timeouts
        - transient network issues
        """
        delays = (
            0,
            0.5,
            2,
            5,
        )

        last_exception = None

        for delay in delays:
            if delay:
                time.sleep(
                    delay
                )

            try:
                conn = (
                    self.pool.get_connection()
                )

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
                    "Database connection failed "
                    "health check; retrying with "
                    "fresh connection"
                )

        raise last_exception

    def get_cursor(
        self,
        conn,
    ):
        """
        Default VV cursor:
        - Unicode/UTF-8
        - Buffered

        dbSNP loader will override this independently.
        """
        try:
            cursor = conn.cursor(
                buffered=True,
            )
        except Exception:
            self.init_db()

            conn = self.get_conn()

            cursor = conn.cursor(
                buffered=True,
            )

        return cursor


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License version 3 (or at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
