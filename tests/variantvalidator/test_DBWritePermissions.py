from __future__ import annotations

import os
from typing import Any

import pytest
from unittest.mock import MagicMock

from VariantValidator.modules.vvDBInit import (
    Mixin,
)

from VariantValidator.modules.vvDBInsert import (
    Mixin as DBInsertMixin,
)


# ----------------------------------------------------------------------
# Live database configuration
#
# Defaults match the developer workstation configuration.
#
# A deployment can override these without changing this test file:
#
#   VV_TEST_DB_HOST
#   VV_TEST_DB_PORT
#   VV_TEST_DB_NAME
#   VV_TEST_DB_USER
#   VV_TEST_DB_PASSWORD
# ----------------------------------------------------------------------

DB_HOST = os.getenv(
    "VV_TEST_DB_HOST",
    "127.0.0.1",
)

DB_PORT = int(
    os.getenv(
        "VV_TEST_DB_PORT",
        "3306",
    )
)

DB_NAME = os.getenv(
    "VV_TEST_DB_NAME",
    "validator",
)

DB_USER = os.getenv(
    "VV_TEST_DB_USER",
    "root",
)

DB_PASSWORD = os.getenv(
    "VV_TEST_DB_PASSWORD",
    "",
)


# ----------------------------------------------------------------------
# Grant parser tests
# ----------------------------------------------------------------------

def test_parse_global_write_privileges() -> None:
    grants = [
        (
            "GRANT SELECT, INSERT, UPDATE, DELETE "
            "ON *.* TO `root`@`localhost` "
            "WITH GRANT OPTION"
        ),
    ]

    permissions = (
        Mixin._parse_write_permissions(
            grants,
            "validator",
        )
    )

    assert permissions == {
        "INSERT": True,
        "UPDATE": True,
        "DELETE": True,
    }


def test_parse_schema_write_privileges() -> None:
    grants = [
        (
            "GRANT SELECT, INSERT, UPDATE, DELETE "
            "ON `validator`.* "
            "TO `root`@`localhost` "
            "WITH GRANT OPTION"
        ),
    ]

    permissions = (
        Mixin._parse_write_permissions(
            grants,
            "validator",
        )
    )

    assert permissions == {
        "INSERT": True,
        "UPDATE": True,
        "DELETE": True,
    }


def test_partial_revoke_disables_schema_writes() -> None:
    grants = [
        (
            "GRANT SELECT, INSERT, UPDATE, DELETE "
            "ON *.* TO `root`@`localhost` "
            "WITH GRANT OPTION"
        ),
        (
            "REVOKE INSERT, UPDATE, DELETE "
            "ON `validator`.* "
            "FROM `root`@`localhost`"
        ),
    ]

    permissions = (
        Mixin._parse_write_permissions(
            grants,
            "validator",
        )
    )

    assert permissions == {
        "INSERT": False,
        "UPDATE": False,
        "DELETE": False,
    }


def test_partial_revoke_can_disable_only_one_privilege() -> None:
    grants = [
        (
            "GRANT SELECT, INSERT, UPDATE, DELETE "
            "ON *.* TO `root`@`localhost` "
            "WITH GRANT OPTION"
        ),
        (
            "REVOKE DELETE "
            "ON `validator`.* "
            "FROM `root`@`localhost`"
        ),
    ]

    permissions = (
        Mixin._parse_write_permissions(
            grants,
            "validator",
        )
    )

    assert permissions == {
        "INSERT": True,
        "UPDATE": True,
        "DELETE": False,
    }


def test_select_only_account_is_read_only() -> None:
    grants = [
        (
            "GRANT SELECT "
            "ON `validator`.* "
            "TO `vv`@`localhost`"
        ),
    ]

    permissions = (
        Mixin._parse_write_permissions(
            grants,
            "validator",
        )
    )

    assert permissions == {
        "INSERT": False,
        "UPDATE": False,
        "DELETE": False,
    }


def test_all_privileges_allow_database_writes() -> None:
    grants = [
        (
            "GRANT ALL PRIVILEGES "
            "ON `validator`.* "
            "TO `root`@`localhost` "
            "WITH GRANT OPTION"
        ),
    ]

    permissions = (
        Mixin._parse_write_permissions(
            grants,
            "validator",
        )
    )

    assert permissions == {
        "INSERT": True,
        "UPDATE": True,
        "DELETE": True,
    }


def test_global_all_privileges_allow_database_writes() -> None:
    grants = [
        (
            "GRANT ALL PRIVILEGES "
            "ON *.* "
            "TO `root`@`localhost` "
            "WITH GRANT OPTION"
        ),
    ]

    permissions = (
        Mixin._parse_write_permissions(
            grants,
            "validator",
        )
    )

    assert permissions == {
        "INSERT": True,
        "UPDATE": True,
        "DELETE": True,
    }


def test_unrelated_database_revoke_does_not_affect_validator() -> None:
    grants = [
        (
            "GRANT SELECT, INSERT, UPDATE, DELETE "
            "ON *.* TO `root`@`localhost` "
            "WITH GRANT OPTION"
        ),
        (
            "REVOKE INSERT, UPDATE, DELETE "
            "ON `other_database`.* "
            "FROM `root`@`localhost`"
        ),
    ]

    permissions = (
        Mixin._parse_write_permissions(
            grants,
            "validator",
        )
    )

    assert permissions == {
        "INSERT": True,
        "UPDATE": True,
        "DELETE": True,
    }


def test_empty_grants_are_read_only() -> None:
    permissions = (
        Mixin._parse_write_permissions(
            [],
            "validator",
        )
    )

    assert permissions == {
        "INSERT": False,
        "UPDATE": False,
        "DELETE": False,
    }


# ----------------------------------------------------------------------
# can_write tests
# ----------------------------------------------------------------------

def test_can_write_uses_detected_permissions() -> None:
    db = Mixin.__new__(Mixin)

    db.write_permissions = {
        "INSERT": True,
        "UPDATE": False,
        "DELETE": True,
    }

    assert (
        db.can_write("INSERT")
        is True
    )

    assert (
        db.can_write("UPDATE")
        is False
    )

    assert (
        db.can_write("DELETE")
        is True
    )


def test_can_write_unknown_privilege_is_false() -> None:
    db = Mixin.__new__(Mixin)

    db.write_permissions = {
        "INSERT": True,
        "UPDATE": True,
        "DELETE": True,
    }

    assert (
        db.can_write("CREATE")
        is False
    )


# ----------------------------------------------------------------------
# DBInsert read-only behaviour
# ----------------------------------------------------------------------

@pytest.mark.parametrize(
    (
        "method_name",
        "args",
    ),
    [
        (
            "insert",
            (
                "NM_000001.1",
                [
                    None,
                    "description",
                    "NM_000001.1:c.1A>G",
                    "1",
                    "GENE1",
                    "UTA1",
                ],
                "transcript_info",
            ),
        ),
        (
            "insert_refseq_gene_data",
            (
                [
                    "NG_000001.1",
                    "NC_000001.11",
                    "GRCh38",
                    "1",
                    "10",
                    "+",
                    "10",
                    "1",
                    "1",
                    "1",
                    "GENE1",
                ],
            ),
        ),
        (
            "insert_refseq_gene_id_from_lrg_id",
            (
                [
                    "LRG_1",
                    "GENE1",
                    "NG_000001.1",
                    "public",
                ],
            ),
        ),
        (
            "insert_lrg_transcript_data",
            (
                [
                    "LRG_1t1",
                    "NM_000001.1",
                ],
            ),
        ),
        (
            "insert_lrg_protein_data",
            (
                "LRG_1p1",
                "NP_000001.1",
            ),
        ),
        (
            "insert_gene_stable_ids",
            (
                {
                    "hgnc_id": "HGNC:1",
                    "hgnc_symbol": "GENE1",
                    "entrez_id": "1",
                    "ensembl_gene_id": "ENSG000001",
                    "omim_id": "1",
                    "ucsc_id": "uc001",
                    "vega_id": "OTTHUMG000001",
                    "ccds_id": "CCDS1",
                },
            ),
        ),
        (
            "update",
            (
                "NM_000001.1",
                [
                    None,
                    "description",
                    "NM_000001.1:c.1A>G",
                    "1",
                    "GENE1",
                    "UTA1",
                ],
            ),
        ),
        (
            "update_refseq_gene_data",
            (
                [
                    "NG_000001.1",
                    "NC_000001.11",
                    "GRCh38",
                    "1",
                    "10",
                    "+",
                    "10",
                    "1",
                    "1",
                    "1",
                    "GENE1",
                ],
            ),
        ),
        (
            "update_gene_stable_ids",
            (
                {
                    "hgnc_id": "HGNC:1",
                    "hgnc_symbol": "GENE1",
                    "entrez_id": "1",
                    "ensembl_gene_id": "ENSG000001",
                    "omim_id": "1",
                    "ucsc_id": "uc001",
                    "vega_id": "OTTHUMG000001",
                    "ccds_id": "CCDS1",
                },
            ),
        ),
        (
            "update_db_version",
            (
                "vvdb_2026_7",
            ),
        ),
    ],
)
def test_dbinsert_write_methods_skip_when_read_only(
    method_name: str,
    args: tuple[Any, ...],
) -> None:
    db = DBInsertMixin.__new__(
        DBInsertMixin,
    )

    db.dbConfig = {
        "database": "validator",
    }

    db.write_permissions = {
        "INSERT": False,
        "UPDATE": False,
        "DELETE": False,
    }

    db.pool = MagicMock()

    method = getattr(
        db,
        method_name,
    )

    result = method(
        *args,
    )

    assert result == "false"

    db.pool.get_connection.assert_not_called()


def test_dbinsert_insert_executes_when_insert_is_allowed() -> None:
    db = DBInsertMixin.__new__(
        DBInsertMixin,
    )

    db.dbConfig = {
        "database": "validator",
    }

    db.write_permissions = {
        "INSERT": True,
        "UPDATE": False,
        "DELETE": False,
    }

    conn = MagicMock()
    cursor = MagicMock()

    cursor.lastrowid = 123

    db.get_conn = MagicMock(
        return_value=conn,
    )

    db.get_cursor = MagicMock(
        return_value=cursor,
    )

    result = db.insert(
        "NM_000001.1",
        [
            None,
            "description",
            "NM_000001.1:c.1A>G",
            "1",
            "GENE1",
            "UTA1",
        ],
        "transcript_info",
    )

    assert result == "true"

    db.get_conn.assert_called_once()
    db.get_cursor.assert_called_once_with(
        conn,
    )

    cursor.execute.assert_called_once()

    conn.commit.assert_called_once()

    cursor.close.assert_called_once()
    conn.close.assert_called_once()


def test_dbinsert_update_executes_when_update_is_allowed() -> None:
    db = DBInsertMixin.__new__(
        DBInsertMixin,
    )

    db.dbConfig = {
        "database": "validator",
    }

    db.write_permissions = {
        "INSERT": False,
        "UPDATE": True,
        "DELETE": False,
    }

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(
        return_value=conn,
    )

    db.get_cursor = MagicMock(
        return_value=cursor,
    )

    result = db.update(
        "NM_000001.1",
        [
            None,
            "description",
            "NM_000001.1:c.1A>G",
            "1",
            "GENE1",
            "UTA1",
        ],
    )

    assert result == "true"

    db.get_conn.assert_called_once()
    db.get_cursor.assert_called_once_with(
        conn,
    )

    cursor.execute.assert_called_once()

    conn.commit.assert_called_once()

    cursor.close.assert_called_once()
    conn.close.assert_called_once()


def test_dbinsert_insert_does_not_report_success_when_no_rowid() -> None:
    db = DBInsertMixin.__new__(
        DBInsertMixin,
    )

    db.dbConfig = {
        "database": "validator",
    }

    db.write_permissions = {
        "INSERT": True,
        "UPDATE": False,
        "DELETE": False,
    }

    conn = MagicMock()
    cursor = MagicMock()

    cursor.lastrowid = None

    db.get_conn = MagicMock(
        return_value=conn,
    )

    db.get_cursor = MagicMock(
        return_value=cursor,
    )

    result = db.insert(
        "NM_000001.1",
        [
            None,
            "description",
            "NM_000001.1:c.1A>G",
            "1",
            "GENE1",
            "UTA1",
        ],
        "transcript_info",
    )

    assert result == "Unknown error"

    cursor.execute.assert_called_once()
    conn.commit.assert_called_once()


# ----------------------------------------------------------------------
# Live database detection
# ----------------------------------------------------------------------

def live_db() -> Mixin:
    """
    Connect to the database using the same basic configuration used by
    the local VariantValidator installation.

    A connection failure skips the live capability tests rather than
    turning an unavailable database into a false application failure.
    """

    config: dict[str, Any] = {
        "host": DB_HOST,
        "port": DB_PORT,
        "database": DB_NAME,
        "user": DB_USER,
        "password": DB_PASSWORD,
    }

    try:
        return Mixin(config)
    except Exception as exc:
        pytest.skip(
            "Live database capability test skipped: "
            "unable to connect to "
            f"{DB_USER}@{DB_HOST}:{DB_PORT}/{DB_NAME}: "
            f"{exc}"
        )


def test_live_database_reports_write_mode() -> None:
    db = live_db()

    mode = (
        "WRITABLE"
        if db.write_enabled
        else "READ-ONLY"
    )

    print(
        "\n"
        "============================================================\n"
        "VariantValidator database write capability\n"
        "============================================================\n"
        f"Host       : {DB_HOST}\n"
        f"Port       : {DB_PORT}\n"
        f"Database   : {DB_NAME}\n"
        f"User       : {DB_USER}\n"
        f"Mode       : {mode}\n"
        f"INSERT     : {db.can_write('INSERT')}\n"
        f"UPDATE     : {db.can_write('UPDATE')}\n"
        f"DELETE     : {db.can_write('DELETE')}\n"
        "============================================================"
    )

    assert isinstance(
        db.write_enabled,
        bool,
    )

    assert set(
        db.write_permissions
    ) == {
        "INSERT",
        "UPDATE",
        "DELETE",
    }


def test_live_writable_database() -> None:
    db = live_db()

    if not db.write_enabled:
        pytest.skip(
            "Database is currently READ-ONLY; "
            "writable-database branch is not applicable."
        )

    assert (
        db.can_write("INSERT")
        is True
    )

    assert (
        db.can_write("UPDATE")
        is True
    )

    assert (
        db.can_write("DELETE")
        is True
    )

    print(
        "\n"
        "LIVE DATABASE MODE: WRITABLE\n"
        "INSERT/UPDATE/DELETE are available."
    )


def test_live_read_only_database() -> None:
    db = live_db()

    if db.write_enabled:
        pytest.skip(
            "Database currently has WRITE access; "
            "read-only branch is not applicable."
        )

    assert (
        db.can_write("INSERT")
        is False
    )

    assert (
        db.can_write("UPDATE")
        is False
    )

    assert (
        db.can_write("DELETE")
        is False
    )

    print(
        "\n"
        "LIVE DATABASE MODE: READ-ONLY\n"
        "INSERT/UPDATE/DELETE are unavailable."
    )


def test_live_read_access_remains_available() -> None:
    db = live_db()

    conn = db.get_conn()

    try:
        cursor = db.get_cursor(
            conn,
        )

        try:
            cursor.execute(
                "SELECT 1"
            )

            result = cursor.fetchone()

        finally:
            cursor.close()

    finally:
        conn.close()

    assert result == (1,)

    print(
        "\n"
        "LIVE DATABASE READ ACCESS: OK"
    )

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 (or at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
