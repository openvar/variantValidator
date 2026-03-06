"""Unit tests for the SQLite database backend.

These tests use a temporary SQLite file — no MySQL, no network required.
"""
import importlib
import importlib.util
import sqlite3
import os
import sys
import pytest


def _load_vvdbinit():
    """Import vvDBInit directly, bypassing VariantValidator/__init__.py
    which would pull in vvhgvs and other heavy optional dependencies."""
    spec = importlib.util.spec_from_file_location(
        "vvDBInit",
        os.path.join(os.path.dirname(__file__), "..", "VariantValidator", "modules", "vvDBInit.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def vvdbinit_module():
    return _load_vvdbinit()


@pytest.fixture
def sqlite_db_path(tmp_path):
    """Return path to a temporary SQLite database file."""
    return str(tmp_path / "test_vvdb.sqlite")


@pytest.fixture
def sqlite_mixin(sqlite_db_path, vvdbinit_module):
    """Instantiate SQLiteDBInit with a temp database."""
    return vvdbinit_module.SQLiteDBInit(sqlite_db_path)


def test_sqlite_db_init_creates_file(sqlite_db_path, vvdbinit_module):
    """SQLiteDBInit must create the SQLite file if it does not exist."""
    assert not os.path.exists(sqlite_db_path)
    vvdbinit_module.SQLiteDBInit(sqlite_db_path)
    assert os.path.exists(sqlite_db_path)


def test_sqlite_db_init_creates_tables(sqlite_mixin):
    """All required tables must be created on init."""
    conn = sqlite3.connect(sqlite_mixin.db_path)
    cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = {row[0] for row in cursor.fetchall()}
    conn.close()
    expected = {
        'transcript_info', 'refSeqGene_loci', 'LRG_RSG_lookup',
        'LRG_transcripts', 'LRG_proteins', 'stableGeneIds', 'version'
    }
    assert tables == expected, f"Unexpected tables: {tables ^ expected}"


def test_sqlite_get_conn_returns_connection(sqlite_mixin):
    """get_conn() must return a usable sqlite3.Connection."""
    conn = sqlite_mixin.get_conn()
    assert conn is not None
    conn.close()


def test_sqlite_thread_safety(sqlite_db_path, vvdbinit_module):
    """Connection opened with check_same_thread=False must be usable from another thread."""
    import threading

    db = vvdbinit_module.SQLiteDBInit(sqlite_db_path)
    conn = db.get_conn()

    result = []
    error = []

    def worker():
        try:
            cursor = conn.cursor()
            cursor.execute("SELECT 1")
            result.append(cursor.fetchone()[0])
        except Exception as e:
            error.append(str(e))

    t = threading.Thread(target=worker)
    t.start()
    t.join()

    conn.close()
    assert not error, f"Cross-thread connection error: {error}"
    assert result == [1]
