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


# ── Task 3: SQLiteDBGet ──────────────────────────────────────────────────────

def _load_vvdbget():
    """Import vvDBGet directly, bypassing VariantValidator/__init__.py.

    Relative imports (from . import vvDBInit) require the parent package to
    be present in sys.modules.  We register lightweight stubs for the package
    hierarchy and for heavy sibling modules we don't actually need.
    """
    import types

    repo_root = os.path.join(os.path.dirname(__file__), "..")
    modules_dir = os.path.join(repo_root, "VariantValidator", "modules")

    # ── 1. Register stub parent packages so relative imports work ──────────
    for pkg_name in ("VariantValidator", "VariantValidator.modules"):
        if pkg_name not in sys.modules:
            stub = types.ModuleType(pkg_name)
            stub.__path__ = [os.path.join(repo_root, *pkg_name.split("."))]
            stub.__package__ = pkg_name
            sys.modules[pkg_name] = stub

    def _load_file(name, filename):
        if name in sys.modules:
            return sys.modules[name]
        spec = importlib.util.spec_from_file_location(
            name,
            os.path.join(modules_dir, filename),
            submodule_search_locations=[],
        )
        mod = importlib.util.module_from_spec(spec)
        mod.__package__ = "VariantValidator.modules"
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        # Also expose as attribute on the parent package stub
        short_name = name.split(".")[-1]
        parent = sys.modules.get("VariantValidator.modules")
        if parent is not None:
            setattr(parent, short_name, mod)
        return mod

    # ── 2. Load vvDBInit (no heavy deps) ───────────────────────────────────
    vvdbinit_mod = _load_file("VariantValidator.modules.vvDBInit", "vvDBInit.py")

    # ── 3. Stub utils — only handleCursor is needed (it's a no-op wrapper) ─
    if "VariantValidator.modules.utils" not in sys.modules:
        import functools
        utils_stub = types.ModuleType("VariantValidator.modules.utils")
        utils_stub.__package__ = "VariantValidator.modules"

        def handleCursor(func):
            @functools.wraps(func)
            def wrapper(self, *args, **kwargs):
                return func(self, *args, **kwargs)
            return wrapper

        utils_stub.handleCursor = handleCursor
        sys.modules["VariantValidator.modules.utils"] = utils_stub
        setattr(sys.modules["VariantValidator.modules"], "utils", utils_stub)

    # ── 4. Load vvDBGet ─────────────────────────────────────────────────────
    return _load_file("VariantValidator.modules.vvDBGet", "vvDBGet.py")


@pytest.fixture(scope="session")
def vvdbget_module():
    return _load_vvdbget()


@pytest.fixture
def sqlite_get(sqlite_db_path, vvdbget_module):
    """SQLiteDBGet with an empty temp database."""
    return vvdbget_module.SQLiteDBGet(sqlite_db_path)


def test_get_uta_returns_none_when_empty(sqlite_get):
    result = sqlite_get.get_uta("BRCA1")
    assert result[0] == 'none'


def test_get_hgnc_returns_none_when_empty(sqlite_get):
    result = sqlite_get.get_hgnc("BRCA1")
    assert result[0] == 'none'


def test_execute_returns_none_when_empty(sqlite_get):
    result = sqlite_get.execute("SELECT refSeqID FROM transcript_info WHERE refSeqID = ?", ("NM_007294.3",))
    assert result == ['none', 'No data']


def test_execute_all_returns_none_when_empty(sqlite_get):
    result = sqlite_get.execute_all("SELECT refSeqID FROM transcript_info WHERE hgncSymbol = ?", ("BRCA1",))
    assert result == ['none', 'No data']


def test_get_transcript_info_for_gene_returns_empty_when_missing(sqlite_get):
    result = sqlite_get.get_transcript_info_for_gene("NOTAREALgene")
    assert result == ['none', 'No data']
