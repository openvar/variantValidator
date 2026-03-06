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


# ── Task 4: SQLiteDBInsert ───────────────────────────────────────────────────

def _load_vvdbinsert():
    """Import vvDBInsert directly, reusing the stubs set up by _load_vvdbget."""
    import types

    repo_root = os.path.join(os.path.dirname(__file__), "..")
    modules_dir = os.path.join(repo_root, "VariantValidator", "modules")

    # Ensure parent package stubs exist (idempotent — _load_vvdbget may have
    # already registered them in an earlier fixture call).
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
        short_name = name.split(".")[-1]
        parent = sys.modules.get("VariantValidator.modules")
        if parent is not None:
            setattr(parent, short_name, mod)
        return mod

    # Ensure utils stub exists
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

    # Load dependency chain
    _load_file("VariantValidator.modules.vvDBInit", "vvDBInit.py")
    _load_file("VariantValidator.modules.vvDBGet", "vvDBGet.py")
    return _load_file("VariantValidator.modules.vvDBInsert", "vvDBInsert.py")


@pytest.fixture(scope="session")
def vvdbinsert_module():
    return _load_vvdbinsert()


@pytest.fixture
def sqlite_insert(sqlite_db_path, vvdbinsert_module):
    return vvdbinsert_module.SQLiteDBInsert(sqlite_db_path)


# -- transcript_info ----------------------------------------------------------

def test_insert_transcript_info_roundtrip(sqlite_insert):
    """Insert a transcript row then verify it can be retrieved."""
    result = sqlite_insert.insert_transcript_info(
        refseq_id="NM_007294.3",
        description="BRCA1",
        transcript_variant="NM_007294",
        current_version="3",
        hgnc_symbol="BRCA1",
        uta_symbol="BRCA1",
    )
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    rows = conn.execute(
        "SELECT refSeqID FROM transcript_info WHERE refSeqID = ?",
        ("NM_007294.3",),
    ).fetchall()
    conn.close()
    assert len(rows) == 1
    assert rows[0][0] == "NM_007294.3"


def test_insert_via_insert_method(sqlite_insert):
    """The generic insert() method should also write transcript_info rows."""
    data = [None, "TP53 transcript", "NM_000546", "4", "TP53", "TP53"]
    result = sqlite_insert.insert("NM_000546.4", data, table="transcript_info")
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    rows = conn.execute(
        "SELECT refSeqID, hgncSymbol FROM transcript_info WHERE refSeqID = ?",
        ("NM_000546.4",),
    ).fetchall()
    conn.close()
    assert len(rows) == 1
    assert rows[0][1] == "TP53"


def test_update_transcript_info(sqlite_insert):
    """update() should change an existing transcript_info row."""
    # First insert
    sqlite_insert.insert_transcript_info(
        refseq_id="NM_001234.1",
        description="original desc",
        transcript_variant="NM_001234",
        current_version="1",
        hgnc_symbol="GENE1",
        uta_symbol="GENE1",
    )
    # Then update
    data = [None, "updated desc", "NM_001234_updated", "2", "GENE1_UPD", "GENE1_UPD"]
    result = sqlite_insert.update("NM_001234.1", data)
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    row = conn.execute(
        "SELECT description, currentVersion FROM transcript_info WHERE refSeqID = ?",
        ("NM_001234.1",),
    ).fetchone()
    conn.close()
    assert row[0] == "updated desc"
    assert row[1] == "2"


# -- refSeqGene_loci ----------------------------------------------------------

def test_insert_refseq_gene_data_roundtrip(sqlite_insert):
    """INSERT OR REPLACE into refSeqGene_loci."""
    rsg_data = [
        "NG_007524.2", "NC_000017.11", "GRCh38",
        "43044295", "43125370", "+", "81076",
        "chr17:43044295-43125370", "1-81076",
        "672", "BRCA1",
    ]
    result = sqlite_insert.insert_refseq_gene_data(rsg_data)
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    row = conn.execute(
        "SELECT refSeqGeneID, hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = ?",
        ("NG_007524.2",),
    ).fetchone()
    conn.close()
    assert row is not None
    assert row[1] == "BRCA1"


def test_update_refseq_gene_data(sqlite_insert):
    """update_refseq_gene_data() should update hgncSymbol."""
    rsg_data = [
        "NG_009999.1", "NC_000001.11", "GRCh38",
        "1000", "2000", "+", "1001",
        "chr1:1000-2000", "1-1001",
        "9999", "OLDSYM",
    ]
    sqlite_insert.insert_refseq_gene_data(rsg_data)

    rsg_data[10] = "NEWSYM"
    result = sqlite_insert.update_refseq_gene_data(rsg_data)
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    row = conn.execute(
        "SELECT hgncSymbol FROM refSeqGene_loci WHERE refSeqGeneID = ?",
        ("NG_009999.1",),
    ).fetchone()
    conn.close()
    assert row[0] == "NEWSYM"


# -- LRG tables ---------------------------------------------------------------

def test_insert_lrg_rsg_lookup_roundtrip(sqlite_insert):
    """INSERT OR REPLACE into LRG_RSG_lookup."""
    lrg_rs_lookup = ["LRG_1", "BRCA1", "NG_007524.2", "public"]
    result = sqlite_insert.insert_refseq_gene_id_from_lrg_id(lrg_rs_lookup)
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    row = conn.execute(
        "SELECT RefSeqGeneID FROM LRG_RSG_lookup WHERE lrgID = ?",
        ("LRG_1",),
    ).fetchone()
    conn.close()
    assert row[0] == "NG_007524.2"


def test_insert_lrg_transcript_data_roundtrip(sqlite_insert):
    """INSERT OR REPLACE into LRG_transcripts."""
    result = sqlite_insert.insert_lrg_transcript_data(["LRG_1t1", "NM_007294.3"])
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    row = conn.execute(
        "SELECT RefSeqTranscriptID FROM LRG_transcripts WHERE LRGtranscriptID = ?",
        ("LRG_1t1",),
    ).fetchone()
    conn.close()
    assert row[0] == "NM_007294.3"


def test_insert_lrg_protein_data_roundtrip(sqlite_insert):
    """INSERT OR REPLACE into LRG_proteins."""
    result = sqlite_insert.insert_lrg_protein_data("LRG_1p1", "NP_009225.1")
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    row = conn.execute(
        "SELECT RefSeqProteinID FROM LRG_proteins WHERE LRGproteinID = ?",
        ("LRG_1p1",),
    ).fetchone()
    conn.close()
    assert row[0] == "NP_009225.1"


# -- stableGeneIds ------------------------------------------------------------

def test_insert_gene_stable_ids_roundtrip(sqlite_insert):
    """INSERT OR REPLACE into stableGeneIds."""
    data = {
        "hgnc_id": "HGNC:1100",
        "hgnc_symbol": "BRCA1",
        "entrez_id": "672",
        "ensembl_gene_id": "ENSG00000012048",
        "omim_id": "113705",
        "ucsc_id": "uc002icp.4",
        "vega_id": "OTTHUMG00000022110",
        "ccds_id": "CCDS14762",
    }
    result = sqlite_insert.insert_gene_stable_ids(data)
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    row = conn.execute(
        "SELECT hgnc_symbol, ensembl_gene_id FROM stableGeneIds WHERE hgnc_id = ?",
        ("HGNC:1100",),
    ).fetchone()
    conn.close()
    assert row[0] == "BRCA1"
    assert row[1] == "ENSG00000012048"


def test_update_gene_stable_ids(sqlite_insert):
    """update_gene_stable_ids() should modify an existing stableGeneIds row."""
    data = {
        "hgnc_id": "HGNC:5555",
        "hgnc_symbol": "ORIG",
        "entrez_id": "5555",
        "ensembl_gene_id": "ENSG00000000001",
        "omim_id": "000001",
        "ucsc_id": "uc001aaa.1",
        "vega_id": "OTTHUMG00000000001",
        "ccds_id": "CCDS00001",
    }
    sqlite_insert.insert_gene_stable_ids(data)

    updated = dict(data)
    updated["hgnc_symbol"] = "UPDATED"
    updated["entrez_id"] = "9999"
    result = sqlite_insert.update_gene_stable_ids(updated)
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    row = conn.execute(
        "SELECT hgnc_symbol, entrez_id FROM stableGeneIds WHERE hgnc_id = ?",
        ("HGNC:5555",),
    ).fetchone()
    conn.close()
    assert row[0] == "UPDATED"
    assert row[1] == "9999"


# -- version ------------------------------------------------------------------

def test_update_db_version(sqlite_insert):
    """update_db_version() should write into the version table."""
    # Seed the version table with one row first
    conn = sqlite_insert.get_conn()
    conn.execute("INSERT OR IGNORE INTO version(current_version) VALUES(?)", ("0.0.0",))
    conn.commit()
    conn.close()

    result = sqlite_insert.update_db_version("1.2.3")
    assert result == 'true'

    conn = sqlite_insert.get_conn()
    row = conn.execute("SELECT current_version FROM version").fetchone()
    conn.close()
    assert row[0] == "1.2.3"


# ── Task 5: SQLiteDatabase ───────────────────────────────────────────────────

def _load_vvdatabase():
    """Import vvDatabase directly, reusing the stubs set up by _load_vvdbinsert."""
    import types

    repo_root = os.path.join(os.path.dirname(__file__), "..")
    modules_dir = os.path.join(repo_root, "VariantValidator", "modules")

    # Ensure parent package stubs exist (idempotent)
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
        short_name = name.split(".")[-1]
        parent = sys.modules.get("VariantValidator.modules")
        if parent is not None:
            setattr(parent, short_name, mod)
        return mod

    # Ensure utils stub exists
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

    # Stub heavy siblings imported by vvDatabase at module level
    for stub_name in ("vvhgvs", "vvhgvs.exceptions"):
        if stub_name not in sys.modules:
            stub = types.ModuleType(stub_name)
            if stub_name == "vvhgvs.exceptions":
                stub.HGVSDataNotAvailableError = Exception
            sys.modules[stub_name] = stub

    # Load the full dependency chain
    _load_file("VariantValidator.modules.vvDBInit", "vvDBInit.py")
    _load_file("VariantValidator.modules.vvDBGet", "vvDBGet.py")
    _load_file("VariantValidator.modules.vvDBInsert", "vvDBInsert.py")
    return _load_file("VariantValidator.modules.vvDatabase", "vvDatabase.py")


@pytest.fixture(scope="session")
def vvdatabase_module():
    return _load_vvdatabase()


def test_sqlite_database_instantiates(sqlite_db_path, vvdatabase_module):
    """SQLiteDatabase must be importable and instantiable."""
    db = vvdatabase_module.SQLiteDatabase(sqlite_db_path)
    assert db.db_path == sqlite_db_path


def test_sqlite_database_inherits_get_methods(sqlite_db_path, vvdatabase_module):
    """SQLiteDatabase inherits GET methods from SQLiteDBGet."""
    db = vvdatabase_module.SQLiteDatabase(sqlite_db_path)
    result = db.get_uta("BRCA1")
    assert result[0] == 'none'  # empty DB, but method exists and returns sentinel


def test_sqlite_database_inherits_insert_methods(sqlite_db_path, vvdatabase_module):
    """SQLiteDatabase inherits INSERT methods from SQLiteDBInsert."""
    db = vvdatabase_module.SQLiteDatabase(sqlite_db_path)
    db.insert_transcript_info(
        refseq_id="NM_000059.4",
        description="BRCA2",
        transcript_variant="NM_000059",
        current_version="4",
        hgnc_symbol="BRCA2",
        uta_symbol="BRCA2",
    )
    conn = db.get_conn()
    rows = conn.execute(
        "SELECT refSeqID FROM transcript_info WHERE hgncSymbol = ?", ("BRCA2",)
    ).fetchall()
    conn.close()
    assert len(rows) == 1
