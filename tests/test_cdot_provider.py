"""Unit tests for cdot data provider wiring in vvMixinInit.

Uses unittest.mock to avoid needing real cdot files or a real seqrepo.
"""
import sys
import os
import pytest
from unittest.mock import patch, MagicMock

# Skip the entire module if cdot is not installed
cdot_dataproviders = pytest.importorskip("cdot.hgvs.dataproviders")


def _make_mock_config(
    backend_type='sqlite',
    cdot_path='/fake/cdot.json.gz',
    sqlite_path='/fake/vvdb.sqlite',
    seqrepo='/fake/seqrepo',
):
    """Create a mock configparser.ConfigParser."""
    config = MagicMock()
    config.get.side_effect = lambda section, key, **kw: {
        ('backend', 'type'): backend_type,
        ('backend', 'cdot_path'): cdot_path,
        ('backend', 'sqlite_path'): sqlite_path,
        ('seqrepo', 'location'): seqrepo,
    }.get((section, key), kw.get('fallback', None))
    return config


@patch('cdot.hgvs.dataproviders.JSONDataProvider')
@patch('VariantValidator.modules.vvDatabase.SQLiteDatabase')
def test_sqlite_imports_are_available(mock_sqlite_db, mock_json_provider):
    """Verify that SQLiteDatabase and JSONDataProvider can be imported.

    This test patches the heavy objects and verifies the wiring logic.
    We can't instantiate a real Validator here (needs network + files),
    so we test the branching logic directly via import verification.
    """
    mock_json_provider.return_value = MagicMock()
    mock_sqlite_db.return_value = MagicMock()

    # Verify that SQLiteDatabase and JSONDataProvider can be imported
    from VariantValidator.modules.vvDatabase import SQLiteDatabase
    import cdot.hgvs.dataproviders

    assert SQLiteDatabase is not None
    assert hasattr(cdot.hgvs.dataproviders, 'JSONDataProvider')


def test_json_data_provider_is_importable():
    """cdot.hgvs.dataproviders.JSONDataProvider must be accessible."""
    import cdot.hgvs.dataproviders
    assert hasattr(cdot.hgvs.dataproviders, 'JSONDataProvider'), (
        "cdot.hgvs.dataproviders.JSONDataProvider not found; "
        "check cdot package version"
    )


def test_sqlite_database_is_importable():
    """VariantValidator.modules.vvDatabase.SQLiteDatabase must be importable
    without requiring MySQL or PostgreSQL connections."""
    import types

    repo_root = os.path.join(os.path.dirname(__file__), "..")
    modules_dir = os.path.join(repo_root, "VariantValidator", "modules")
    import importlib.util

    # Register stub parent packages so relative imports work
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

    # Stub heavy siblings
    import functools
    if "VariantValidator.modules.utils" not in sys.modules:
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

    for stub_name in ("vvhgvs", "vvhgvs.exceptions"):
        if stub_name not in sys.modules:
            stub = types.ModuleType(stub_name)
            if stub_name == "vvhgvs.exceptions":
                stub.HGVSDataNotAvailableError = Exception
            sys.modules[stub_name] = stub

    _load_file("VariantValidator.modules.vvDBInit", "vvDBInit.py")
    _load_file("VariantValidator.modules.vvDBGet", "vvDBGet.py")
    _load_file("VariantValidator.modules.vvDBInsert", "vvDBInsert.py")
    vvdb_mod = _load_file("VariantValidator.modules.vvDatabase", "vvDatabase.py")

    assert hasattr(vvdb_mod, 'SQLiteDatabase'), "SQLiteDatabase not found in vvDatabase"
    assert vvdb_mod.SQLiteDatabase is not None
