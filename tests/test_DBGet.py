from unittest.mock import MagicMock

import pytest

from VariantValidator import settings
from VariantValidator.modules.vvDBGet import (
    DB_GET_CACHE,
    Mixin,
    _CACHE_MISS,
    _get_cached,
    _set_cached,
    clear_get_cache,
)


def make_db():
    db = Mixin.__new__(Mixin)
    return db


@pytest.fixture(autouse=True)
def reset_db_get_cache(monkeypatch):
    clear_get_cache()
    monkeypatch.setattr(settings, "vvDB_GET_CACHE", True)
    monkeypatch.setattr(settings, "vvDB_GET_CACHE_SIZE", 10000)

    yield

    clear_get_cache()


def test_execute_fetchone_success():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    cursor.fetchone.return_value = ("A", "B")

    result = db.execute(
        "SELECT 1",
        (),
    )

    cursor.execute.assert_called_once_with(
        "SELECT 1",
        (),
    )

    cursor.close.assert_called_once()
    conn.close.assert_called_once()

    assert result == ("A", "B")


def test_execute_fetchone_none():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    cursor.fetchone.return_value = None

    result = db.execute(
        "SELECT 1",
        (),
    )

    assert result == ["none", "No data"]


def test_execute_all_success():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    cursor.fetchall.return_value = [
        ("A",),
        ("B",),
    ]

    result = db.execute_all(
        "SELECT 1",
        (),
    )

    assert result == [
        ("A",),
        ("B",),
    ]


def test_execute_all_empty():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    cursor.fetchall.return_value = []

    result = db.execute_all(
        "SELECT 1",
        (),
    )

    assert result == [["none", "No data"]]


def test_execute_retry_then_success():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    cursor.execute.side_effect = [
        Exception("boom"),
        None,
    ]
    cursor.fetchone.return_value = ("OK",)

    result = db.execute(
        "SELECT 1",
        (),
    )

    assert cursor.execute.call_count == 2
    conn.reconnect.assert_called_once_with(
        attempts=1,
        delay=0,
    )
    assert result == ("OK",)


def test_execute_retry_then_fail():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    cursor.execute.side_effect = Exception("boom")

    with pytest.raises(Exception):
        db.execute(
            "SELECT 1",
            (),
        )

    assert cursor.execute.call_count == 3
    assert conn.reconnect.call_count == 2


def test_execute_all_retry_then_success():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    cursor.execute.side_effect = [
        Exception("boom"),
        None,
    ]
    cursor.fetchall.return_value = [("A",)]

    result = db.execute_all(
        "SELECT 1",
        (),
    )

    assert cursor.execute.call_count == 2
    conn.reconnect.assert_called_once_with(
        attempts=1,
        delay=0,
    )
    assert result == [("A",)]


def test_execute_all_retry_then_fail():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    cursor.execute.side_effect = Exception("boom")

    with pytest.raises(Exception):
        db.execute_all(
            "SELECT 1",
            (),
        )

    assert cursor.execute.call_count == 3
    assert conn.reconnect.call_count == 2


def test_get_uta():
    db = make_db()

    db.execute = MagicMock(return_value=["UTA"])

    assert db.get_uta("GENE1") == ["UTA"]

    db.execute.assert_called_once_with(
        "SELECT utaSymbol FROM transcript_info WHERE hgncSymbol = %s",
        ("GENE1",),
    )


def test_get_hgnc():
    db = make_db()

    db.execute = MagicMock(return_value=["HGNC"])

    assert db.get_hgnc("UTA1") == ["HGNC"]

    db.execute.assert_called_once_with(
        "SELECT hgncSymbol FROM transcript_info WHERE utaSymbol = %s",
        ("UTA1",),
    )


def test_get_db_version():
    db = make_db()

    db.execute = MagicMock(return_value=["1.2"])

    assert db.get_db_version() == ["1.2"]

    db.execute.assert_called_once_with(
        "SELECT current_version FROM version",
    )


def test_get_transcript_description():
    db = make_db()

    db.execute = MagicMock(return_value=["description"])

    assert (
        db.get_transcript_description("NM_000001.1")
        == "description"
    )

    db.execute.assert_called_once_with(
        "SELECT description FROM transcript_info WHERE refSeqID = %s",
        ("NM_000001.1",),
    )


def test_get_transcript_annotation():
    db = make_db()

    db.execute = MagicMock(return_value=["annotation"])

    assert (
        db.get_transcript_annotation("NM_000001.1")
        == "annotation"
    )

    db.execute.assert_called_once_with(
        "SELECT transcriptVariant FROM transcript_info WHERE refSeqID = %s",
        ("NM_000001.1",),
    )


def test_get_gene_symbol_from_transcript_id():
    db = make_db()

    db.execute = MagicMock(return_value=["GENE1"])

    assert (
        db.get_gene_symbol_from_transcript_id("NM_000001.1")
        == "GENE1"
    )

    db.execute.assert_called_once_with(
        "SELECT hgncSymbol FROM transcript_info WHERE refSeqID = %s",
        ("NM_000001.1",),
    )


def test_get_refseq_data_by_refseq_id():
    db = make_db()

    db.execute = MagicMock(return_value=["NG_000001"])

    assert (
        db.get_refseq_data_by_refseq_id(
            "NG_000001.1",
            "GRCh38",
        )
        == ["NG_000001"]
    )

    db.execute.assert_called_once_with(
        "SELECT refSeqGeneID, refSeqChromosomeID, genomeBuild, startPos, endPos, "
        "orientation, totalLength, chrPos, rsgPos, entrezID, hgncSymbol "
        "FROM refSeqGene_loci WHERE refSeqGeneID = %s AND genomeBuild = %s",
        ("NG_000001.1", "GRCh38"),
    )


def test_cache_set_and_get():
    key = ("test", "A")

    assert _get_cached(key) is _CACHE_MISS

    _set_cached(key, "value")

    assert _get_cached(key) == "value"


def test_cache_disabled(monkeypatch):
    monkeypatch.setattr(settings, "vvDB_GET_CACHE", False)

    key = ("test", "A")

    _set_cached(key, "value")

    assert key not in DB_GET_CACHE
    assert _get_cached(key) is _CACHE_MISS


def test_cache_clear():
    _set_cached(("test", "A"), "A")
    _set_cached(("test", "B"), "B")

    assert len(DB_GET_CACHE) == 2

    clear_get_cache()

    assert len(DB_GET_CACHE) == 0


def test_cache_respects_maximum_size(monkeypatch):
    monkeypatch.setattr(settings, "vvDB_GET_CACHE_SIZE", 2)

    _set_cached(("test", "A"), "A")
    _set_cached(("test", "B"), "B")
    _set_cached(("test", "C"), "C")

    assert len(DB_GET_CACHE) == 2
    assert ("test", "A") not in DB_GET_CACHE
    assert ("test", "B") in DB_GET_CACHE
    assert ("test", "C") in DB_GET_CACHE


def test_cache_hit_refreshes_lru_order(monkeypatch):
    monkeypatch.setattr(settings, "vvDB_GET_CACHE_SIZE", 2)

    key_a = ("test", "A")
    key_b = ("test", "B")
    key_c = ("test", "C")

    _set_cached(key_a, "A")
    _set_cached(key_b, "B")

    assert _get_cached(key_a) == "A"

    _set_cached(key_c, "C")

    assert key_a in DB_GET_CACHE
    assert key_b not in DB_GET_CACHE
    assert key_c in DB_GET_CACHE


@pytest.mark.parametrize(
    ("method_name", "argument", "db_result", "expected"),
    [
        (
            "get_transcript_description",
            "NM_CACHE.1",
            ["description"],
            "description",
        ),
        (
            "get_transcript_annotation",
            "NM_CACHE.1",
            ["annotation"],
            "annotation",
        ),
        (
            "get_gene_symbol_from_transcript_id",
            "NM_CACHE.1",
            ["GENE1"],
            "GENE1",
        ),
        (
            "get_lrg_id_from_refseq_gene_id",
            "NG_CACHE.1",
            ["LRG_1", "public"],
            ["LRG_1", "public"],
        ),
        (
            "get_lrg_protein_id_from_ref_seq_protein_id",
            "NP_CACHE.1",
            ["LRG_1p1"],
            "LRG_1p1",
        ),
        (
            "get_lrg_data_from_lrg_id",
            "LRG_1",
            ["data"],
            ["data"],
        ),
        (
            "get_stable_gene_id_info",
            "GENE1",
            ["stable"],
            ["stable"],
        ),
    ],
)
def test_cached_getters_only_query_database_once(
        method_name,
        argument,
        db_result,
        expected,
):
    db = make_db()
    db.execute = MagicMock(return_value=db_result)

    method = getattr(db, method_name)

    assert method(argument) == expected
    assert method(argument) == expected

    db.execute.assert_called_once()


def test_cached_getter_queries_again_when_cache_disabled(
        monkeypatch,
):
    monkeypatch.setattr(settings, "vvDB_GET_CACHE", False)

    db = make_db()
    db.execute = MagicMock(return_value=["annotation"])

    assert (
        db.get_transcript_annotation("NM_CACHE.1")
        == "annotation"
    )
    assert (
        db.get_transcript_annotation("NM_CACHE.1")
        == "annotation"
    )

    assert db.execute.call_count == 2

def test_execute_write_success():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    db.execute_write(
        "ALTER TABLE transcript_info DROP INDEX refSeqID_index"
    )

    cursor.execute.assert_called_once_with(
        "ALTER TABLE transcript_info DROP INDEX refSeqID_index"
    )
    conn.commit.assert_called_once()
    cursor.close.assert_called_once()
    conn.close.assert_called_once()


def test_execute_write_with_parameters():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    db.execute_write(
        "UPDATE transcript_info SET hgncSymbol = %s WHERE refSeqID = %s",
        ("GENE1", "NM_000001.1"),
    )

    cursor.execute.assert_called_once_with(
        "UPDATE transcript_info SET hgncSymbol = %s WHERE refSeqID = %s",
        ("GENE1", "NM_000001.1"),
    )
    conn.commit.assert_called_once()
    cursor.close.assert_called_once()
    conn.close.assert_called_once()


def test_execute_write_failure_closes_connection():
    db = make_db()

    conn = MagicMock()
    cursor = MagicMock()

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    cursor.execute.side_effect = Exception("boom")

    with pytest.raises(Exception, match="boom"):
        db.execute_write(
            "ALTER TABLE transcript_info DROP INDEX refSeqID_index"
        )

    conn.commit.assert_not_called()
    cursor.close.assert_called_once()
    conn.close.assert_called_once()


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
