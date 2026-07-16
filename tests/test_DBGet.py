from unittest.mock import MagicMock
import pytest

from VariantValidator.modules.vvDBGet import Mixin


def make_db():
    db = Mixin.__new__(Mixin)
    return db


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
