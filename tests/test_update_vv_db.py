import pytest
from unittest.mock import patch, MagicMock
import VariantValidator.update_vv_db as uv

# -----------------------------
# Fixtures
# -----------------------------

@pytest.fixture
def fake_db():
    """A fake database object with the methods used by update_vv_db."""
    db = MagicMock()
    db.conn = MagicMock()
    db.conn.commit = MagicMock()
    db.execute = MagicMock()
    db.update_refseqgene_loci = MagicMock()
    db.update_lrg_rs_lookup = MagicMock()
    db.update_lrgt_rst = MagicMock()
    db.update_lrg_p_rs_p_lookup = MagicMock()
    db.get_gene_symbol_from_refseq_id = MagicMock(return_value="none")
    return db

# -----------------------------
# Tests
# -----------------------------

def test_connect_returns_database_object():
    dummy_config = {
        'user': 'u', 'password': 'p', 'host': 'h', 'port': '3306', 'database': 'db'
    }
    with patch('VariantValidator.update_vv_db.ConfigParser') as mock_cp_class, \
         patch('VariantValidator.update_vv_db.vvDatabase.Database') as mock_db_class:
        mock_cp = MagicMock()
        mock_cp.__getitem__.return_value = dummy_config
        mock_cp_class.return_value = mock_cp

        mock_db_instance = MagicMock()
        mock_db_class.return_value = mock_db_instance

        db = uv.connect()
        assert db == mock_db_instance
        mock_db_class.assert_called_once()


def test_delete_executes_all_statements(fake_db):
    uv.delete.__globals__['connect'] = MagicMock(return_value=fake_db)

    uv.delete()
    expected_sql = [
        'DELETE FROM transcript_info',
        'DELETE FROM refSeqGene_loci',
        'DELETE FROM LRG_transcripts',
        'DELETE FROM LRG_proteins',
        'DELETE FROM LRG_RSG_lookup',
    ]
    executed_sql = [call[0][0] for call in fake_db.execute.call_args_list]
    for sql in expected_sql:
        assert sql in executed_sql
    fake_db.conn.commit.assert_called_once()


@patch('VariantValidator.update_vv_db.update_refseq')
@patch('VariantValidator.update_vv_db.update_lrg')
def test_update_calls_refseq_and_lrg(mock_update_lrg, mock_update_refseq):
    fake_db = MagicMock()
    uv.connect = MagicMock(return_value=fake_db)

    uv.update()

    uv.connect.assert_called_once()
    mock_update_refseq.assert_called_once_with(fake_db)
    mock_update_lrg.assert_called_once_with(fake_db)


@patch('VariantValidator.update_vv_db.requests.get')
def test_update_lrg_executes_all_updates(mock_requests, fake_db):
    mock_requests.side_effect = [
        MagicMock(text="#comment\nLRG_1 GENE1 RSG123 tx1 rst1"),
        MagicMock(text="#comment\nLRG_1 something active"),
        MagicMock(text="#comment\nLRG_1P RS_1P"),
    ]
    uv.update_lrg(fake_db)

    assert fake_db.update_lrg_rs_lookup.called
    assert fake_db.update_lrgt_rst.called
    assert fake_db.update_lrg_p_rs_p_lookup.called


def test_map_line_builds_list_correctly():
    line = {'rsg_id': 'RSG123', 'chr_id': 'NC_000001', 'rsg_start': '100', 'rsg_end': '200', 'ori': '+'}
    rsg_info = [{'rsg_id': 'RSG123', 'symbol': 'GENE1', 'gene_id': '1234'}]
    ml = uv.map_line(line, 'GRCh38', rsg_info)

    assert ml[0] == 'RSG123'
    assert ml[1] == 'NC_000001'
    assert ml[2] == 'GRCh38'
    assert ml[6] == 'GENE1'
    assert ml[7] == '1234'

def test_drop_core_indexes_success(fake_db):
    uv.drop_core_indexes(fake_db)

    assert fake_db.execute.call_count == 6


def test_drop_core_indexes_index_missing(fake_db):
    fake_db.execute.side_effect = Exception("missing")

    # Should not raise
    uv.drop_core_indexes(fake_db)

    assert fake_db.execute.call_count == 6


def test_rebuild_core_indexes_success(fake_db):
    uv.rebuild_core_indexes(fake_db)

    assert fake_db.execute.call_count == 6


def test_rebuild_core_indexes_failure(fake_db):
    fake_db.execute.side_effect = Exception("boom")

    # Should not raise
    uv.rebuild_core_indexes(fake_db)

    assert fake_db.execute.call_count == 6


def test_count_ng_nc_non_mapping():
    assert uv.count_ng_nc("hello world") is None


def test_count_ng_nc_comment():
    line = "# NC_000001 NG_000001"
    assert uv.count_ng_nc(line) == "failed"


def test_count_ng_nc_rejected():
    line = (
        "NC_000001\tRefSeq\tmatch\t1\t2\t.\t+\t.\t"
        "ID=x;Target=NG_000001.1 1 2 +"
    )
    assert uv.count_ng_nc(line) == "rejected"


def test_count_ng_nc_bad_columns():
    line = (
        "NC_000001\tRefSeq\tmatch\t1\t2\t.\t+\t."
        "\tID=x"
    )
    assert uv.count_ng_nc(line) == None


def test_count_ng_nc_success():
    line = (
        "NC_000001\tRefSeq\tmatch\t100\t200\t.\t+\t.\t"
        "ID=x;Target=NG_000001.1 1 100 +;gap_count=0"
    )

    result = uv.count_ng_nc(line)

    assert result == {
        "rsg_id": "NG_000001.1",
        "chr_id": "NC_000001",
        "rsg_start": "100",
        "rsg_end": "200",
        "ori": "+"
    }


def test_map_line_no_matching_rsg():
    line = {
        "rsg_id": "NG_1",
        "chr_id": "NC_1",
        "rsg_start": "1",
        "rsg_end": "2",
        "ori": "+"
    }

    result = uv.map_line(line, "GRCh38", [])

    assert result == [
        "NG_1",
        "NC_1",
        "GRCh38",
        "1",
        "2",
        "+"
    ]


@patch("VariantValidator.update_vv_db.connect")
def test_delete_connect_called(mock_connect, fake_db):
    mock_connect.return_value = fake_db

    uv.delete()

    mock_connect.assert_called_once()


@patch("VariantValidator.update_vv_db.connect")
def test_update_connect_called(mock_connect):
    fake_db = MagicMock()
    mock_connect.return_value = fake_db

    with patch("VariantValidator.update_vv_db.update_refseq") as mock_refseq, \
         patch("VariantValidator.update_vv_db.update_lrg") as mock_lrg:

        uv.update()

        mock_connect.assert_called_once()
        mock_refseq.assert_called_once_with(fake_db)
        mock_lrg.assert_called_once_with(fake_db)

def test_drop_core_indexes_calls_all(fake_db):
    uv.drop_core_indexes(fake_db)

    assert fake_db.execute.call_count == 6


def test_rebuild_core_indexes_calls_all(fake_db):
    uv.rebuild_core_indexes(fake_db)

    assert fake_db.execute.call_count == 6


def test_drop_core_indexes_ignores_errors(fake_db):
    fake_db.execute.side_effect = Exception("boom")

    # Should not raise
    uv.drop_core_indexes(fake_db)

    assert fake_db.execute.call_count == 6


def test_rebuild_core_indexes_ignores_errors(fake_db):
    fake_db.execute.side_effect = Exception("boom")

    # Should not raise
    uv.rebuild_core_indexes(fake_db)

    assert fake_db.execute.call_count == 6


def test_count_ng_nc_returns_none():
    assert uv.count_ng_nc("hello world") is None


def test_count_ng_nc_comment_b():
    assert uv.count_ng_nc("# NC_000001 NG_000001") == "failed"


def test_count_ng_nc_bad_columns_b():
    line = (
        "NC_000001\tRefSeq\tmatch\t1\t2\t.\t+\t.\t"
        "ID=x;Target=NG_000001.1 1 2 +;gap_count=0"
        "\textra"
    )

    assert uv.count_ng_nc(line) == "failed"


def test_map_line_no_lookup():
    line = {
        "rsg_id": "NG_1",
        "chr_id": "NC_1",
        "rsg_start": "1",
        "rsg_end": "2",
        "ori": "+",
    }

    result = uv.map_line(
        line,
        "GRCh38",
        [],
    )

    assert result == [
        "NG_1",
        "NC_1",
        "GRCh38",
        "1",
        "2",
        "+",
    ]


@patch("VariantValidator.update_vv_db.requests.get")
def test_update_lrg_skips_comments_and_short_lines(mock_requests, fake_db):
    mock_requests.side_effect = [
        MagicMock(
            text="#comment\nshort\nLRG_1 GENE1 RSG1 tx1 NM_1"
        ),
        MagicMock(
            text="#comment\nLRG_1 x active"
        ),
        MagicMock(
            text="#comment\nLRG_1p NP_1"
        ),
    ]

    uv.update_lrg(fake_db)

    fake_db.update_lrg_rs_lookup.assert_called_once()
    fake_db.update_lrgt_rst.assert_called_once()
    fake_db.update_lrg_p_rs_p_lookup.assert_called_once()

@patch("VariantValidator.update_vv_db.requests.get")
def test_update_refseq_latest_assembly_http_error(mock_get, fake_db):
    """HTTP failure when retrieving the latest assembly listing."""

    rsg = MagicMock()
    rsg.text = ""

    latest = MagicMock()
    latest.raise_for_status.side_effect = RuntimeError("HTTP error")

    mock_get.side_effect = [rsg, latest]

    with pytest.raises(RuntimeError):
        uv.update_refseq(fake_db)


@patch("VariantValidator.update_vv_db.requests.get")
def test_update_refseq_missing_gcf_directory(mock_get, fake_db):
    """No GCF assembly directory is present."""

    rsg = MagicMock()
    rsg.text = ""

    latest = MagicMock()
    latest.raise_for_status.return_value = None
    latest.text = "<html></html>"

    mock_get.side_effect = [rsg, latest]

    with pytest.raises(IndexError):
        uv.update_refseq(fake_db)


@patch("VariantValidator.update_vv_db.requests.get")
def test_update_refseq_missing_genomic_gff(mock_get, fake_db):
    """Assembly directory exists but contains no genomic gff."""

    rsg = MagicMock()
    rsg.text = ""

    latest = MagicMock()
    latest.raise_for_status.return_value = None
    latest.text = '<a href="GCF_test/">'

    assembly = MagicMock()
    assembly.raise_for_status.return_value = None
    assembly.text = "<html></html>"

    mock_get.side_effect = [
        rsg,
        latest,
        assembly,
    ]

    with pytest.raises(IndexError):
        uv.update_refseq(fake_db)

@patch("VariantValidator.update_vv_db.gzip.decompress")
@patch("VariantValidator.update_vv_db.requests.get")
def test_update_refseq_single_valid_mapping(mock_get, mock_decompress, fake_db):
    """Run update_refseq() with one genuine NG->NC mapping."""

    fake_db.get_gene_symbol_from_refseq_id.return_value = "none"

    mock_get.side_effect = [
        # gene_RefSeqGene
        MagicMock(
            text="0\t1234\tGENE1\tNG_000001.1"
        ),

        # latest_assembly_versions
        MagicMock(
            text='href="GCF_000001405.40_GRCh38.p14/"',
            raise_for_status=MagicMock(),
        ),

        # assembly directory
        MagicMock(
            text='href="GCF_000001405.40_GRCh38.p14_genomic.gff.gz"',
            raise_for_status=MagicMock(),
        ),

        # downloaded file
        MagicMock(
            content=b"dummy",
            raise_for_status=MagicMock(),
        ),
    ]

    mock_decompress.return_value = (
        b"NC_000001.11\tRefSeq\tmatch\t100\t200\t.\t+\t."
        b"\tID=x;Target=NG_000001.1 1 100 +;gap_count=0"
    )

    uv.update_refseq(fake_db)

    fake_db.update_refseqgene_loci.assert_called_once()

    written = fake_db.update_refseqgene_loci.call_args.args[0]

    assert written[0] == "NG_000001.1"
    assert written[1] == "NC_000001.11"
    assert written[2] == "GRCh38"
    assert written[9] == "1234"
    assert written[10] == "GENE1"


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
