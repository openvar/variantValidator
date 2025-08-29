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

from . import configure
from . import logger
from .version import __version__

from .validator import Validator

__all__ = ["Validator"]

# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
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
