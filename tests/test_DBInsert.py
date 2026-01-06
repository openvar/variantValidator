# tests/test_DBInsert.py
import pytest
from unittest.mock import MagicMock, patch
from VariantValidator.modules import vvDBInsert


# ------------------- Fixtures -------------------

def dummy_init(self, db_config):
    """Patch vvDBGet.Mixin.__init__ to avoid creating real DB pool."""
    self.pool = None
    self.db_config = db_config


@pytest.fixture
def mixin_instance():
    with patch('VariantValidator.modules.vvDBGet.Mixin.__init__', dummy_init):
        return vvDBInsert.Mixin(db_config={})


@pytest.fixture
def mock_db_objects():
    mock_cursor = MagicMock()
    mock_conn = MagicMock()
    return mock_conn, mock_cursor


def patch_db_methods(mixin_instance, mock_conn, mock_cursor):
    return patch.object(mixin_instance, 'get_conn', return_value=mock_conn), \
           patch.object(mixin_instance, 'get_cursor', return_value=mock_cursor)


# ------------------- Tests -------------------

@pytest.mark.parametrize("lastrowid", [1, None])
def test_insert_methods(mixin_instance, mock_db_objects, lastrowid):
    mock_conn, mock_cursor = mock_db_objects
    mock_cursor.lastrowid = lastrowid

    with patch_db_methods(mixin_instance, mock_conn, mock_cursor)[0], \
         patch_db_methods(mixin_instance, mock_conn, mock_cursor)[1]:

        # insert transcript_info
        entry = 'NM_000000.1'
        data = ['NM_000000.1', 'desc', 'variant', 'v1', 'HGNC1', 'UTA1']
        expected = 'true' if lastrowid else 'Unknown error'
        assert mixin_instance.insert(entry, data, 'transcript_info') == expected

        # insert_refseq_gene_data
        rsg_data = [1, 'chr1', 'GRCh38', 100, 200, '+', 100, 'chr1:100-200', 'rsg:100-200', 101, 'HGNC1']
        expected = 'true' if lastrowid else 'Unknown error'
        assert mixin_instance.insert_refseq_gene_data(rsg_data) == expected

        # insert_refseq_gene_id_from_lrg_id
        lrg_rs_lookup = ['LRG_1', 'HGNC1', 'RSG_1', 'active']
        expected = 'true' if lastrowid else 'Unknown error'
        assert mixin_instance.insert_refseq_gene_id_from_lrg_id(lrg_rs_lookup) == expected

        # insert_lrg_transcript_data
        lrgtx_to_rst_id = ['LRGTx1', 'RST1']
        expected = 'true' if lastrowid else 'Unknown error'
        assert mixin_instance.insert_lrg_transcript_data(lrgtx_to_rst_id) == expected

        # insert_lrg_protein_data
        expected = 'true' if lastrowid else 'Unknown error'
        assert mixin_instance.insert_lrg_protein_data('LRGP1', 'RSP1') == expected

        # insert_gene_stable_ids
        gene_data = {
            'hgnc_id': 'HGNC:1',
            'hgnc_symbol': 'HGNC1',
            'entrez_id': 101,
            'ensembl_gene_id': 'ENSG000001',
            'omim_id': 'OMIM1',
            'ucsc_id': 'UC001',
            'vega_id': 'VEGA1',
            'ccds_id': 'CCDS1'
        }
        # This method returns lowercase 'unknown error'
        expected = 'true' if lastrowid else 'unknown error'
        assert mixin_instance.insert_gene_stable_ids(gene_data) == expected


@pytest.mark.parametrize("lastrowid", [1])
def test_update_methods(mixin_instance, mock_db_objects, lastrowid):
    mock_conn, mock_cursor = mock_db_objects
    mock_cursor.lastrowid = lastrowid

    with patch_db_methods(mixin_instance, mock_conn, mock_cursor)[0], \
         patch_db_methods(mixin_instance, mock_conn, mock_cursor)[1]:

        # update transcript_info
        entry = 'NM_000000.1'
        data = ['NM_000000.1', 'desc', 'variant', 'v2', 'HGNC1', 'UTA1']
        assert mixin_instance.update(entry, data) == 'true'

        # update_refseq_gene_data
        rsg_data = [1, 'chr1', 'GRCh38', 100, 200, '+', 100, 'chr1:100-200', 'rsg:100-200', 101, 'HGNC2']
        assert mixin_instance.update_refseq_gene_data(rsg_data) == 'true'

        # update_gene_stable_ids
        gene_stable_ids = {
            'hgnc_id': 'HGNC:1',
            'hgnc_symbol': 'HGNC2',
            'entrez_id': 102,
            'ensembl_gene_id': 'ENSG000002',
            'omim_id': 'OMIM2',
            'ucsc_id': 'UC002',
            'vega_id': 'VEGA2',
            'ccds_id': 'CCDS2'
        }
        assert mixin_instance.update_gene_stable_ids(gene_stable_ids) == 'true'

        # update_db_version
        assert mixin_instance.update_db_version('v2.0') == 'true'

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
