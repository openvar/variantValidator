import unittest
from unittest.mock import Mock, patch
from VariantValidator.modules.hgvs_utils import pvcf_to_hgvs, PseudoVCF2HGVSError

class TestPVCFtoHGVS(unittest.TestCase):

    def setUp(self):
        # Mock validator object with .hn.normalize
        self.mock_validator = Mock()
        self.mock_validator.hn.normalize = Mock(side_effect=lambda x: x)
        self.mock_validator.db.get_refseq_id_from_lrg_id = Mock(side_effect=lambda x: "NM_000000.1")

        # Mock reverse_normalizer
        self.mock_reverse = Mock()
        self.mock_reverse.normalize = Mock(side_effect=lambda x: x)

    @patch('VariantValidator.modules.hgvs_utils.seq_data.to_accession')
    @patch('VariantValidator.modules.hgvs_utils.hgvs_delins_parts_to_hgvs_obj')
    def test_simple_substitution(self, mock_hgvs_obj, mock_to_accession):
        # Setup mocks
        mock_to_accession.return_value = "NM_000000.1"
        mock_hgvs_obj.return_value = "HGVS_OBJ"

        # Provide a simple pVCF string
        query = "chr1-123-A-T"

        # Call the function
        result = pvcf_to_hgvs(query, selected_assembly="GRCh38", normalization_direction=3,
                              reverse_normalizer=self.mock_reverse, validator=self.mock_validator)

        # Assertions
        self.assertEqual(result, "HGVS_OBJ")
        self.mock_validator.hn.normalize.assert_called_once()
        mock_to_accession.assert_called_once_with("1", "GRCh38")
        mock_hgvs_obj.assert_called()  # At least called once

    def test_unsupported_format(self):
        query = "chr1-123-A"  # missing ALT
        with self.assertRaises(PseudoVCF2HGVSError):
            pvcf_to_hgvs(query, selected_assembly="GRCh38", normalization_direction=3,
                          reverse_normalizer=self.mock_reverse, validator=self.mock_validator)

if __name__ == '__main__':
    unittest.main()

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
