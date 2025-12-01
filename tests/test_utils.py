from unittest import TestCase
from VariantValidator.modules import utils
import vvhgvs.parser
import json
from unittest.mock import patch, MagicMock
import requests


class TestHGNCRest(TestCase):
    """Test the hgnc_rest function"""

    def test_non_json(self):
        path = ''
        with self.assertRaises(json.decoder.JSONDecodeError):
            utils.hgnc_rest(path)

    def test_symbol(self):
        path = '/fetch/symbol/NANOG'
        output = utils.hgnc_rest(path)
        self.assertIsInstance(output, dict)
        self.assertEqual(list(output.keys()), ['record', 'error'])
        self.assertNotEqual(output['record'], '')
        self.assertEqual(output['error'], 'false')
        self.assertIsInstance(output['record'], dict)
        self.assertListEqual(list(output['record'].keys()), ['responseHeader', 'response'])
        self.assertGreater(output['record']['response']['numFound'], 0)

    def test_symbol_wrong(self):
        path = '/fetch/symbol/IAMNOTAGENE'
        output = utils.hgnc_rest(path)
        self.assertIsInstance(output, dict)
        self.assertEqual(list(output.keys()), ['record', 'error'])
        self.assertNotEqual(output['record'], '')
        self.assertEqual(output['error'], 'false')
        self.assertIsInstance(output['record'], dict)
        self.assertListEqual(list(output['record'].keys()), ['responseHeader', 'response'])
        self.assertEqual(output['record']['response']['numFound'], 0)


class TestValStr(TestCase):
    """Test the valstr function"""

    def setUp(self):
        self.hp = vvhgvs.parser.Parser()

    def test_string(self):
        var = ''
        with self.assertRaises(AttributeError):
            utils.valstr(var)

    def test_variant_sub(self):
        """ Will test that reference isn't removed """
        stringvar = 'NM_015120.4:c.34C>T'
        var = self.hp.parse_hgvs_variant(stringvar)
        output = utils.valstr(var)
        self.assertEqual(var.posedit.edit.type, 'sub')
        self.assertEqual(output, stringvar)

    def test_variant_identity(self):
        """ Will test that the reference is removed """
        stringvar = 'NM_015120.4:c.34CG='
        var = self.hp.parse_hgvs_variant(stringvar)
        output = utils.valstr(var)
        self.assertEqual(var.posedit.edit.type, 'identity')
        self.assertEqual(output, 'NM_015120.4:c.34=')

    def test_variant_identity2(self):
        """ Will test that the reference is not removed """
        stringvar = 'NM_015120.4:c.34='
        var = self.hp.parse_hgvs_variant(stringvar)
        output = utils.valstr(var)
        self.assertEqual(var.posedit.edit.type, 'identity')
        self.assertEqual(output, 'NM_015120.4:c.34=')


class TestProteinInv(TestCase):
    """Test the pro_inv_info function"""

    def test_empty(self):
        pro1 = ''
        pro2 = ''
        output = utils.pro_inv_info(pro1, pro2)
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'identity')
        self.assertEqual(output['terminate'], 'false')
        self.assertEqual(output['ter_pos'], 0)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '')
        self.assertEqual(output['prot_ins_seq'], '')
        self.assertEqual(output['edit_start'], 0)
        self.assertEqual(output['edit_end'], 0)

    def test_equal(self):
        output = utils.pro_inv_info('MTACGP', 'MTACGP')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'identity')
        self.assertEqual(output['terminate'], 'false')
        self.assertEqual(output['ter_pos'], 0)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '')
        self.assertEqual(output['prot_ins_seq'], '')
        self.assertEqual(output['edit_start'], 0)
        self.assertEqual(output['edit_end'], 0)

    def test_equal_with_ter(self):
        output = utils.pro_inv_info('MTACGP*', 'MTACGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'identity')
        self.assertEqual(output['terminate'], 'false')
        self.assertEqual(output['ter_pos'], 0)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '')
        self.assertEqual(output['prot_ins_seq'], '')
        self.assertEqual(output['edit_start'], 0)
        self.assertEqual(output['edit_end'], 0)

    def test_unequal(self):
        output = utils.pro_inv_info('MTACGP', 'MGCATP')
        self.assertIsNone(output)

    def test_ref_has_ter(self):
        output = utils.pro_inv_info('MTACGP*', 'MTACGPAL')
        self.assertIsNone(output)

    def test_has_ter(self):
        output = utils.pro_inv_info('MTACGP', 'MTACGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 7)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '*')
        self.assertEqual(output['prot_ins_seq'], '**')
        self.assertEqual(output['edit_start'], 7)
        self.assertEqual(output['edit_end'], 7)

    def test_has_ter_inv(self):
        output = utils.pro_inv_info('MTATGLCGP*', 'MTALGTCGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 10)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], 'TGL')
        self.assertEqual(output['prot_ins_seq'], 'LGT')
        self.assertEqual(output['edit_start'], 4)
        self.assertEqual(output['edit_end'], 6)

    def test_has_ter_sub(self):
        output = utils.pro_inv_info('MTATCGP*', 'MTACCGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 8)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], 'T')
        self.assertEqual(output['prot_ins_seq'], 'C')
        self.assertEqual(output['edit_start'], 4)
        self.assertEqual(output['edit_end'], 4)

    def test_has_ter_del(self):
        output = utils.pro_inv_info('MTATCGP*', 'MTACGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 7)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], 'TCGP')
        self.assertEqual(output['prot_ins_seq'], 'CGP*')
        self.assertEqual(output['edit_start'], 4)
        self.assertEqual(output['edit_end'], 7)

    def test_has_ter_ins(self):
        output = utils.pro_inv_info('MTACGP*', 'MTATCGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 8)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '*')
        self.assertEqual(output['prot_ins_seq'], 'T*')
        self.assertEqual(output['edit_start'], 4)
        self.assertEqual(output['edit_end'], 4)


class TestProteinDelIns(TestCase):
    """Test the pro_delins_info function"""

    def test_empty(self):
        pro1 = ''
        pro2 = ''
        output = utils.pro_delins_info(pro1, pro2)
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'identity')
        self.assertEqual(output['terminate'], 'false')
        self.assertEqual(output['ter_pos'], 0)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '')
        self.assertEqual(output['prot_ins_seq'], '')
        self.assertEqual(output['edit_start'], 0)
        self.assertEqual(output['edit_end'], 0)

    def test_equal(self):
        output = utils.pro_delins_info('MTACGP', 'MTACGP')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'identity')
        self.assertEqual(output['terminate'], 'false')
        self.assertEqual(output['ter_pos'], 0)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '')
        self.assertEqual(output['prot_ins_seq'], '')
        self.assertEqual(output['edit_start'], 0)
        self.assertEqual(output['edit_end'], 0)

    def test_equal_with_ter(self):
        output = utils.pro_delins_info('MTACGP*', 'MTACGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'identity')
        self.assertEqual(output['terminate'], 'false')
        self.assertEqual(output['ter_pos'], 0)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '')
        self.assertEqual(output['prot_ins_seq'], '')
        self.assertEqual(output['edit_start'], 0)
        self.assertEqual(output['edit_end'], 0)

    def test_unequal(self):
        output = utils.pro_delins_info('MTACGP', 'MGCATP')
        self.assertIsNone(output)

    def test_ref_has_ter(self):
        output = utils.pro_delins_info('MTACGP*', 'MTACGPAL')
        self.assertIsNone(output, dict)

    def test_has_ter(self):
        output = utils.pro_delins_info('MTACGP', 'MTACGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 7)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '')
        self.assertEqual(output['prot_ins_seq'], '*')
        self.assertEqual(output['edit_start'], 7)
        self.assertEqual(output['edit_end'], 6)

    def test_has_ter_inv(self):
        output = utils.pro_delins_info('MTATGLCGP*', 'MTALGTCGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 10)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], 'TGL')
        self.assertEqual(output['prot_ins_seq'], 'LGT')
        self.assertEqual(output['edit_start'], 4)
        self.assertEqual(output['edit_end'], 6)

    def test_has_ter_sub(self):
        output = utils.pro_delins_info('MTATCGP*', 'MTACCGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 8)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], 'T')
        self.assertEqual(output['prot_ins_seq'], 'C')
        self.assertEqual(output['edit_start'], 4)
        self.assertEqual(output['edit_end'], 4)

    def test_has_ter_del(self):
        output = utils.pro_delins_info('MTATCGP*', 'MTACGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 7)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], 'TCGP')
        self.assertEqual(output['prot_ins_seq'], 'CGP*')
        self.assertEqual(output['edit_start'], 4)
        self.assertEqual(output['edit_end'], 7)

    def test_has_ter_ins(self):
        output = utils.pro_delins_info('MTACGP*', 'MTATCGP*')
        self.assertIsInstance(output, dict)
        print(output)
        self.assertEqual(output['variant'], 'true')
        self.assertEqual(output['terminate'], 'true')
        self.assertEqual(output['ter_pos'], 8)
        self.assertEqual(output['error'], 'false')
        self.assertEqual(output['prot_del_seq'], '')
        self.assertEqual(output['prot_ins_seq'], 'T')
        self.assertEqual(output['edit_start'], 4)
        self.assertEqual(output['edit_end'], 3)


class TestProteinSwap(TestCase):
    """ Test one_to_three function """

    def test_empty(self):
        output = utils.one_to_three('')
        self.assertEqual(output, '')

    def test_wrong(self):
        with self.assertRaises(TypeError):
            utils.one_to_three('Z')

    def test_single(self):
        output = utils.one_to_three('A')
        self.assertEqual(output, 'Ala')

    def test_wrong_pair(self):
        with self.assertRaises(TypeError):
            utils.one_to_three('AZ')

    def test_all(self):
        output = utils.one_to_three('RKDEQNHSTYCWMAILFVPG*')
        self.assertEqual(output, 'ArgLysAspGluGlnAsnHisSerThrTyrCysTrpMetAlaIleLeuPheValProGlyTer')


class TestNInversion(TestCase):
    """ Test n_inversion function. To be honest this looks more like an del+ins """

    def test_empty(self):
        output = utils.n_inversion('', '', '', 0, 0)
        self.assertEqual(output, '')

    def test_empty2(self):
        """
        Warning this output might need checking.
        Passing in 0 as first integer becomes -1 which has meaning!
        """
        output = utils.n_inversion('ATGGAC', '', '', 0, 0)
        self.assertEqual(output, 'ATGGAATGGAC')

    def test_empty3(self):
        output = utils.n_inversion('ATGGAC', '', '', 1, 0)
        self.assertEqual(output, 'ATGGAC')

    def test_correct(self):
        output = utils.n_inversion('ATGGAC', 'GG', 'AA', 3, 4)
        self.assertEqual(output, 'ATAAAC')

    def test_del_incorrect(self):
        output = utils.n_inversion('ATGGAC', 'GC', 'AA', 3, 4)
        self.assertEqual(output, 'error')

    def test_start_incorrect(self):
        output = utils.n_inversion('ATGGAC', 'GG', 'AA', 2, 4)
        self.assertEqual(output, 'error')

    def test_end_incorrect(self):
        output = utils.n_inversion('ATGGAC', 'GG', 'AA', 3, 3)
        self.assertEqual(output, 'error')

    def test_types(self):
        with self.assertRaises(TypeError):
            utils.n_inversion('ATGGAC', 'GG', 'AA', '3', 3)

    def test_types2(self):
        with self.assertRaises(TypeError):
            utils.n_inversion('ATGGAC', '', 0, 0, 0)

    def test_types3(self):
        with self.assertRaises(TypeError):
            utils.n_inversion(0, 0, 0, 0, 0)


class TestHGVSdup2indel(TestCase):
    """ Will test the hgvs_dup2indel function"""

    def setUp(self):
        self.hp = vvhgvs.parser.Parser()

    def test_empty(self):
        with self.assertRaises(AttributeError):
            utils.hgvs_dup2indel('')

    def test_sub(self):
        stringseq = 'NM_015120.4:c.34C>T'
        hgvsseq = self.hp.parse_hgvs_variant(stringseq)
        output = utils.hgvs_dup2indel(hgvsseq)
        self.assertIsInstance(output, str)
        self.assertEqual(output, 'NM_015120.4:c.34_34delCinsCC')

    def test_del(self):
        stringseq = 'NM_015120.4:c.34del'
        hgvsseq = self.hp.parse_hgvs_variant(stringseq)
        output = utils.hgvs_dup2indel(hgvsseq)
        self.assertIsInstance(output, str)
        self.assertEqual(output, 'NM_015120.4:c.34_34delins')

    def test_dup(self):
        stringseq = 'NM_015120.4:c.34dupG'
        hgvsseq = self.hp.parse_hgvs_variant(stringseq)
        output = utils.hgvs_dup2indel(hgvsseq)
        self.assertIsInstance(output, str)
        self.assertEqual(output, 'NM_015120.4:c.34_34delGinsGG')

    def test_dup_pair(self):
        stringseq = 'NM_015120.4:c.34dupGA'
        hgvsseq = self.hp.parse_hgvs_variant(stringseq)
        output = utils.hgvs_dup2indel(hgvsseq)
        self.assertIsInstance(output, str)
        self.assertEqual(output, 'NM_015120.4:c.34_34delGAinsGAGA')

class TestEnsemblRest(TestCase):

    @patch('VariantValidator.modules.utils.requests.get')
    def test_successful_request_grch38(self, mock_get):
        # Mock a successful JSON response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {'id': 'ENST00000380152', 'gene': 'BRCA2'}
        mock_get.return_value = mock_response

        result = utils.ensembl_rest('ENST00000380152.4', '/lookup/id/', 'GRCh38')
        self.assertEqual(result['record'], {'id': 'ENST00000380152', 'gene': 'BRCA2'})
        self.assertEqual(result['error'], 'false')

        # Ensure the request URL is correctly formed without options
        mock_get.assert_called_once()
        called_url = mock_get.call_args[0][0]
        self.assertTrue(called_url.startswith('https://rest.ensembl.org/lookup/id/ENST00000380152'))

    @patch('VariantValidator.modules.utils.requests.get')
    def test_successful_request_grch37_with_options(self, mock_get):
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {'id': 'ENST00000380152', 'gene': 'BRCA2'}
        mock_get.return_value = mock_response

        result = utils.ensembl_rest(
            'ENST00000380152.4',
            '/lookup/id/',
            'GRCh37',
            options='expand=1'
        )
        self.assertEqual(result['record'], {'id': 'ENST00000380152', 'gene': 'BRCA2'})
        self.assertEqual(result['error'], 'false')
        mock_get.assert_called_once()
        called_url = mock_get.call_args[0][0]
        self.assertIn('expand=1', called_url)
        self.assertTrue(called_url.startswith('https://grch37.rest.ensembl.org/lookup/id/ENST00000380152'))

    @patch('VariantValidator.modules.utils.requests.get')
    def test_request_error(self, mock_get):
        # Simulate a failed request
        mock_response = MagicMock()
        mock_response.status_code = 500
        mock_get.return_value = mock_response

        result = utils.ensembl_rest('ENST00000380152.4', '/lookup/id/', 'GRCh38')
        self.assertTrue(result['error'].startswith('Problem encountered while connecting Ensembl REST'))
        self.assertEqual(result['record'], '')

class TestEnsemblTark(TestCase):

    @patch('VariantValidator.modules.utils.requests.get')
    def test_successful_request(self, mock_get):
        # Mock a successful JSON response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {'id': 'ENST00000380152', 'gene': 'BRCA2'}
        mock_get.return_value = mock_response

        result = utils.ensembl_tark('ENST00000380152.4', '/lookup/id/')
        self.assertEqual(result['record'], {'id': 'ENST00000380152', 'gene': 'BRCA2'})
        self.assertEqual(result['error'], 'false')

        # Ensure the request URL is correctly formed
        mock_get.assert_called_once()
        called_url = mock_get.call_args[0][0]
        self.assertTrue(called_url.startswith('https://tark.ensembl.org/lookup/id/'))
        self.assertIn('stable_id_with_version=ENST00000380152.4', called_url)

    @patch('VariantValidator.modules.utils.requests.get')
    def test_request_error_status(self, mock_get):
        # Simulate a failed request (non-200)
        mock_response = MagicMock()
        mock_response.status_code = 500
        mock_get.return_value = mock_response

        result = utils.ensembl_tark('ENST00000380152.4', '/lookup/id/')
        self.assertTrue(result['error'].startswith('Problem encountered while connecting Ensembl REST'))
        self.assertEqual(result['record'], '')

    @patch('VariantValidator.modules.utils.requests.get')
    def test_request_invalid_schema(self, mock_get):
        # Simulate requests raising InvalidSchema
        mock_get.side_effect = requests.exceptions.InvalidSchema("Invalid schema error")

        result = utils.ensembl_tark('ENST00000380152.4', '/lookup/id/')
        self.assertEqual(result['error'], "Invalid schema error")
        self.assertEqual(result['record'], '')

class TestRemoveReferenceString(TestCase):

    def test_deletion(self):
        input_var = 'NM_000123.4:c.34delG'
        output = utils.remove_reference_string(input_var)
        self.assertEqual(output, 'NM_000123.4:c.34del')

    def test_inversion(self):
        input_var = 'NM_000123.4:c.34invAT'
        output = utils.remove_reference_string(input_var)
        self.assertEqual(output, 'NM_000123.4:c.34inv')

    def test_no_change(self):
        input_var = 'NM_000123.4:c.34C>T'
        output = utils.remove_reference_string(input_var)
        self.assertEqual(output, input_var)

    def test_multiple_patterns(self):
        # Only last matching pattern should be changed
        input_var = 'NM_000123.4:c.34delAinsT'
        output = utils.remove_reference_string(input_var)
        self.assertEqual(output, 'NM_000123.4:c.34delinsT')

class TestUserInput(TestCase):

    def test_g_variant_no_brackets(self):
        query = "NM_000123.4:g.123A>T"
        result = utils.user_input(query)
        self.assertEqual(result, {'variant': query, 'type': ':g.'})

    def test_g_variant_with_brackets(self):
        query = "NM_000123.4:g.(gene)123A>T"
        result = utils.user_input(query)
        self.assertEqual(result, {'variant': "NM_000123.4:g.123A>T", 'type': ':g.'})

    def test_c_variant(self):
        query = "NM_000123.4:c.34C>T"
        result = utils.user_input(query)
        self.assertEqual(result, {'variant': query, 'type': ':c.'})

    def test_r_variant(self):
        query = "NM_000123.4:r.34C>U"
        result = utils.user_input(query)
        self.assertEqual(result, {'variant': query, 'type': ':r.'})

    def test_n_variant(self):
        query = "NM_000123.4:n.34C>T"
        result = utils.user_input(query)
        self.assertEqual(result, {'variant': query, 'type': ':n.'})

    def test_p_variant(self):
        query = "NP_000123.1:p.Val12Gly"
        result = utils.user_input(query)
        self.assertEqual(result, {'variant': query, 'type': ':p.'})

    def test_m_variant(self):
        query = "NC_000123.1:m.123A>G"
        result = utils.user_input(query)
        self.assertEqual(result, {'variant': query, 'type': ':m.'})

    def test_est_variant(self):
        query = "1:234"
        result = utils.user_input(query)
        self.assertEqual(result, {'variant': query, 'type': 'est'})

    def test_invalid_variant(self):
        query = "unknown_variant"
        result = utils.user_input(query)
        self.assertEqual(result, 'invalid')

    def test_whitespace_stripping(self):
        query = "  NM_000123.4:c.34C>T  "
        result = utils.user_input(query)
        self.assertEqual(result, {'variant': "NM_000123.4:c.34C>T", 'type': ':c.'})

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
