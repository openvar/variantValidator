import unittest
import re
import VariantValidator


class TestGene2Transcripts(unittest.TestCase):
    """
    This class will test the gene2transcripts method of the validator
    """

    @classmethod
    def setUpClass(cls):
        cls.vv = VariantValidator.Validator()

    def test_empty(self):
        output = self.vv.gene2transcripts('')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'Please enter HGNC gene name or transcript identifier (NM_, NR_, or ENST)')

    def test_nonsense(self):
        output = self.vv.gene2transcripts('nonsense')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'Unable to recognise gene symbol NONSENSE')

    def test_nonsense_NM(self):
        output = self.vv.gene2transcripts('NM_nonsense')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'No transcript definition for (tx_ac=NM_NONSENSE)')

    def test_nonsense_NR(self):
        output = self.vv.gene2transcripts('nonNR_sense')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'No transcript definition for (tx_ac=NONNR_SENSE)')

    def test_nonsense_NM_dot(self):
        output = self.vv.gene2transcripts('NM_nonsens.e')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'No transcript definition for (tx_ac=NM_NONSENS)')

    def test_nonsense_NM_dot_orf(self):
        output = self.vv.gene2transcripts('NM_nonsense.1ORF2')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'No transcript definition for (tx_ac=NM_NONSENSE)')

    def test_nonsense_LRG(self):
        output = self.vv.gene2transcripts('LRG_nonsense')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'Unable to recognise gene symbol LRG_NONSENSE')

    def test_nonsense_LRGT(self):
        output = self.vv.gene2transcripts('LRGT_nonsense')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'Unable to recognise gene symbol LRGT_NONSENSE')

    def test_NM(self):
        output = self.vv.gene2transcripts('NM_024865.3')
        print(output)
        self.assertEqual(list(output), ['current_symbol', 'previous_symbol', 'current_name',
                                        'hgnc', 'transcripts'])
        self.assertEqual(output['current_symbol'], 'NANOG')
        self.assertTrue(len(output['transcripts']) > 2)
        for transcript in output['transcripts']:
            self.assertTrue(
                    re.match('NM_024865.',transcript['reference'])
                    or
                    re.match('NM_001297698.',transcript['reference'])
                    )

    def test_NM_noversion(self):
        output = self.vv.gene2transcripts('NM_024865')
        print(output)
        self.assertEqual(list(output), ['current_symbol', 'previous_symbol', 'current_name',
                                        'hgnc', 'transcripts'])
        self.assertEqual(output['current_symbol'], 'NANOG')
        self.assertTrue(len(output['transcripts']) > 2)
        for transcript in output['transcripts']:
            self.assertTrue(
                    re.match('NM_024865.',transcript['reference'])
                    or
                    re.match('NM_001297698.',transcript['reference'])
                    )

    def test_sym(self):
        output = self.vv.gene2transcripts('NANOG')
        print(output)
        self.assertEqual(list(output), ['current_symbol', 'previous_symbol', 'current_name',
                                        'hgnc', 'transcripts'])
        self.assertEqual(output['current_symbol'], 'NANOG')
        self.assertTrue(len(output['transcripts']) > 2)
        for transcript in output['transcripts']:
            self.assertTrue(
                    re.match('NM_024865.', transcript['reference'])
                    or
                    re.match('NM_001297698.', transcript['reference'])
                    )

    def test_old_sym(self):
        output = self.vv.gene2transcripts('OTF3')
        print(output)
        self.assertEqual(list(output), ['current_symbol', 'previous_symbol', 'current_name',
                                        'hgnc', 'transcripts'])
        self.assertEqual(output['current_symbol'], 'POU5F1')
        self.assertTrue(len(output['transcripts']) > 7)
        for transcript in output['transcripts']:
            ref = re.sub('\.\d+', '.', transcript['reference'])
            self.assertTrue(ref in [
                'NM_203289.','NM_001173531.','NM_002701.','NM_001285986.','NM_001285987.'
                ])

    def test_ens(self):
        output = self.vv.gene2transcripts('ENSG00000204531')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'Unable to recognise gene symbol ENSG00000204531')

    def test_current_single_previous(self):
        symbol = 'P3H1'
        results = self.vv.gene2transcripts(symbol)
        print(results)

        assert results["current_name"] == "prolyl 3-hydroxylase 1"
        assert results["current_symbol"] == "P3H1"
        assert results["hgnc"] == "HGNC:19316"
        assert results["previous_symbol"] == "LEPRE1"

    def test_previous_single_previous(self):
        symbol = 'LEPRE1'
        results = self.vv.gene2transcripts(symbol)
        print(results)

        assert results["current_name"] == "prolyl 3-hydroxylase 1"
        assert results["current_symbol"] == "P3H1"
        assert results["hgnc"] == "HGNC:19316"
        assert results["previous_symbol"] == "LEPRE1"

    def test_current_symbol_also_previous_symbol(self):
        symbol = 'HTT'  # http://rest.genenames.org/search/prev_symbol/HTT (SLC6A4)
        results = self.vv.gene2transcripts(symbol)
        print(results)

        assert results["current_name"] == "huntingtin"
        assert results["current_symbol"] == "HTT"
        assert results["hgnc"] == "HGNC:4851"
        assert results["previous_symbol"] == "HD"

    def test_current_symbol_multi_previous_symbol(self):
        symbol = 'SLC6A4'
        results = self.vv.gene2transcripts(symbol)
        print(results)

        assert results["current_name"] == "solute carrier family 6 member 4"
        assert results["current_symbol"] == "SLC6A4"
        assert results["hgnc"] == "HGNC:11050"
        assert results["previous_symbol"] == "HTT,OCD1"


# <LICENSE>
# Copyright (C) 2016-2021 VariantValidator Contributors
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
