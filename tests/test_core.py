import unittest
import VariantValidator


class TestValidator(unittest.TestCase):
    """
    Going to test the Validator function with a series of different inputs/situations that aren't covered in
    test_inputs.py
    """

    @classmethod
    def setUpClass(cls):
        cls.vv = VariantValidator.Validator()

    def test_transcript_seq_nonsense(self):
        var = 'NM_015120.4:c.34C>T'
        with self.assertRaises(Exception):
            self.vv.validate(var, 'GRCh37', 'all', transcript_set='nonsense')

    def test_transcript_seq_ensembl(self):
        var = 'NM_015120.4:c.34C>T'
        with self.assertRaises(Exception):
            self.vv.validate(var, 'GRCh37', 'all', transcript_set='ensembl')

        self.assertEqual(self.vv.alt_aln_method, 'genebuild')

    def test_transcript_list(self):
        var = 'NM_015120.4:c.34C>T'

        output = self.vv.validate(var, 'GRCh37', 'Trans1').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'empty_result')

    def test_transcript_list_realid(self):
        var = 'NM_015120.4:c.34C>T'

        output = self.vv.validate(var, 'GRCh37', 'NM_015120.4').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_transcript_list_real_pair(self):
        var = 'NM_015120.4:c.34C>T'

        output = self.vv.validate(var, 'GRCh37', 'NM_015120.4|NM_015120.5').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_transcript_list_lrg(self):
        var = 'NM_015120.4:c.34C>T'

        output = self.vv.validate(var, 'GRCh37', 'LRG1').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'empty_result')

    def test_non_ascii(self):
        var = 'NM_015120.4:c.34C>T\202'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        self.assertIn('Submitted variant description contains an invalid character',
                      str(output['validation_warning_1']['validation_warnings']))

    def test_assembly_hg19(self):
        var = 'NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'hg19', 'all')
        for variant in out.output_list:
            self.assertEqual(variant.primary_assembly, 'GRCh37')
        output = out.format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_assembly_hg38(self):
        var = 'NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'hg38', 'all')
        for variant in out.output_list:
            self.assertEqual(variant.primary_assembly, 'GRCh38')
        output = out.format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_assembly_grch(self):
        var = 'NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'grch37', 'all')
        for variant in out.output_list:
            self.assertEqual(variant.primary_assembly, 'GRCh37')
        output = out.format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])

    def test_assembly_invalid(self):
        var = 'NM_015120.4:c.34C>T'

        out = self.vv.validate(var, 'nonsense', 'all')
        for variant in out.output_list:
            self.assertEqual(variant.primary_assembly, 'GRCh38')
        output = out.format_as_dict()
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', 'NM_015120.4:c.34C>T', 'metadata'])
        self.assertIn('Invalid genome build has been specified',
                      str(output['NM_015120.4:c.34C>T']['validation_warnings']))

    def test_variant_invalid(self):
        var = 'NM_015120.4c.34C>T'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        self.assertIn('Unable to identify a colon (:) in the variant description',
                      str(output['validation_warning_1']['validation_warnings']))

    def test_variant_invalid_2(self):
        var = 'NM_015120.4:c34C>T'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        self.assertIn('lacks the . character between <type> and <position>',
                      str(output['validation_warning_1']['validation_warnings']))

    def test_variant_invalid_3(self):
        var = 'nonsense'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        self.assertIn('Variant description nonsense is not in an accepted format',
                      str(output['validation_warning_1']['validation_warnings']))

    def test_variant_con(self):
        var = 'NM_015120.4:c.34con'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'warning')
        self.assertIn('Gene conversions currently unsupported',
                      str(output['validation_warning_1']['validation_warnings']))

    def test_variant_RNA(self):
        # TODO: This situation needs looking at as I'm sure it shouldn't be returning an empty string.
        var = 'NM_015120.4:r.34DEL'

        output = self.vv.validate(var, 'GRCh37', 'all').format_as_dict()
        print(output)
        self.assertEqual(output['flag'], 'gene_variant')
        self.assertEqual(list(output), ['flag', '', 'metadata'])


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
        self.assertEqual(output['error'], 'No transcript definition for (tx_ac=NM_NONSENS.E)')

    def test_nonsense_NM_dot_orf(self):
        output = self.vv.gene2transcripts('NM_nonsense.1ORF2')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'No transcript definition for (tx_ac=NM_NONSENSE.1orf2)')

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
                                        'previous_name', 'transcripts'])
        self.assertEqual(output['current_symbol'], 'NANOG')
        self.assertEqual(len(output['transcripts']), 3)

    def test_NM_noversion(self):
        output = self.vv.gene2transcripts('NM_024865')
        print(output)
        self.assertEqual(list(output), ['current_symbol', 'previous_symbol', 'current_name',
                                        'previous_name', 'transcripts'])
        self.assertEqual(output['current_symbol'], 'NANOG')
        self.assertEqual(len(output['transcripts']), 3)

    def test_sym(self):
        output = self.vv.gene2transcripts('NANOG')
        print(output)
        self.assertEqual(list(output), ['current_symbol', 'previous_symbol', 'current_name',
                                        'previous_name', 'transcripts'])
        self.assertEqual(output['current_symbol'], 'NANOG')
        self.assertEqual(len(output['transcripts']), 3)

    def test_old_sym(self):
        output = self.vv.gene2transcripts('OTF3')
        print(output)
        self.assertEqual(list(output), ['current_symbol', 'previous_symbol', 'current_name',
                                        'previous_name', 'transcripts'])
        self.assertEqual(output['current_symbol'], 'POU5F1')
        self.assertEqual(len(output['transcripts']), 8)

    def test_ens(self):
        output = self.vv.gene2transcripts('ENSG00000204531')
        print(output)
        self.assertEqual(list(output), ['error'])
        self.assertEqual(output['error'], 'Unable to recognise gene symbol ENSG00000204531')


class TestHGVS2Ref(unittest.TestCase):
    """
    class will test the inputs for the hgvs2ref method of the validator()
    """

    @classmethod
    def setUpClass(cls):
        cls.vv = VariantValidator.Validator()

    def test_empty(self):
        output = self.vv.hgvs2ref('')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertEqual(output['error'], ': char 1: end of input')

    def test_nonsense(self):
        output = self.vv.hgvs2ref('nonsense')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertEqual(output['error'], 'nonsense: char 9: end of input')

    def test_nonsense_colon(self):
        output = self.vv.hgvs2ref('non:sense')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertEqual(output['error'],
                         'non:sense: char 4: expected one of \'c\', \'g\', \'m\', \'n\', \'p\', or \'r\'')

    def test_nonsense_hgvs(self):
        output = self.vv.hgvs2ref('nonsense:c.34C>T')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertTrue('Failed to fetch nonsense from SeqRepo' in output['error'])

    def test_valid_c(self):
        output = self.vv.hgvs2ref('NM_015120.4:c.34C>T')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertEqual(output['error'], '')
        self.assertEqual(output['start_position'], '34')
        self.assertEqual(output['sequence'], 'C')

    def test_valid_g(self):
        output = self.vv.hgvs2ref('NM_015120.4:g.34C>T')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertEqual(output['error'], '')
        self.assertEqual(output['start_position'], '34')
        self.assertEqual(output['sequence'], 'A')

    def test_valid_n(self):
        output = self.vv.hgvs2ref('NM_015120.4:n.34C>T')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertEqual(output['error'], '')
        self.assertEqual(output['start_position'], '34')
        self.assertEqual(output['sequence'], 'A')

    def test_valid_p(self):
        output = self.vv.hgvs2ref('NM_015120.4:p.Thr34=')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertEqual(output['error'], '')
        self.assertEqual(output['start_position'], 'Thr34')
        self.assertEqual(output['sequence'], 'A')

    def test_valid_m(self):
        output = self.vv.hgvs2ref('NM_015120.4:m.34C>T')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertEqual(output['error'], '')
        self.assertEqual(output['start_position'], '')
        self.assertEqual(output['sequence'], '')

    def test_valid_r(self):
        output = self.vv.hgvs2ref('NM_015120.4:r.34C>U')
        print(output)
        self.assertEqual(list(output), ['variant', 'start_position', 'end_position', 'warning', 'sequence', 'error'])
        self.assertEqual(output['error'], '')
        self.assertEqual(output['start_position'], '')
        self.assertEqual(output['sequence'], '')

# <LICENSE>
# Copyright (C) 2019 VariantValidator Contributors
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