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
