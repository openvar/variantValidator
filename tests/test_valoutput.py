import json
from unittest import TestCase
import VariantValidator
from VariantValidator.modules.valoutput import ValOutput
from VariantValidator.modules.variant import Variant


class TestValOutput(TestCase):
    """
    class will test the ValOutput object and methods
    """
    @classmethod
    def setUpClass(cls):
        cls.vv = VariantValidator.Validator()

    def test_creation_empty(self):
        obj = ValOutput([], self.vv)

        self.assertIsInstance(obj, ValOutput)
        self.assertEqual(obj.output_list, [])
        self.assertEqual(obj.validator, self.vv)

    def test_creation_str(self):
        obj = ValOutput('hello', 'validator')

        self.assertEqual(obj.output_list, 'hello')
        self.assertEqual(obj.validator, 'validator')
        with self.assertRaises(AttributeError):
            obj.format_as_dict()

    def test_output_meta(self):
        obj = ValOutput([], self.vv)
        res = obj.format_as_dict(with_meta=True)
        self.assertIsInstance(res, dict)
        self.assertEqual(list(res), ['flag', 'metadata'])
        self.assertIsInstance(res['metadata'], dict)
        self.assertEqual(list(res['metadata']), ['variantvalidator_version', 'variantvalidator_hgvs_version',
                                                 'vvta_version', 'vvseqrepo_db', 'vvdb_version'])

        res2 = obj.format_as_table(with_meta=True)
        self.assertIsInstance(res2, list)
        self.assertTrue(res2[0].startswith('#'))
        self.assertTrue('variantvalidator_version' in res2[0])

    def test_dict_no_variants(self):
        obj = ValOutput([], self.vv)
        res = obj.format_as_dict(with_meta=False)
        self.assertIsInstance(res, dict)
        self.assertEqual(res, {'flag': 'empty_result'})

    def test_dict_one_variant(self):
        var = Variant('')
        obj = ValOutput([var], self.vv)
        res = obj.format_as_dict(with_meta=False)
        self.assertIsInstance(res, dict)
        self.assertEqual(list(res), ['flag', 'validation_warning_1'])
        self.assertEqual(res['flag'], 'warning')

    def test_dict_two_variants(self):
        var1 = Variant('var1')
        var2 = Variant('var2')
        obj = ValOutput([var1, var2], self.vv)
        res = obj.format_as_dict(with_meta=False)
        self.assertIsInstance(res, dict)
        self.assertEqual(list(res), ['flag', 'validation_warning_1', 'validation_warning_2'])
        self.assertEqual(res['flag'], 'warning')
        self.assertEqual(res['validation_warning_1']['submitted_variant'], 'var1')
        self.assertEqual(res['validation_warning_2']['submitted_variant'], 'var2')

    def test_dict_one_good_variant(self):
        var1 = Variant('var1')
        var1.output_type_flag = 'gene'
        var2 = Variant('var2')
        obj = ValOutput([var1, var2], self.vv)
        res = obj.format_as_dict(with_meta=False)
        print(res)
        self.assertIsInstance(res, dict)
        self.assertEqual(list(res), ['flag', 'None', 'validation_warning_1'])
        self.assertEqual(res['flag'], 'gene_variant')
        self.assertEqual(res['None']['submitted_variant'], 'var1')
        self.assertEqual(res['validation_warning_1']['submitted_variant'], 'var2')

    def test_dict_one_intergenic(self):
        var1 = Variant('var1')
        var1.output_type_flag = 'intergenic'
        var2 = Variant('var2')
        obj = ValOutput([var1, var2], self.vv)
        res = obj.format_as_dict(with_meta=False)
        print(res)
        self.assertIsInstance(res, dict)
        self.assertEqual(list(res), ['flag', 'intergenic_variant_1', 'validation_warning_1'])
        self.assertEqual(res['flag'], 'intergenic')
        self.assertEqual(res['intergenic_variant_1']['submitted_variant'], 'var1')
        self.assertEqual(res['validation_warning_1']['submitted_variant'], 'var2')

    def test_dict_one_intergenic_and_one_gene(self):
        var1 = Variant('var1')
        var1.output_type_flag = 'intergenic'
        var2 = Variant('var2')
        var2.output_type_flag = 'gene'
        obj = ValOutput([var1, var2], self.vv)
        res = obj.format_as_dict(with_meta=False)
        print(res)
        self.assertIsInstance(res, dict)
        self.assertEqual(list(res), ['flag', 'intergenic_variant_1', 'None'])
        self.assertEqual(res['flag'], 'gene_variant')
        self.assertEqual(res['intergenic_variant_1']['submitted_variant'], 'var1')
        self.assertEqual(res['None']['submitted_variant'], 'var2')

    def test_dict_one_intergenic_and_one_gene_reversed(self):
        var1 = Variant('var1')
        var1.output_type_flag = 'intergenic'
        var2 = Variant('var2')
        var2.output_type_flag = 'gene'
        obj = ValOutput([var2, var1], self.vv)
        res = obj.format_as_dict(with_meta=False)
        print(res)
        self.assertIsInstance(res, dict)
        self.assertEqual(list(res), ['flag', 'None', 'intergenic_variant_1'])
        self.assertEqual(res['flag'], 'intergenic')
        self.assertEqual(res['intergenic_variant_1']['submitted_variant'], 'var1')
        self.assertEqual(res['None']['submitted_variant'], 'var2')

    def test_dict_each_with_error_in_warnings(self):
        var1 = Variant('var1')
        var1.warnings = ['Validation error']
        var1.output_type_flag = 'gene'
        var2 = Variant('var2')
        var2.warnings = ['Validation error']
        var2.output_type_flag = 'intergenic'
        var3 = Variant('var3')
        var3.warnings = ['Validation error']

        obj = ValOutput([var1, var2, var3], self.vv)
        res = obj.format_as_dict(with_meta=False)
        print(res)

        self.assertIsInstance(res, dict)
        self.assertEqual(list(res), ['flag', 'validation_error_1', 'intergenic_variant_1', 'validation_error_2'])
        self.assertEqual(res['flag'], 'intergenic')
        self.assertEqual(res['validation_error_1']['submitted_variant'], 'var1')
        self.assertEqual(res['intergenic_variant_1']['submitted_variant'], 'var2')
        self.assertEqual(res['validation_error_2']['submitted_variant'], 'var3')

    def test_dict_each_obsolete(self):
        var1 = Variant('var1')
        var1.warnings = ['obsolete']
        var1.output_type_flag = 'gene'
        var2 = Variant('var2')
        var2.warnings = ['obsolete']
        var2.output_type_flag = 'intergenic'
        var3 = Variant('var3')
        var3.warnings = ['obsolete']
        var4 = Variant('var4')
        var4.warnings = ['obsolete']
        var4.output_type_flag = 'gene'
        var4.hgvs_transcript_variant = ''

        obj = ValOutput([var1, var2, var3, var4], self.vv)
        res = obj.format_as_dict(with_meta=False)
        print(res)

        self.assertIsInstance(res, dict)
        self.assertEqual(list(res), ['flag', 'None', 'intergenic_variant_1', 'obsolete_record_1', 'obsolete_record_2'])
        self.assertEqual(res['flag'], 'gene_variant')
        self.assertEqual(res['None']['submitted_variant'], 'var1')
        self.assertEqual(res['intergenic_variant_1']['submitted_variant'], 'var2')
        self.assertEqual(res['obsolete_record_1']['submitted_variant'], 'var3')
        self.assertEqual(res['obsolete_record_2']['submitted_variant'], 'var4')

    def test_json(self):
        var = Variant('')
        obj = ValOutput([var], self.vv)
        res = obj.format_as_json(with_meta=False)
        self.assertIsInstance(res, str)
        self.assertIn('"flag": "warning"', res)
        self.assertIn('"validation_warning_1": {"selected_assembly": false, "submitted_variant": ""', res)
        self.assertEqual(json.loads(res), obj.format_as_dict(with_meta=False))

    def test_table_empty(self):
        obj = ValOutput([], self.vv)
        res = obj.format_as_table(with_meta=False)
        print(res)
        self.assertIsInstance(res, list)
        self.assertEqual(res, [['Input', 'Warnings', 'HGVS_transcript', 'HGVS_intronic_chr_context',
                                'HGVS_intronic_rsg_context', 'HGVS_RefSeqGene', 'HGVS_LRG',
                                'HGVS_LRG_transcript', 'HGVS_Predicted_Protein', 'HGVS_Genomic_GRCh37', 'GRCh37_CHR',
                                'GRCh37_POS', 'GRCh37_ID', 'GRCh37_REF', 'GRCh37_ALT', 'HGVS_Genomic_GRCh38',
                                'GRCh38_CHR', 'GRCh38_POS', 'GRCh38_ID', 'GRCh38_REF', 'GRCh38_ALT',
                                'Gene_Symbol', 'HGNC_Gene_ID', 'Transcript_description', 'Alt_genomic_loci']])

    def test_table_one(self):
        var1 = Variant('var1')
        obj = ValOutput([var1], self.vv)
        res = obj.format_as_table(with_meta=False)
        print(res)
        self.assertIsInstance(res, list)
        self.assertEqual(res, [['Input', 'Warnings', 'HGVS_transcript', 'HGVS_intronic_chr_context',
                                'HGVS_intronic_rsg_context', 'HGVS_RefSeqGene', 'HGVS_LRG',
                                'HGVS_LRG_transcript', 'HGVS_Predicted_Protein', 'HGVS_Genomic_GRCh37', 'GRCh37_CHR',
                                'GRCh37_POS', 'GRCh37_ID', 'GRCh37_REF', 'GRCh37_ALT', 'HGVS_Genomic_GRCh38',
                                'GRCh38_CHR', 'GRCh38_POS', 'GRCh38_ID', 'GRCh38_REF', 'GRCh38_ALT',
                                'Gene_Symbol', 'HGNC_Gene_ID', 'Transcript_description', 'Alt_genomic_loci'],
                               ['var1', '', None, None, None, None, None, None, '', '', '', '', '', '', '', '', '',
                                '', '', '', '', '', '', '', '']])

    def test_table_one_gene(self):
        var1 = Variant('var1')
        var1.output_type_flag = 'gene'
        obj = ValOutput([var1], self.vv)
        res = obj.format_as_table(with_meta=False)
        print(res)
        self.assertIsInstance(res, list)
        self.assertEqual(res[0], ['Input', 'Warnings', 'HGVS_transcript',
                                  'HGVS_intronic_chr_context',
                                  'HGVS_intronic_rsg_context', 'HGVS_RefSeqGene', 'HGVS_LRG',
                                  'HGVS_LRG_transcript', 'HGVS_Predicted_Protein', 'HGVS_Genomic_GRCh37', 'GRCh37_CHR',
                                  'GRCh37_POS', 'GRCh37_ID', 'GRCh37_REF', 'GRCh37_ALT', 'HGVS_Genomic_GRCh38',
                                  'GRCh38_CHR', 'GRCh38_POS', 'GRCh38_ID', 'GRCh38_REF', 'GRCh38_ALT',
                                  'Gene_Symbol', 'HGNC_Gene_ID', 'Transcript_description', 'Alt_genomic_loci'])
        self.assertEqual(res[1], ['var1', '', None, None, None, None, None, None, '', '', '', '', '', '', '', '', '', '',
                                  '', '', '', '', '', '', ''])
        self.assertEqual(len(res), 2)

    def test_table_intergenic(self):
        var1 = Variant('var1')
        var1.output_type_flag = 'intergenic'
        obj = ValOutput([var1], self.vv)
        res = obj.format_as_table(with_meta=False)
        print(res)
        self.assertIsInstance(res, list)
        self.assertEqual(res[1], ['var1', '', None, None, None, None, None, None,'', '', '', '', '', '', '', '', '',
                                  '', '', '', '', '', '', '', ''])
        self.assertEqual(len(res), 2)

    def test_table_gene_warnings(self):
        var1 = Variant('var1')
        var1.output_type_flag = 'gene'
        var1.warnings = ['Validation error']
        var2 = Variant('var2')
        var2.output_type_flag = 'gene'
        var2.warnings = ['obsolete']
        var3 = Variant('var3')
        var3.output_type_flag = 'gene'
        var3.warnings = ['obsolete']
        var3.hgvs_transcript_variant = ''

        obj = ValOutput([var1, var2, var3], self.vv)
        res = obj.format_as_table(with_meta=False)

        self.assertIsInstance(res, list)
        self.assertEqual(res[1], ['var1', 'Validation error', None, None, None, None, None, None, '', '', '', '', '',
                                  '', '', '', '', '', '', '', '', '', '', '', ''])
        self.assertEqual(res[2], ['var2', 'obsolete', None, None, None, None, None, None, '', '', '', '', '', '', '',
                                  '', '', '', '', '', '', '', '', '', ''])
        self.assertEqual(res[3], ['var3', 'obsolete', '', None, None, None, None, None, '', '', '', '', '', '', '', '',
                                  '', '', '', '', '', '', '', '', ''])
        self.assertEqual(len(res), 4)

    def test_table_intergenic_warnings(self):
        var1 = Variant('var1')
        var1.output_type_flag = 'intergenic'
        var1.warnings = ['Validation error']
        var2 = Variant('var2')
        var2.output_type_flag = 'intergenic'
        var2.warnings = ['obsolete']
        var3 = Variant('var3')
        var3.output_type_flag = 'intergenic'
        var3.warnings = ['obsolete']
        var3.hgvs_transcript_variant = ''

        obj = ValOutput([var1, var2, var3], self.vv)
        res = obj.format_as_table(with_meta=False)
        self.assertIsInstance(res, list)
        self.assertEqual(res[1], ['var1', 'Validation error', None, None, None, None, None, None, '', '', '', '', '',
                                  '', '', '', '', '', '', '', '', '', '', '', ''])
        self.assertEqual(res[2], ['var2', 'obsolete', None, None, None, None, None, None, '', '', '', '', '', '', '',
                                  '', '', '', '', '', '', '', '', '', ''])
        self.assertEqual(res[3], ['var3', 'obsolete', '', None, None, None, None, None, '', '', '', '', '', '', '',
                                  '', '', '', '', '', '', '', '', '', ''])
        self.assertEqual(len(res), 4)

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
