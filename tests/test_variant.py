from unittest import TestCase
from VariantValidator.modules.variant import Variant
from VariantValidator.modules.utils import VariantValidatorError
TestCase.maxDiff = None


class TestCreation(TestCase):
    """
    Test the creation and attributes of the Variant obj.
    """

    def test_create(self):
        var = Variant('NM_015120.4:c.34=')

        self.assertIsInstance(var, Variant)
        self.assertEqual(var.original, 'NM_015120.4:c.34=')

    def test_create_all_init(self):
        var = Variant('NM_015120.4:c.34=', quibble='NM_015120.4:c.34=', warnings=['Got a warning'], write=False,
                      primary_assembly='GRCh37', order=1)

        self.assertEqual(var.quibble, 'NM_015120.4:c.34=')
        self.assertEqual(var.warnings, ['Got a warning'])
        self.assertFalse(var.write)
        self.assertEqual(var.primary_assembly, 'GRCh37')
        self.assertEqual(var.order, 1)

    def test_quibble_type(self):
        var = Variant('NM_015120.4:c.34=', quibble=0)
        self.assertEqual(var.quibble, 0)

    def test_quibble_not_set(self):
        var = Variant('NM_015120.4:c.34=')
        self.assertEqual(var.quibble, 'NM_015120.4:c.34=')

    def test_warnings_type(self):
        var = Variant('NM_015120.4:c.34=', warnings='string')
        self.assertEqual(var.warnings, ['string'])

    def test_warnings_type_int(self):
        var = Variant('NM_015120.4:c.34=', warnings=0)
        self.assertEqual(var.warnings, [0])

    def test_warnings_not_set(self):
        var = Variant('NM_015120.4:c.34=')
        self.assertEqual(var.warnings, [])

    def test_write_type(self):
        var = Variant('NM_015120.4:c.34=', write='banana')
        self.assertEqual(var.write, 'banana')

    def test_write_not_set(self):
        var = Variant('NM_015120.4:c.34=')
        self.assertEqual(var.write, True)

    def test_primary_assembly_not_set(self):
        var = Variant('NM_015120.4:c.34=')
        self.assertEqual(var.primary_assembly, False)

    def test_order_not_set(self):
        var = Variant('NM_015120.4:c.34=')
        self.assertEqual(var.order, False)

    def test_all_defaults(self):
        var = Variant('NM_015120.4:c.34=')

        self.assertEqual(var.hgvs_formatted, None)
        self.assertEqual(var.hgvs_genomic, None)
        self.assertEqual(var.hgvs_coding, None)
        self.assertEqual(var.post_format_conversion, None)
        self.assertEqual(var.pre_RNA_conversion, None)
        self.assertEqual(var.input_parses, None)

        self.assertEqual(var.description, '')
        self.assertEqual(var.coding, '')
        self.assertEqual(var.coding_g, '')
        self.assertEqual(var.genomic_r, '')
        self.assertEqual(var.genomic_g, '')
        self.assertEqual(var.protein, '')
        self.assertEqual(var.output_type_flag, 'warning')
        self.assertEqual(var.gene_symbol, '')

        self.assertEqual(var.timing, {})

        self.assertEqual(var.refsource, None)
        self.assertEqual(var.reftype, None)

        # Normalizers
        self.assertEqual(var.hn, None)
        self.assertEqual(var.reverse_normalizer, None)
        self.assertEqual(var.evm, None)
        self.assertEqual(var.no_norm_evm, None)
        self.assertEqual(var.min_evm, None)
        self.assertEqual(var.lose_vm,  None)

        # Required for output
        self.assertEqual(var.hgvs_transcript_variant, None)
        self.assertEqual(var.genome_context_intronic_sequence, None)
        self.assertEqual(var.refseqgene_context_intronic_sequence, None)
        self.assertEqual(var.hgvs_refseqgene_variant, None)
        self.assertEqual(var.hgvs_predicted_protein_consequence, None)
        self.assertEqual(var.hgvs_lrg_transcript_variant, None)
        self.assertEqual(var.hgvs_lrg_variant, None)
        self.assertEqual(var.alt_genomic_loci, None)
        self.assertEqual(var.primary_assembly_loci, None)
        self.assertEqual(var.reference_sequence_records, None)
        self.assertEqual(var.validated, False)


class TestMethods(TestCase):
    """ Test each method in the Variant Obj"""

    def setUp(self):
        self.var = Variant('NM_015120.4:c.34=')

    def test_is_ascii(self):
        self.assertTrue(self.var.is_ascii())

    def test_is_ascii_false(self):
        self.var.quibble = 'NM_015120.4:c.34=\u0086'
        self.assertFalse(self.var.is_ascii())

    def test_is_ascii_false2(self):
        self.var.quibble = 'NM_015120.4:c.34='
        self.assertFalse(self.var.is_ascii())

    def test_get_non_ascii(self):
        chars, pos = self.var.get_non_ascii()
        self.assertEqual(chars, [])
        self.assertEqual(pos, [])

    def test_get_non_ascii_encoded(self):
        self.var.quibble = 'NM_\u0086015120.4:c.34='
        chars, pos = self.var.get_non_ascii()
        self.assertEqual(chars, [''])
        self.assertEqual(pos, [4])

    def test_get_non_ascii_decoded(self):
        self.var.quibble = 'NM_015120.4:c.34='
        chars, pos = self.var.get_non_ascii()
        self.assertEqual(chars, [''])
        self.assertEqual(pos, [18])

    def test_get_non_ascii_pair(self):
        self.var.quibble = 'NM_\u0086015120.4:c.34='
        chars, pos = self.var.get_non_ascii()
        self.assertEqual(chars, ['', ''])
        self.assertEqual(pos, [4, 12])

    def test_remove_whitespace(self):
        self.var.remove_whitespace()
        self.assertEqual(self.var.quibble, 'NM_015120.4:c.34=')

    def test_remove_whitespace_space(self):
        self.var.quibble = 'NM_015120  .4:c. 34='
        self.var.remove_whitespace()
        self.assertEqual(self.var.quibble, 'NM_015120.4:c.34=')

    def test_remove_whitespace_tab(self):
        self.var.quibble = 'NM_015120.\t4:c.34  ='
        self.var.remove_whitespace()
        self.assertEqual(self.var.quibble, 'NM_015120.4:c.34=')

    def test_remove_whitespace_newline(self):
        self.var.quibble = 'NM_015120.4:c\n.34='
        self.var.remove_whitespace()
        self.assertEqual(self.var.quibble, 'NM_015120.4:c.34=')

    def test_format_quibble(self):
        output = self.var.format_quibble()
        self.assertEqual(self.var.quibble, 'NM_015120.4:c.34=')
        self.assertFalse(output)

    def test_format_quibble_brackets(self):
        self.var.quibble = 'NM_0151(REMOVE)20.4:c.34='
        output = self.var.format_quibble()
        self.assertEqual(self.var.quibble, 'NM_015120.4:c.34=')
        self.assertFalse(output)

    def test_format_quibble_source_fail(self):
        self.var.quibble = 'F_015120.4:c.34='
        output = self.var.format_quibble()
        self.assertTrue(output)
        self.assertEqual(self.var.quibble, 'F_015120.4:c.34=')

    def test_format_quibble_type_fail(self):
        self.var.quibble = 'NM_015120.4:w.34='
        output = self.var.format_quibble()
        self.assertTrue(output)
        self.assertEqual(self.var.quibble, 'NM_015120.4:w.34=')

    def test_set_reftype(self):
        self.var.set_reftype()
        self.assertEqual(self.var.reftype, ':c.')

    def test_set_reftype_rna(self):
        self.var.quibble = 'NM_015120.4:r.34='
        self.var.set_reftype()
        self.assertEqual(self.var.reftype, ':r.')

    def test_set_reftype_nucl(self):
        self.var.quibble = 'NM_015120.4:n.34='
        self.var.set_reftype()
        self.assertEqual(self.var.reftype, ':n.')

    def test_set_reftype_mito(self):
        self.var.quibble = 'NM_015120.4:m.34='
        self.var.set_reftype()
        self.assertEqual(self.var.reftype, ':m.')

    def test_set_reftype_genome(self):
        self.var.quibble = 'NM_015120.4:g.34='
        self.var.set_reftype()
        self.assertEqual(self.var.reftype, ':g.')

    def test_set_reftype_prot(self):
        self.var.quibble = 'NM_015120.4:p.34='
        self.var.set_reftype()
        self.assertEqual(self.var.reftype, ':p.')

    def test_set_reftype_est(self):
        self.var.quibble = 'NM_015120.4:3.34='
        self.var.set_reftype()
        self.assertEqual(self.var.reftype, 'est')

    def test_set_reftype_none(self):
        self.var.quibble = 'NM_015120.4:.34='
        with self.assertRaises(VariantValidatorError):
            self.var.set_reftype()

    def test_set_source(self):
        self.var.set_refsource()
        self.assertEqual(self.var.refsource, 'RefSeq')

    def test_set_source_refseq_min(self):
        self.var.quibble = 'N'
        self.var.set_refsource()
        self.assertEqual(self.var.refsource, 'RefSeq')

    def test_set_source_lrg(self):
        self.var.quibble = 'LRG_015120.4:c.34='
        self.var.set_refsource()
        self.assertEqual(self.var.refsource, 'LRG')

    def test_set_source_lrg_min(self):
        self.var.quibble = 'LRG'
        self.var.set_refsource()
        self.assertEqual(self.var.refsource, 'LRG')

    def test_set_source_ens(self):
        self.var.quibble = 'ENSG_015120.4:c.34='
        self.var.set_refsource()
        self.assertEqual(self.var.refsource, 'ENS')

    def test_set_source_ens_min(self):
        self.var.quibble = 'ENS'
        self.var.set_refsource()
        self.assertEqual(self.var.refsource, 'ENS')

    def test_set_source_none(self):
        self.var.quibble = 'SOMETHING ELSE'
        with self.assertRaises(VariantValidatorError):
            self.var.set_refsource()

    def test_set_quibble(self):
        self.var.set_quibble('New:c.var')
        self.assertEqual(self.var.quibble, 'New:c.var')
        self.assertEqual(self.var.refsource, 'RefSeq')
        self.assertEqual(self.var.reftype, ':c.')

    def test_is_obsolete(self):
        self.assertFalse(self.var.is_obsolete())

    def test_is_obsolete_false(self):
        self.var.warnings = ['Nearly obso', 'lete']
        self.assertFalse(self.var.is_obsolete())

    def test_is_obsolete_true(self):
        self.var.warnings = ['obsoleteANDother']
        self.assertTrue(self.var.is_obsolete())

    def test_process_warnings(self):
        output = self.var.process_warnings()
        self.assertIsInstance(output, list)
        self.assertEqual(output, [])

    def test_process_warnings_sub(self):
        self.var.warnings = ['variantdelATGCTAGCTA']
        output = self.var.process_warnings()
        self.assertEqual(output, ['variantdel'])

    def test_process_warnings_sub_not(self):
        self.var.warnings = ['variantdelATG']
        output = self.var.process_warnings()
        self.assertEqual(output, ['variantdelATG'])

    def test_process_warnings_strip(self):
        self.var.warnings = [' warning ']
        output = self.var.process_warnings()
        self.assertEqual(output, ['warning'])

    def test_process_warnings_replace(self):
        self.var.warnings = ['\'warning\'']
        output = self.var.process_warnings()
        self.assertEqual(output, ['warning'])

    def test_process_warnings_unique(self):
        self.var.warnings = ['one', 'two', 'one']
        output = self.var.process_warnings()
        self.assertEqual(output, ['one', 'two'])

    def test_output_dict_empty(self):
        output = self.var.output_dict()
        self.assertIsInstance(output, dict)
        self.assertEqual(output, {
            'submitted_variant': 'NM_015120.4:c.34=',
            'gene_ids': None,
            'gene_symbol': '',
            'transcript_description': '',
            'annotations': '',
            'hgvs_transcript_variant': None,
            'genome_context_intronic_sequence': None,
            'refseqgene_context_intronic_sequence': None,
            'hgvs_refseqgene_variant': None,
            'hgvs_predicted_protein_consequence': None,
            'validation_warnings': [],
            'hgvs_lrg_transcript_variant': None,
            'hgvs_lrg_variant': None,
            'alt_genomic_loci': None,
            'primary_assembly_loci': None,
            'reference_sequence_records': None,
            'selected_assembly': False
        })

    def test_output_dict_set(self):
        self.var.gene_symbol = 'Symbol'
        self.var.annotations = 'annotated'
        self.var.description = 'Desc'
        self.var.stable_gene_ids = 'My_id'
        self.var.hgvs_transcript_variant = 'hgvsvar'
        self.var.genome_context_intronic_sequence = 'gintronic'
        self.var.refseqgene_context_intronic_sequence = 'rintronic'
        self.var.hgvs_refseqgene_variant = 'refseq'
        self.var.hgvs_predicted_protein_consequence = 'prot'
        self.var.warnings = ['warning']
        self.var.hgvs_lrg_transcript_variant = 'lrgT'
        self.var.hgvs_lrg_variant = 'lrg'
        self.var.alt_genomic_loci = 'alt'
        self.var.primary_assembly_loci = 'primary'
        self.var.reference_sequence_records = 'records'
        self.var.selected_assembly = 'assembly'
        output = self.var.output_dict()
        self.assertIsInstance(output, dict)
        self.assertEqual(output, {
            'submitted_variant': 'NM_015120.4:c.34=',
            'gene_symbol': 'Symbol',
            'annotations': 'annotated',
            'gene_ids': 'My_id',
            'transcript_description': 'Desc',
            'hgvs_transcript_variant': 'hgvsvar',
            'genome_context_intronic_sequence': 'gintronic',
            'refseqgene_context_intronic_sequence': 'rintronic',
            'hgvs_refseqgene_variant': 'refseq',
            'hgvs_predicted_protein_consequence': 'prot',
            'validation_warnings': ['warning'],
            'hgvs_lrg_transcript_variant': 'lrgT',
            'hgvs_lrg_variant': 'lrg',
            'selected_assembly': 'assembly',
            'alt_genomic_loci': 'alt',
            'primary_assembly_loci': 'primary',
            'reference_sequence_records': 'records',
        })

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
