import re
from . import utils as fn


class Variant(object):
    """
    This Variant object will contain the original input, the processed variant description and any other data that's
    relevant to what kind of variant it is.
    """

    def __init__(self, original, quibble=None, warnings=None, write=True, primary_assembly=False, order=False,
                 selected_assembly=False):
        self.original = original
        if quibble is None:
            self.quibble = original
        else:
            self.quibble = quibble
        self.hgvs_formatted = None
        self.hgvs_genomic = None
        self.hgvs_coding = None
        self.post_format_conversion = None  # Used for first gapped_mapping function
        self.pre_RNA_conversion = None
        self.input_parses = None  # quibble as hgvs variant object

        if warnings is None:
            self.warnings = []
        else:
            if isinstance(warnings, list):
                self.warnings = warnings
            else:
                self.warnings = [warnings]
        self.description = ''  # hgnc_gene_info variable
        self.annotations = ''
        self.coding = ''
        self.coding_g = ''
        self.genomic_r = ''
        self.genomic_g = ''
        self.protein = ''
        self.write = write
        self.primary_assembly = primary_assembly
        self.selected_assembly = selected_assembly
        self.order = order
        self.output_type_flag = 'warning'
        self.gene_symbol = ''

        self.timing = {}

        self.refsource = None
        self.reftype = None

        # Normalizers
        self.hn = None
        self.reverse_normalizer = None
        self.evm = None
        self.no_norm_evm = None
        self.min_evm = None
        self.lose_vm = None

        # Required for output
        self.stable_gene_ids = None
        self.hgvs_transcript_variant = None  # variant.coding but edited
        self.genome_context_intronic_sequence = None
        self.refseqgene_context_intronic_sequence = None
        self.hgvs_refseqgene_variant = None  # genomic_r but edited
        self.hgvs_predicted_protein_consequence = None
        self.hgvs_lrg_transcript_variant = None
        self.hgvs_lrg_variant = None  # Same as hgvs_refseqgene_variant but with LRG accession
        self.alt_genomic_loci = None
        self.primary_assembly_loci = None
        self.reference_sequence_records = None
        self.validated = False

    def is_ascii(self):
        """
        Instead of the previous test for unicode rich text characters.
        Now going to test that all characters are within the ascii alphabet
        """
        try:
            self.quibble.encode('ascii')
            return True
        except UnicodeEncodeError or UnicodeDecodeError:
            # Will catch errors raised by python 2 and python 3
            return False

    def get_non_ascii(self):
        """
        Will return non ascii character positions within variant description
        :return:
        """
        chars = []
        positions = []

        for i, c in enumerate(self.quibble):
            try:
                c.encode('ascii')
            except UnicodeEncodeError or UnicodeDecodeError:
                chars.append(c)
                positions.append(i+1)

        return chars, positions

    def remove_whitespace(self):
        """
        Will remove all whitespace from quibble
        :return:
        """
        prev = self.quibble
        self.quibble = ''.join(self.quibble.split())
        if self.quibble != prev:
            caution = 'Whitespace removed from variant description %s' % self.quibble
            self.warnings.append(caution)

    def remove_quotes(self):
        if self.quibble.startswith('"') or self.quibble.startswith("'"):
            self.quibble = self.quibble[1:]
        if self.quibble.endswith('"') or self.quibble.endswith("'"):
            self.quibble = self.quibble[:-1]

    def format_quibble(self):
        """
        Removes whitespace from the ends of the string
        Removes anything in brackets
        Identifies variant type (p. c. etc)
        Accepts c, g, n, r currently. And now P also 15.07.15
        """
        # Set regular expressions for if statements
        pat_gene = re.compile(r'\(.+?\)')  # Pattern looks for (....)

        if pat_gene.search(self.quibble):
            self.quibble = pat_gene.sub('', self.quibble)

        try:
            self.set_refsource()
        except fn.VariantValidatorError:
            return True

        try:
            self.set_reftype()
        except fn.VariantValidatorError:
            return True

        return False

    def set_reftype(self):
        """
        Method will set the reftype based on the quibble
        :return:
        """
        pat_est = re.compile(r'\d:\d')

        if ':g.' in self.quibble:
            self.reftype = ':g.'
        elif ':r.' in self.quibble:
            self.reftype = ':r.'
        elif ':n.' in self.quibble:
            self.reftype = ':n.'
        elif ':c.' in self.quibble:
            self.reftype = ':c.'
        elif ':p.' in self.quibble:
            self.reftype = ':p.'
        elif ':m.' in self.quibble:
            self.reftype = ':m.'
        elif pat_est.search(self.quibble):
            self.reftype = 'est'
        else:
            raise fn.VariantValidatorError("Unable to identity reference type from %s" % self.quibble)

    def set_refsource(self):
        """
        Method will set the refsource based on the quibble
        :return:
        """
        if self.quibble.startswith('LRG'):
            self.refsource = 'LRG'
        elif self.quibble.startswith('ENS'):
            self.refsource = 'ENS'
        elif self.quibble.startswith('N'):
            self.refsource = 'RefSeq'
        else:
            raise fn.VariantValidatorError("Unable to identify reference source from %s" % self.quibble)

    def set_quibble(self, newval):
        """
        Method will set the quibble and reset the refsource and reftype
        :param newval:
        :return:
        """
        self.quibble = newval
        self.set_refsource()
        self.set_reftype()

    def output_dict(self, test=False):
        """
        Method will return the output values as a dictionary
        :return: dict
        """
        if test is True:
            try:
                del self.stable_gene_ids['ensembl_gene_id']
                del self.stable_gene_ids['ccds_ids']
            except KeyError:
                pass
            try:
                del self.hgvs_predicted_protein_consequence['lrg_tlr']
                del self.hgvs_predicted_protein_consequence['lrg_slr']
            except KeyError:
                pass

        dict_out = {
            'selected_assembly': self.selected_assembly,
            'submitted_variant': self.original,
            'gene_symbol': self.gene_symbol,
            'gene_ids': self.stable_gene_ids,
            'annotations': self.annotations,
            'transcript_description': self.description,
            'hgvs_transcript_variant': self.hgvs_transcript_variant,
            'genome_context_intronic_sequence': self.genome_context_intronic_sequence,
            'refseqgene_context_intronic_sequence': self.refseqgene_context_intronic_sequence,
            'hgvs_refseqgene_variant': self.hgvs_refseqgene_variant,
            'hgvs_predicted_protein_consequence': self.hgvs_predicted_protein_consequence,
            'validation_warnings': self.process_warnings(),
            'hgvs_lrg_transcript_variant': self.hgvs_lrg_transcript_variant,
            'hgvs_lrg_variant': self.hgvs_lrg_variant,
            'alt_genomic_loci': self.alt_genomic_loci,
            'primary_assembly_loci': self.primary_assembly_loci,
            'reference_sequence_records': self.reference_sequence_records,
        }
        return dict_out

    def is_obsolete(self):
        """
        Checks whether the keyword 'obsolete' appears within the validation warnings
        :return:
        """
        return any('obsolete' in warning for warning in self.warnings)

    def process_warnings(self):
        refined = []
        for warning in self.warnings:
            warning = re.sub('del[GATC][GATC][GATC][GATC]+', 'del', warning)
            warning = warning.strip()
            warning = warning.replace("'", "")
            if warning == '':
                continue
            if warning not in refined:
                refined.append(warning)
        return refined

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
