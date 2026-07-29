import re

from . import utils as fn
from .transcript_map_data import TranscriptMapData


class Variant(object):
    """
    This Variant object will contain the original input, the processed variant
    description and any other data that's relevant to what kind of variant it is.
    """

    def __init__(
            self,
            original,
            quibble=None,
            warnings=None,
            write=True,
            primary_assembly=False,
            order=False,
            selected_assembly=False,
            reformat_output=False,
            expanded_repeat=None
    ):
        self.original = original
        self.quibble = original if quibble is None else quibble

        self.hgvs_formatted = None
        self.hgvs_genomic = None
        self.hgvs_coding = None
        self.post_format_conversion = None
        self.pre_RNA_conversion = None
        self.input_parses = None
        self.transcript_type = None
        self.lovd_syntax_check = None
        self.shorthand_vcf = None
        self.lovd_messages = None
        self.lovd_corrections = None

        # Placeholder for alt_reference
        self.genomic_context_ac = None

        if warnings is None:
            self.warnings = []
        elif isinstance(warnings, list):
            self.warnings = warnings
        else:
            self.warnings = [warnings]

        self.description = ''
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
        self.expanded_repeat = expanded_repeat

        # Set reformat options
        self.reformat_output = reformat_output

        # Normalizers and mappers
        self.hn = None
        self.reverse_normalizer = None
        self.cross_hn = None
        self.evm = None
        self.no_norm_evm = None
        self.min_evm = None
        self.lose_vm = None
        self.no_replace_vm = None
        self.map_dat = TranscriptMapData()

        # Required for output
        self.stable_gene_ids = None
        self.hgvs_transcript_variant = None
        self.genome_context_intronic_sequence = None
        self.refseqgene_context_intronic_sequence = None
        self.hgvs_refseqgene_variant = None
        self.hgvs_predicted_protein_consequence = None
        self.hgvs_lrg_transcript_variant = None
        self.hgvs_lrg_variant = None
        self.alt_genomic_loci = None
        self.primary_assembly_loci = None
        self.reference_sequence_records = None
        self.validated = False
        self.exonic_positions = None
        self.rna_data = None

    def is_ascii(self):
        """
        Test that all characters in quibble are ASCII.
        """
        try:
            self.quibble.encode('ascii')
            return True
        except (UnicodeEncodeError, UnicodeDecodeError):
            return False

    def get_non_ascii(self):
        """
        Return non-ASCII characters and their positions within the variant
        description.
        """
        chars = []
        positions = []

        for i, char in enumerate(self.quibble):
            try:
                char.encode('ascii')
            except (UnicodeEncodeError, UnicodeDecodeError):
                chars.append(char)
                positions.append(i + 1)

        return chars, positions

    def remove_whitespace(self):
        """
        Remove all whitespace from quibble.
        """
        previous = self.quibble
        self.quibble = ''.join(self.quibble.split())

        if self.quibble != previous:
            caution = (
                'VariantSyntaxError: Whitespace removed from variant '
                f'description {self.original}'
            )
            self.warnings.append(caution)

    def remove_quotes(self):
        if self.quibble.startswith(('"', "'")):
            self.quibble = self.quibble[1:]

        if self.quibble.endswith(('"', "'")):
            self.quibble = self.quibble[:-1]

    def non_alphanum_start(self):
        """
        Check for an invalid leading character after removing whitespace
        and surrounding quotes.
        """
        if not re.search(r'^\w', self.original):
            self.remove_whitespace()
            self.remove_quotes()

            if not re.search(r'^\w', self.quibble):
                return True

        return False

    def format_quibble(self):
        """
        Remove formatting errors and identify the reference source/type.
        """
        try:
            self.set_refsource()
        except fn.VariantValidatorError:
            return True

        try:
            self.set_reftype()
        except fn.VariantValidatorError:
            return True

        # Upper case characters in edit types, e.g. Ins or dUP.
        # delins must precede del.
        edit_type_patterns = ('delins', 'dup', 'ins', 'del')

        for pattern in edit_type_patterns:
            matches = re.findall(pattern, self.quibble, re.IGNORECASE)

            for match in matches:
                if match != match.lower():
                    replacement = match.lower()
                    self.warnings.append(
                        f'Edit type {match} should be in the lower case, '
                        f'i.e. {replacement}'
                    )
                    self.quibble = self.quibble.replace(
                        match,
                        replacement
                    )

        return False

    def set_reftype(self):
        """
        Set the reference type based on quibble.
        """
        if not isinstance(self.quibble, str):
            reftype = self.quibble.type

            if reftype in ('g', 'r', 'n', 'c', 'p', 'm'):
                self.reftype = f':{reftype}.'
                return True

            raise fn.VariantValidatorError(
                "Unable to identity reference type from "
                f"{self.quibble}"
            )

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
        elif re.search(r'\d:\d', self.quibble):
            self.reftype = 'est'
        else:
            raise fn.VariantValidatorError(
                "Unable to identity reference type from "
                f"{self.quibble}"
            )

    def set_refsource(self):
        """
        Set the reference source based on quibble.
        """
        if isinstance(self.quibble, str):
            ac_testval = self.quibble
        else:
            ac_testval = self.quibble.ac

        if ac_testval.startswith('LRG'):
            self.refsource = 'LRG'
        elif ac_testval.startswith('ENS'):
            self.refsource = 'ENS'
        elif ac_testval.startswith('N'):
            self.refsource = 'RefSeq'
        else:
            raise fn.VariantValidatorError(
                "Unable to identify reference source from "
                f"{self.quibble}"
            )

    def set_quibble(self, newval):
        """
        Set quibble and reset the reference source and reference type.
        """
        self.quibble = newval
        self.set_refsource()
        self.set_reftype()

    def output_dict(self, test=False):
        """
        Return the output values as a dictionary.
        """
        if test is True:
            try:
                del self.stable_gene_ids['ensembl_gene_id']
                del self.stable_gene_ids['ccds_ids']
            except (KeyError, TypeError):
                pass

            try:
                del self.hgvs_predicted_protein_consequence['lrg_tlr']
                del self.hgvs_predicted_protein_consequence['lrg_slr']
            except (KeyError, TypeError):
                pass

        return {
            'selected_assembly': self.selected_assembly,
            'submitted_variant': self.original,
            'gene_symbol': self.gene_symbol,
            'gene_ids': self.stable_gene_ids,
            'annotations': self.annotations,
            'transcript_description': self.description,
            'hgvs_transcript_variant': self.hgvs_transcript_variant,
            'rna_variant_descriptions': self.rna_data,
            'genome_context_intronic_sequence':
                self.genome_context_intronic_sequence,
            'refseqgene_context_intronic_sequence':
                self.refseqgene_context_intronic_sequence,
            'hgvs_refseqgene_variant': self.hgvs_refseqgene_variant,
            'hgvs_predicted_protein_consequence':
                self.hgvs_predicted_protein_consequence,
            'validation_warnings': self.process_warnings(),
            'lovd_messages': self.lovd_messages,
            'lovd_corrections': self.lovd_corrections,
            'hgvs_lrg_transcript_variant':
                self.hgvs_lrg_transcript_variant,
            'hgvs_lrg_variant': self.hgvs_lrg_variant,
            'alt_genomic_loci': self.alt_genomic_loci,
            'primary_assembly_loci': self.primary_assembly_loci,
            'variant_exonic_positions': self.exonic_positions,
            'reference_sequence_records': self.reference_sequence_records
        }

    def is_obsolete(self):
        """
        Check whether 'obsolete' appears in the validation warnings.
        """
        return any(
            'obsolete' in warning
            for warning in self.warnings
        )

    def process_warnings(self, string_all=False):
        """
        Remove duplicate warnings and normalise warning strings.
        """
        refined = []

        for warning in self.warnings:
            if isinstance(warning, dict) and string_all is False:
                processed_warning = warning
            else:
                processed_warning = re.sub(
                    r'del[GATC]{4,}',
                    'del',
                    str(warning)
                )
                processed_warning = processed_warning.strip()
                processed_warning = processed_warning.replace("'", "")

                if not processed_warning:
                    continue

            if processed_warning not in refined:
                refined.append(processed_warning)

        return refined

    def remove_typos(self):
        """
        Remove an expanding list of common typos from the variant description.
        """
        # Double or multiple colons.
        if re.search(r'::+[cgpnr]\.', self.quibble):
            self.warnings.append(
                "VariantSyntaxError: Multiple colons found in variant "
                "description"
            )
            self.quibble = re.sub(r'::+', ':', self.quibble)

        # Missing colon.
        if (
                re.search(r'[gcrnpmo]\.', self.quibble)
                and not re.search(r':[gcrnpmo]\.', self.quibble)
        ):
            error = (
                "VariantSyntaxError: Unable to identify a colon (:) in the "
                f"variant description {self.quibble}. A colon is required in "
                "HGVS variant descriptions to separate the reference "
                "accession from the reference type i.e. <accession>:<type>. "
                "e.g. :c."
            )
            self.warnings.append(error)
            self.quibble = re.sub(
                r'([gcrnpmo])\.',
                r':\1.',
                self.quibble
            )


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
